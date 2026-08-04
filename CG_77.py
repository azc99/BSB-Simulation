import os
import subprocess
import sys

# =============================================================================
# GPU/WSL ENVIRONMENT
# =============================================================================

# Physical NVIDIA GPU selected before CUDA_VISIBLE_DEVICES is applied.
GPU_ID = int(os.environ.get("ION_GPU_ID", "0"))

# Number of independent configuration processes allowed to share the GPU.
# Start with 1. Benchmark 2 or more only after the single-worker version works.
MAX_GPU_WORKERS = int(os.environ.get("GPU_WORKERS", "7"))

# The host only launches GPU kernels. The expensive OpenMP TARGET loops run
# on the GPU, so each process normally needs only one host OpenMP thread.
HOST_OMP_THREADS = 1

if GPU_ID < 0:
    raise ValueError("ION_GPU_ID must be nonnegative.")
if MAX_GPU_WORKERS < 1:
    raise ValueError("GPU_WORKERS must be at least 1.")

_RESTART_FLAG = "ION_FITTER_GPU_RESTARTED"
_REQUIRED_ENV = {
    "CUDA_VISIBLE_DEVICES": str(GPU_ID),
    "OMP_TARGET_OFFLOAD": "MANDATORY",
    "OMP_NUM_THREADS": str(HOST_OMP_THREADS),
    "OMP_DYNAMIC": "FALSE",
    "OMP_DISPLAY_ENV": "FALSE",
    "OMP_DISPLAY_AFFINITY": "FALSE",
    "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "NUMEXPR_NUM_THREADS": "1",
}

needs_restart = any(
    os.environ.get(name) != value
    for name, value in _REQUIRED_ENV.items()
)

if needs_restart:
    if os.environ.get(_RESTART_FLAG) == "1":
        mismatches = {
            name: (os.environ.get(name), value)
            for name, value in _REQUIRED_ENV.items()
            if os.environ.get(name) != value
        }
        raise RuntimeError(
            "The fitter restarted, but the GPU environment was not preserved: "
            f"{mismatches}"
        )

    print(
        "Restarting Python with NVIDIA GPU offloading enabled "
        f"(physical GPU {GPU_ID}, host OMP threads={HOST_OMP_THREADS})...",
        flush=True,
    )

    child_env = os.environ.copy()
    child_env.update(_REQUIRED_ENV)
    child_env[_RESTART_FLAG] = "1"

    script_path = os.path.abspath(sys.argv[0])
    child_command = [sys.executable, script_path, *sys.argv[1:]]

    child = subprocess.Popen(child_command, env=child_env)

    try:
        child_returncode = child.wait()
    except KeyboardInterrupt:
        print("\nStopping restarted fitter...", flush=True)
        child.terminate()
        try:
            child.wait(timeout=3)
        except subprocess.TimeoutExpired:
            child.kill()
            child.wait()
        raise SystemExit(130)

    raise SystemExit(child_returncode)

from pathlib import Path

# =============================================================================
# CONFIGURATION
# =============================================================================

NUM_IONS = 7
NUM_MODES = 7
NUM_POINTS = 32

USE_PARALLEL = MAX_GPU_WORKERS > 1

DATA_FILE = Path("10ns_paper_77.txt")
SIM_LIBRARY = Path("libfortran_model.so")
DATA_LIBRARY = Path("libdata_output_model.so")
TARGET_OBJECTIVE = 3.0
CG_MAXITER = 5000
CG_GTOL = 1.0
CG_EPS = 1e-5

# True LD entries smaller than this are excluded from percentage-error calculations
PERCENT_REFERENCE_FLOOR = 0.015

import signal
from multiprocessing import (
    TimeoutError as MPTimeoutError,
    current_process,
    get_context,
)
from time import perf_counter

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize


# =============================================================================
# MULTIPROCESSING WORKER SETUP
# =============================================================================

_WORKER_STATE = {}


def _init_worker(
    j,
    k,
    mag_c,
    points,
    tau,
    w_access,
    pi_times,
    detune,
    worker_count,
):
    """Configure one spawned process before loading the GPU shared library."""
    # Only the main process handles Ctrl+C. A worker may be inside native
    # Fortran/CUDA code and should not process SIGINT independently.
    signal.signal(signal.SIGINT, signal.SIG_IGN)

    expected_environment = {
        "CUDA_VISIBLE_DEVICES": str(GPU_ID),
        "OMP_TARGET_OFFLOAD": "MANDATORY",
        "OMP_NUM_THREADS": str(HOST_OMP_THREADS),
    }

    for name, expected in expected_environment.items():
        actual = os.environ.get(name)
        if actual != expected:
            raise RuntimeError(
                f"GPU worker started with {name}={actual!r}; "
                f"expected {expected!r}."
            )

    identity = current_process()._identity
    worker_number = identity[0] if identity else 1
    worker_slot = (worker_number - 1) % worker_count

    print(
        f"GPU worker {worker_slot + 1}/{worker_count}: "
        f"CUDA_VISIBLE_DEVICES={os.environ['CUDA_VISIBLE_DEVICES']}, "
        f"OMP_NUM_THREADS={os.environ['OMP_NUM_THREADS']}",
        flush=True,
    )

    # Import only after the process has inherited the GPU/OpenMP environment.
    import fortran_module

    _WORKER_STATE.update(
        {
            "fortran_module": fortran_module,
            "j": int(j),
            "k": int(k),
            "mag_c": np.asfortranarray(mag_c, dtype=np.float64),
            "points": int(points),
            "tau": float(tau),
            "w_access": np.asfortranarray(w_access, dtype=np.float64),
            "pi_times": np.asfortranarray(pi_times, dtype=np.float64),
            "detune": np.asfortranarray(detune, dtype=np.float64),
        }
    )


def _run_single_config(task):
    """Run one configuration inside one persistent worker process."""
    config, ld = task
    state = _WORKER_STATE
    start = perf_counter()

    try:
        result = state["fortran_module"].run_model(
            state["j"],
            state["k"],
            ld,
            state["mag_c"],
            config,
            state["points"],
            state["tau"],
            state["w_access"],
            state["pi_times"],
            state["detune"],
        )
        elapsed = perf_counter() - start
        return config, float(result), elapsed
    except Exception as exc:
        elapsed = perf_counter() - start
        print(f"Error in config {config}: {exc}", flush=True)
        return config, 1e10, elapsed


# =============================================================================
# OPTIMIZER
# =============================================================================


class TargetReached(Exception):
    pass


class IonDynamicsOptimizer:
    def __init__(
        self,
        mag_c,
        num_ions,
        num_modes,
        num_points,
        tau,
        w_access,
        pi_times,
        detune,
        target_objective=TARGET_OBJECTIVE,
    ):
        self.num_ions = int(num_ions)
        self.num_modes = int(num_modes)
        self.num_points = int(num_points)
        self.tau = float(tau)

        self.mag_c = np.asfortranarray(mag_c, dtype=np.float64)
        self.w_access = np.asfortranarray(w_access, dtype=np.float64)
        self.pi_times = np.asfortranarray(pi_times, dtype=np.float64)
        self.detune = np.asfortranarray(detune, dtype=np.float64)

        self.n_free = self.num_ions * self.num_modes
        self.configs = list(range(1, self.num_ions + 1))

        self.n_calls = 0
        self.best_sum = None
        self.best_ld = None
        self.history = []
        self.objective_history = []
        self.pool = None
        self.target_objective = target_objective
        self.worker_count = min(MAX_GPU_WORKERS, self.num_ions)

    def _to_ld(self, params):
        params = np.asarray(params, dtype=np.float64)

        if params.size != self.n_free:
            raise ValueError(
                f"Expected {self.n_free} LD parameters, got {params.size}"
            )

        return np.asfortranarray(
            params.reshape(self.num_ions, self.num_modes),
            dtype=np.float64,
        )

    def _objective(self, params):
        call_number = self.n_calls + 1
        call_start = perf_counter()
        ld = self._to_ld(params)
        sums = np.zeros(self.num_ions, dtype=np.float64)
        config_times = np.zeros(self.num_ions, dtype=np.float64)

        if USE_PARALLEL:
            tasks = [(config, ld) for config in self.configs]

            # Poll map_async() so the main process remains responsive to
            # Ctrl+C while workers execute native Fortran/CUDA code.
            pending = self.pool.map_async(_run_single_config, tasks)
            while True:
                try:
                    config_results = pending.get(timeout=0.20)
                    break
                except MPTimeoutError:
                    continue

            for config, value, elapsed in config_results:
                sums[config - 1] = value
                config_times[config - 1] = elapsed
        else:
            # Lazy import so the shared library loads after the GPU/OpenMP settings.
            import fortran_module

            for config in self.configs:
                config_start = perf_counter()
                sums[config - 1] = fortran_module.run_model(
                    self.num_ions,
                    self.num_modes,
                    ld,
                    self.mag_c,
                    config,
                    self.num_points,
                    self.tau,
                    self.w_access,
                    self.pi_times,
                    self.detune,
                )
                config_times[config - 1] = perf_counter() - config_start

        result = float(np.sum(sums))
        call_elapsed = perf_counter() - call_start
        self.n_calls = call_number

        if not np.isfinite(result):
            print(
                f"Call {call_number}: non-finite objective; returning penalty",
                flush=True,
            )
            return 1e10

        if self.best_sum is None or result < self.best_sum:
            self.best_sum = result
            self.best_ld = ld.copy()

        self.history.append(
            {
                "i": call_number,
                "ld": ld.copy(),
                "sums": sums.copy(),
                "config_times": config_times.copy(),
                "total": result,
                "objective": result,
                "best": self.best_sum,
                "elapsed": call_elapsed,
            }
        )
        self.objective_history.append(self.best_sum)

        times_text = ", ".join(
            f"C{config}={config_times[config - 1]:.3f}s"
            for config in self.configs
        )

        print(
            f"Call {call_number}: objective={result:.12e}, "
            f"best={self.best_sum:.12e}, wall={call_elapsed:.3f}s "
            f"[{times_text}]",
            flush=True,
        )

        if (
            self.target_objective is not None
            and self.best_sum <= self.target_objective
        ):
            print(
                "Target objective reached. Stopping optimization.",
                flush=True,
            )
            raise TargetReached

        return result

    def optimize(self, initial_ld):
        initial_ld = np.asarray(initial_ld, dtype=np.float64)

        expected_shape = (self.num_ions, self.num_modes)
        if initial_ld.shape != expected_shape:
            raise ValueError(
                f"initial_ld must have shape {expected_shape}, "
                f"got {initial_ld.shape}"
            )

        initial_params = initial_ld.reshape(-1)
        optimization_start = perf_counter()

        if USE_PARALLEL:
            # CUDA runtimes should not be inherited through fork. "spawn"
            # starts each GPU worker in a fresh Python process.
            context = get_context("spawn")
            self.pool = context.Pool(
                processes=self.worker_count,
                initializer=_init_worker,
                initargs=(
                    self.num_ions,
                    self.num_modes,
                    self.mag_c,
                    self.num_points,
                    self.tau,
                    self.w_access,
                    self.pi_times,
                    self.detune,
                    self.worker_count,
                ),
            )

        try:
            result = minimize(
                self._objective,
                initial_params,
                method="CG",
                callback=None,
                options={
                    "disp": True,
                    "maxiter": CG_MAXITER,
                    "gtol": CG_GTOL,
                    "norm": 2,
                    "eps": CG_EPS,
                },
            )
        except TargetReached:
            print(
                "Optimization terminated early because the target was reached."
            )

            class EarlyResult:
                pass

            result = EarlyResult()
            result.x = self.best_ld.reshape(-1)
            result.fun = self.best_sum
            result.nit = len(self.objective_history)
            result.nfev = self.n_calls
            result.success = True
            result.message = "Stopped early: target objective reached."
        except KeyboardInterrupt:
            print(
                "\nCtrl+C received. Terminating Fortran worker processes...",
                flush=True,
            )
            if self.pool is not None:
                self.pool.terminate()
                self.pool.join()
                self.pool = None
            raise
        finally:
            if self.pool is not None:
                # Normal completion: finish queued work and shut down cleanly.
                self.pool.close()
                self.pool.join()
                self.pool = None

        result.elapsed = perf_counter() - optimization_start
        result.ld_opt = self._to_ld(result.x)
        return result

    def print_results(self, result):
        print(
            f"Iterations: {result.nit}, "
            f"runtime: {result.elapsed:.1f}s, "
            f"objective: {result.fun:.6e}"
        )
        print(f"\nOptimal LD matrix:\n{result.ld_opt}")

        if self.history:
            best = min(self.history, key=lambda item: item["objective"])
            print(
                f"\nBest objective: {best['objective']:.6e}\n"
                f"Configuration sums: {best['sums']}\n"
                f"Configuration runtimes: {best['config_times']}"
            )

    def compare_with_true_ld(
        self,
        true_ld,
        output_dir=".",
        reference_floor=PERCENT_REFERENCE_FLOOR,
    ):
        """Compare the best fitted LD matrix with the known true LD matrix."""
        if self.best_ld is None:
            raise RuntimeError("No best LD matrix is available for comparison.")

        true_ld = np.asarray(true_ld, dtype=np.float64)
        expected_shape = (self.num_ions, self.num_modes)

        if true_ld.shape != expected_shape:
            raise ValueError(
                f"true_ld must have shape {expected_shape}, "
                f"got {true_ld.shape}"
            )

        if reference_floor < 0:
            raise ValueError("reference_floor must be nonnegative.")

        valid = np.abs(true_ld) >= reference_floor
        percentage_difference = np.full(
            expected_shape,
            np.nan,
            dtype=np.float64,
        )

        percentage_difference[valid] = (
            np.abs(self.best_ld[valid] - true_ld[valid])
            / np.abs(true_ld[valid])
            * 100.0
        )

        if not np.any(valid):
            raise ValueError(
                "No true LD entries are large enough for percentage comparison."
            )

        average_percentage_difference = float(
            np.mean(percentage_difference[valid])
        )

        print("\nTrue LD matrix:")
        print(true_ld)
        print("\nBest fitted LD matrix:")
        print(self.best_ld)
        print("\nAbsolute percentage-difference matrix (%):")
        print(percentage_difference)
        print(
            f"Average percentage difference: "
            f"{average_percentage_difference:.6f}%"
        )

        excluded_count = int(np.size(valid) - np.count_nonzero(valid))
        if excluded_count:
            print(
                f"Excluded {excluded_count} near-zero true LD "
                f"{'entry' if excluded_count == 1 else 'entries'} "
                f"with |true_ld| < {reference_floor:.1e}; "
                "these appear as nan in the matrix."
            )

        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        comparison_file = output_path / "CG_ld_percent_difference.txt"

        np.savetxt(
            comparison_file,
            percentage_difference,
            fmt="%.10e",
            header=(
                "Absolute percentage difference between best_ld and true_ld. "
                f"nan indicates |true_ld| < {reference_floor:.1e}. "
                f"Average over included entries = "
                f"{average_percentage_difference:.10e}%"
            ),
        )

        print(f"Percentage differences saved to: {comparison_file.resolve()}")
        return percentage_difference, average_percentage_difference

    def plot_convergence(self, filename="CG_convergence.png"):
        iterations = np.arange(1, len(self.objective_history) + 1)

        plt.figure(figsize=(8, 5))
        plt.plot(
            iterations,
            self.objective_history,
            linewidth=1.5,
            label="best after each objective call",
        )
        plt.legend()
        plt.xlabel("Objective call")
        plt.ylabel("Best objective")
        plt.title("Convergence")
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(filename, dpi=150)
        plt.close()

    def _generate_combined_output(self, ld, combined_path):
        # This import occurs in the main process after the environment setup.
        import data_output_module

        output_files = [
            Path(f"config{config}_total.txt")
            for config in range(self.num_ions)
        ]

        try:
            for config in self.configs:
                data_output_module.run_model(
                    self.num_ions,
                    self.num_modes,
                    ld,
                    config,
                    self.num_points,
                    self.tau,
                    self.w_access,
                    self.pi_times,
                    self.detune,
                )

            arrays = [np.loadtxt(path) for path in output_files]
            combined = arrays[0]

            for array in arrays[1:]:
                combined = np.hstack((combined, array))

            np.savetxt(combined_path, combined)
            return combined
        finally:
            for path in output_files:
                try:
                    path.unlink()
                except FileNotFoundError:
                    pass

    def plot_results(self, filename="CG_after_77.png", ld=None):
        if ld is None:
            if self.best_ld is None:
                raise RuntimeError("No optimized LD matrix is available.")
            plot_ld = self.best_ld
        else:
            plot_ld = np.asfortranarray(ld, dtype=np.float64)
            expected_shape = (self.num_ions, self.num_modes)

            if plot_ld.shape != expected_shape:
                raise ValueError(
                    f"LD matrix must have shape {expected_shape}, "
                    f"got {plot_ld.shape}"
                )

        combined_path = Path("config_total.txt")

        try:
            simulated = self._generate_combined_output(
                plot_ld,
                combined_path,
            )

            experimental = self.mag_c
            time_experimental = experimental[:, 0]
            time_simulated = simulated[:, 0]

            output_path = Path(filename).resolve()
            output_path.parent.mkdir(parents=True, exist_ok=True)

            fig, axes = plt.subplots(
                nrows=self.num_ions,
                ncols=self.num_modes,
                figsize=(24, 20),
                sharex=True,
                sharey=True,
                constrained_layout=True,
            )

            axes = np.asarray(axes)

            if axes.ndim == 1:
                if self.num_ions == 1:
                    axes = axes[np.newaxis, :]
                else:
                    axes = axes[:, np.newaxis]

            for config_index in range(self.num_ions):
                for ion_index in range(self.num_modes):
                    axis = axes[config_index, ion_index]

                    data_column = (
                            self.num_ions * config_index
                            + ion_index
                            + 1
                    )

                    mode = (
                                   self.num_modes
                                   - config_index
                                   + ion_index
                           ) % self.num_modes + 1

                    axis.plot(
                        time_experimental,
                        experimental[:, data_column],
                        linewidth=1.2,
                        label="Experimental",
                    )

                    axis.scatter(
                        time_experimental,
                        experimental[:, data_column],
                        s=14,
                        marker="x",
                    )

                    axis.plot(
                        time_simulated,
                        simulated[:, data_column],
                        linewidth=1.2,
                        label="Simulated",
                    )

                    axis.scatter(
                        time_simulated,
                        simulated[:, data_column],
                        s=10,
                        marker="x",
                    )

                    # Shorter title reduces clutter.
                    axis.set_title(
                        f"Ion {ion_index + 1}, M{mode}",
                        fontsize=10,
                        pad=3,
                    )

                    axis.set_ylim(0.0, 1.0)
                    axis.grid(True, alpha=0.2)
                    axis.tick_params(
                        axis="both",
                        labelsize=8,
                    )

                # One configuration label per row.
                axes[config_index, 0].annotate(
                    f"Config {config_index + 1}",
                    xy=(-0.32, 0.5),
                    xycoords="axes fraction",
                    rotation=90,
                    ha="center",
                    va="center",
                    fontsize=12,
                    fontweight="bold",
                )

            handles, labels = axes[0, 0].get_legend_handles_labels()

            fig.legend(
                handles,
                labels,
                loc="upper center",
                ncol=2,
                bbox_to_anchor=(0.5, 1.01),
                fontsize=12,
            )

            fig.suptitle(
                output_path.stem.replace("_", " "),
                fontsize=17,
            )

            fig.supxlabel(
                "Duration (µs)",
                fontsize=14,
            )

            fig.supylabel(
                "Bright-state probability",
                fontsize=14,
            )

            fig.savefig(
                output_path,
                dpi=250,
                bbox_inches="tight",
            )

            plt.close(fig)

            print(
                f"Plot saved to: {output_path}",
                flush=True,
            )

        finally:
            try:
                combined_path.unlink()
            except FileNotFoundError:
                pass

    def save_results(self, result, output_dir="."):
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        np.savetxt(
            output_path / "CG_ld.txt",
            result.ld_opt,
            fmt="%.10e",
            header=(
                f"obj={result.fun:.6e}, "
                f"t={result.elapsed:.1f}s, "
                f"function evals={result.nfev}"
            ),
        )

        if self.history:
            np.savetxt(
                output_path / "history.txt",
                [
                    [
                        item["i"],
                        item["objective"],
                        item["best"],
                        item["elapsed"],
                    ]
                    for item in self.history
                ],
                fmt="%.10e",
                header="call objective best wall_seconds",
            )

        print(f"Results saved to {output_path.resolve()}")


# =============================================================================
# MAIN
# =============================================================================


def main():
    if os.name == "nt":
        sys.exit(
            "Error: run this fitter with Linux Python inside WSL, not "
            "Windows Python."
        )

    print(
        f"GPU configuration: physical GPU {GPU_ID}, "
        f"{MAX_GPU_WORKERS} configuration worker"
        f"{'s' if MAX_GPU_WORKERS != 1 else ''}, "
        f"{HOST_OMP_THREADS} host OpenMP thread per process"
    )
    print(
        "Startup environment: "
        f"CUDA_VISIBLE_DEVICES={os.environ.get('CUDA_VISIBLE_DEVICES')}, "
        f"OMP_TARGET_OFFLOAD={os.environ.get('OMP_TARGET_OFFLOAD')}, "
        f"OMP_NUM_THREADS={os.environ.get('OMP_NUM_THREADS')}"
    )

    if MAX_GPU_WORKERS > 1:
        print(
            "Warning: multiple processes will share one GPU. Benchmark against "
            "GPU_WORKERS=1 because more workers may be slower.",
            flush=True,
        )

    try:
        gpu_result = subprocess.run(
            [
                "nvidia-smi",
                "--id",
                str(GPU_ID),
                "--query-gpu=name,memory.total,driver_version",
                "--format=csv,noheader",
            ],
            check=True,
            capture_output=True,
            text=True,
        )
    except FileNotFoundError:
        sys.exit("Error: nvidia-smi was not found inside WSL.")
    except subprocess.CalledProcessError as exc:
        message = exc.stderr.strip() or exc.stdout.strip()
        sys.exit(f"Error: GPU {GPU_ID} could not be queried: {message}")

    print(f"GPU: {gpu_result.stdout.strip()}")

    for library in (SIM_LIBRARY, DATA_LIBRARY):
        if not library.exists():
            sys.exit(
                f"Error: {library} was not found. Compile the nvfortran "
                "shared libraries before running the fitter."
            )

    if not DATA_FILE.exists():
        sys.exit(f"Error: {DATA_FILE} was not found.")

    data = np.loadtxt(DATA_FILE, dtype=np.float64)
    expected_shape = (
        NUM_POINTS + 1,
        1 + NUM_IONS * NUM_MODES,
    )

    if data.shape != expected_shape:
        sys.exit(
            f"Error: data shape {data.shape}, expected {expected_shape}"
        )

    tau_us = round(data[-1, 0]) + .0000001
    tau = tau_us * 1e-6
    print(f"tau = {tau_us:.9f} us = {tau:.12e} s")

    initial_ld = np.array([
      [ 0.04641081,  0.06193354, -0.06183822, -0.04480131, -0.02324654, -0.01237491, -0.00358810],
      [ 0.04302856,  0.03665466,  0.00363358,  0.04539584,  0.06761749,  0.04917535,  0.02299433],
      [ 0.03845976,  0.01824714,  0.03790939,  0.04043542, -0.01065526, -0.06281157, -0.05804929],
      [ 0.03718567,  1.20355e-5, 0.05242114,  5.70065e-5, -0.05051431, -9.03943e-5, 0.06706654],
      [ 0.03621037, -0.02000657,  0.03944061, -0.04844214, -0.00984468,  0.05961443, -0.05476451],
      [ 0.03817230, -0.03666898,  0.00371188, -0.04414933,  0.06441588, -0.04710043,  0.02303717],
      [ 0.04576233, -0.06034708, -0.06153428,  0.04562097, -0.02463289,  0.01082162, -0.00351373]
    ])

    # initial_ld = np.array([
      #  [0.0440,  0.0627, -0.0577, -0.0431, -0.0253, -0.0113, -0.00341],
      #  [0.0407,  0.0396,  0.00367, 0.0466,  0.0623,  0.0478,  0.0216 ],
      #  [0.0389,  0.0193,  0.0375,  0.0443, -0.0100, -0.0612, -0.0552 ],
      #  [0.0383,  1.10e-5, 0.0486,  5.40e-5, -0.0538, -9.10e-5, 0.0738 ],
      #  [0.0389, -0.0193,  0.0376, -0.0443, -0.0102,  0.0612, -0.0551 ],
      #  [0.0407, -0.0396,  0.00373, -0.0467, 0.0623, -0.0477,  0.0216 ],
      #  [0.0440, -0.0627, -0.0577,  0.0430, -0.0252,  0.0113, -0.00339]
    #])
    # Known true Lamb-Dicke parameters used to generate the reference data.
    true_ld = np.array(
        [
          [0.0440,  0.0627, -0.0577, -0.0431, -0.0253, -0.0113, -0.00341],
          [0.0407,  0.0396,  0.00367, 0.0466,  0.0623,  0.0478,  0.0216 ],
          [0.0389,  0.0193,  0.0375,  0.0443, -0.0100, -0.0612, -0.0552 ],
          [0.0383,  1.10e-5, 0.0486,  5.40e-5, -0.0538, -9.10e-5, 0.0738 ],
          [0.0389, -0.0193,  0.0376, -0.0443, -0.0102,  0.0612, -0.0551 ],
          [0.0407, -0.0396,  0.00373, -0.0467, 0.0623, -0.0477,  0.0216 ],
          [0.0440, -0.0627, -0.0577,  0.0430, -0.0252,  0.0113, -0.00339]
        ],
        dtype=np.float64,
    )

    w_access = np.zeros(7, dtype=np.float64)
    w_access[:7] = (
        np.array(
            [
                205.351994,
                205.367894,
                205.397094,
                205.431194,
                205.469594,
                205.511094,
                205.554694
            ],
            dtype=np.float64,
        )
        * 2.0
        * np.pi
    )

    pi_times = np.zeros(7, dtype=np.float64)
    pi_times[:7] = [11.32, 9.38, 9.59, 9.88, 8.77, 9.34,9.41]

    # 1e-5 matches the Fortran value 1.d-5.
    detune = np.full((7, 7), 1e-5, dtype=np.float64)

    optimizer = IonDynamicsOptimizer(
        data,
        NUM_IONS,
        NUM_MODES,
        NUM_POINTS,
        tau,
        w_access,
        pi_times,
        detune,
    )

    script_dir = Path(__file__).resolve().parent
    before_plot = script_dir / "CG_before_77.png"
    after_plot = script_dir / "CG_after_77.png"

    print(f"\nGenerating pre-fit plot: {before_plot}")
    optimizer.plot_results(filename=before_plot, ld=initial_ld)

    result = optimizer.optimize(initial_ld)

    print(result.message)
    optimizer.print_results(result)
    optimizer.plot_convergence(filename=script_dir / "CG_convergence.png")
    optimizer.plot_results(filename=after_plot)
    optimizer.save_results(result, output_dir=script_dir)
    optimizer.compare_with_true_ld(true_ld, output_dir=script_dir)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print("Interrupted by user.", flush=True)
        raise SystemExit(130)
