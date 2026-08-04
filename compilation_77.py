import os
import subprocess
import sys
from pathlib import Path

GPU_ID = 0

# The large OpenMP TARGET loops execute on the GPU.
# The host generally only needs one CPU thread to launch kernels.
HOST_OMP_THREADS = 1

required_environment = {
    "OMP_NUM_THREADS": str(HOST_OMP_THREADS),
    "OMP_DYNAMIC": "FALSE",

    # Do not silently run target regions on the CPU.
    "OMP_TARGET_OFFLOAD": "MANDATORY",

    # Make only this GPU visible to the process.
    "CUDA_VISIBLE_DEVICES": str(GPU_ID),
}

# CUDA_VISIBLE_DEVICES and OpenMP variables should be present before
# the Fortran/OpenMP/CUDA runtime is initialized.
needs_restart = any(
    os.environ.get(name) != value
    for name, value in required_environment.items()
)

if needs_restart:
    print(
        "Restarting Python with GPU/OpenMP environment variables...",
        flush=True,
    )

    child_env = os.environ.copy()
    child_env.update(required_environment)

    completed = subprocess.run(
        [sys.executable, *sys.argv],
        env=child_env,
    )
    raise SystemExit(completed.returncode)

# Optional OpenMP diagnostics
os.environ["OMP_DISPLAY_ENV"] = "FALSE"

# Set environment variables before importing NumPy or the wrappers.
import importlib
from time import perf_counter

import numpy as np


# ------------------------------------------------------------
# Files
# ------------------------------------------------------------

SIM_SOURCE = Path("sim_77.f")
SIM_LIBRARY = Path("libfortran_model.so")

DATA_SOURCE = Path("output_77.f")
DATA_LIBRARY = Path("libdata_output_model.so")

DATA_FILE = Path("10ns_paper_77.txt")


# ------------------------------------------------------------
# GPU/compiler checks
# ------------------------------------------------------------

def check_environment() -> bool:
    """Check that the script is running in Linux/WSL with an NVIDIA GPU."""

    if os.name == "nt":
        print(
            "ERROR: This GPU build must be run with Linux Python inside WSL.\n"
            "Open a WSL terminal and run:\n"
            "    python3 build_gpu_modules.py"
        )
        return False

    try:
        version = subprocess.run(
            ["nvfortran", "-V"],
            check=True,
            capture_output=True,
            text=True,
        )
    except FileNotFoundError:
        print("ERROR: nvfortran was not found on PATH")
        return False
    except subprocess.CalledProcessError as exc:
        print("ERROR: nvfortran -V failed")
        if exc.stdout:
            print(exc.stdout)
        if exc.stderr:
            print(exc.stderr)
        return False

    compiler_output = version.stdout or version.stderr
    print("\nNVIDIA Fortran compiler:")
    print(compiler_output.strip())

    try:
        gpu = subprocess.run(
            [
                "nvidia-smi",
                "--query-gpu=index,name,compute_cap,memory.total",
                "--format=csv,noheader",
            ],
            check=True,
            capture_output=True,
            text=True,
        )
    except FileNotFoundError:
        print("ERROR: nvidia-smi was not found")
        return False
    except subprocess.CalledProcessError as exc:
        print("ERROR: nvidia-smi failed")
        if exc.stdout:
            print(exc.stdout)
        if exc.stderr:
            print(exc.stderr)
        return False

    gpu_lines = [
        line.strip()
        for line in gpu.stdout.splitlines()
        if line.strip()
    ]

    if not gpu_lines:
        print("ERROR: No NVIDIA GPU was detected")
        return False

    print("\nDetected NVIDIA GPU(s):")
    for line in gpu_lines:
        print(f"  {line}")

    if GPU_ID >= len(gpu_lines):
        print(
            f"ERROR: GPU_ID={GPU_ID}, but only "
            f"{len(gpu_lines)} GPU(s) were detected"
        )
        return False

    print(f"\nCUDA_VISIBLE_DEVICES={GPU_ID}")
    print(f"OMP_TARGET_OFFLOAD={os.environ['OMP_TARGET_OFFLOAD']}")
    print(f"OMP_NUM_THREADS={os.environ['OMP_NUM_THREADS']}")

    return True


# ------------------------------------------------------------
# Compilation
# ------------------------------------------------------------

def compile_shared_library(source: Path, output: Path) -> bool:
    """
    Compile one fixed-form Fortran file into a Linux shared library
    containing OpenMP GPU target code.
    """

    if not source.exists():
        print(f"ERROR: {source} not found")
        return False

    # Remove an old library so a failed compilation cannot leave us
    # accidentally testing a stale build.
    try:
        output.unlink()
    except FileNotFoundError:
        pass

    cmd = [
        "nvfortran",

        # Create a Linux shared object.
        "-shared",
        "-fpic",

        # Bind calls such as specific_time_ and spin_only_ to the
        # definitions inside this library. Both GPU libraries contain
        # routines with the same Fortran symbol names.
        "-Wl,-Bsymbolic",

        # Optimization.
        "-O3",

        # Compile OpenMP TARGET regions for NVIDIA GPUs.
        "-mp=gpu",

        # RTX 4070 Ti, compute capability 8.9.
        "-gpu=cc89",

        # Allow extended fixed-form source lines.
        "-Mextend",

        # Report which OpenMP regions generated GPU code.
        "-Minfo=mp",

        str(source),
        "-o",
        str(output),
    ]

    print(f"\nCompiling {source} -> {output}")
    print(" ".join(cmd))

    try:
        result = subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True,
        )
    except FileNotFoundError:
        print("ERROR: nvfortran was not found on PATH")
        return False
    except subprocess.CalledProcessError as exc:
        print("\nCompilation failed")

        if exc.stdout:
            print("\nCompiler stdout:")
            print(exc.stdout)

        if exc.stderr:
            print("\nCompiler stderr:")
            print(exc.stderr)

        return False

    # -Minfo=mp normally prints useful information to stderr even when
    # compilation succeeds.
    if result.stdout.strip():
        print("\nCompiler stdout:")
        print(result.stdout)

    if result.stderr.strip():
        print("\nCompiler GPU report:")
        print(result.stderr)

    if not output.exists():
        print(f"ERROR: compiler returned success but {output} was not created")
        return False

    print(f"Created {output.resolve()}")
    return True


# ------------------------------------------------------------
# Test inputs
# ------------------------------------------------------------

def make_test_inputs():
    """Create inputs for the hard-coded 7-ion/7-mode implementation."""

    j = 7
    k = 7
    points = 32
    tau = 800e-6

    w_access = np.zeros(7, dtype=np.float64)
    w_access[:7] = np.abs(
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
    pi_times[:7] = [
        11.32,
        9.38,
        9.59,
        9.88,
        8.77,
        9.34,
        9.41,
    ]

    detune = np.full((7, 7), 10e-6, dtype=np.float64)
    
    ld_test = np.array([
      [0.0440,  0.0627, -0.0577, -0.0431, -0.0253, -0.0113, -0.00341],
      [0.0407,  0.0396,  0.00367, 0.0466,  0.0623,  0.0478,  0.0216 ],
      [0.0389,  0.0193,  0.0375,  0.0443, -0.0100, -0.0612, -0.0552 ],
      [0.0383,  1.10e-5, 0.0486,  5.40e-5, -0.0538, -9.10e-5, 0.0738 ],
      [0.0389, -0.0193,  0.0376, -0.0443, -0.0102,  0.0612, -0.0551 ],
      [0.0407, -0.0396,  0.00373, -0.0467, 0.0623, -0.0477,  0.0216 ],
      [0.0440, -0.0627, -0.0577,  0.0430, -0.0252,  0.0113, -0.00339]
    ])

    return (
        j,
        k,
        points,
        tau,
        w_access,
        pi_times,
        detune,
        ld_test,
    )


# ------------------------------------------------------------
# Import helper
# ------------------------------------------------------------

def import_wrapper(module_name: str, library: Path):
    """
    Import one of the existing Python wrappers.

    The wrapper must load the Linux .so library rather than the old
    Windows .dll library.
    """

    sys.modules.pop(module_name, None)
    importlib.invalidate_caches()

    try:
        return importlib.import_module(module_name)
    except (ImportError, OSError) as exc:
        print(f"\nERROR: could not import {module_name}")
        print(exc)
        print(
            "\nCheck that the Python wrapper loads:\n"
            f"    {library.resolve()}\n"
            "rather than a .dll file."
        )
        raise


# ------------------------------------------------------------
# Simulation test
# ------------------------------------------------------------

def test_module() -> bool:
    if not DATA_FILE.exists():
        print(f"ERROR: {DATA_FILE} not found")
        return False

    try:
        fortran_module = import_wrapper(
            "fortran_module",
            SIM_LIBRARY,
        )
    except (ImportError, OSError):
        return False

    print(
        "\nSuccessfully imported fortran_module\n"
        f"GPU ID: {GPU_ID}\n"
        f"Host OpenMP threads: {HOST_OMP_THREADS}"
    )

    (
        j,
        k,
        points,
        tau,
        w_access,
        pi_times,
        detune,
        ld_test,
    ) = make_test_inputs()

    mag_c_test = np.loadtxt(DATA_FILE, dtype=np.float64)

    expected_shape = (points + 1, j * k + 1)

    if mag_c_test.shape != expected_shape:
        print(
            f"ERROR: {DATA_FILE} must have shape {expected_shape}, "
            f"got {mag_c_test.shape}"
        )
        return False

    print(
        f"System: {j} ions x {k} modes, "
        f"{points + 1} time points, tau={tau:.2e} s"
    )

    total = 0.0
    total_start = perf_counter()

    # These calls are sequential. Each call launches GPU kernels internally.
    for config in range(1, j + 1):
        config_start = perf_counter()

        result = fortran_module.run_model(
            j,
            k,
            ld_test,
            mag_c_test,
            config,
            points,
            tau,
            w_access,
            pi_times,
            detune,
        )

        config_elapsed = perf_counter() - config_start
        total += result

        print(
            f"  Config {config}: {result:.12e} "
            f"| runtime: {config_elapsed:.3f} s"
        )

    total_elapsed = perf_counter() - total_start

    if not np.isfinite(total):
        print(f"ERROR: non-finite objective total: {total}")
        return False

    print(f"Total: {total:.12e}")
    print(f"Total simulation runtime: {total_elapsed:.3f} s")

    return True


# ------------------------------------------------------------
# Data-output test
# ------------------------------------------------------------

def test_data_output() -> bool:
    try:
        data_output_module = import_wrapper(
            "data_output_module",
            DATA_LIBRARY,
        )
    except (ImportError, OSError):
        return False

    print(
        "\nSuccessfully imported data_output_module\n"
        f"GPU ID: {GPU_ID}\n"
        f"Host OpenMP threads: {HOST_OMP_THREADS}"
    )

    (
        j,
        k,
        points,
        tau,
        w_access,
        pi_times,
        detune,
        ld_test,
    ) = make_test_inputs()

    output_files = [
        Path(f"config{config}_total.txt")
        for config in range(j)
    ]

    try:
        for config in range(1, j + 1):
            config_start = perf_counter()

            data_output_module.run_model(
                j,
                k,
                ld_test,
                config,
                points,
                tau,
                w_access,
                pi_times,
                detune,
            )

            elapsed = perf_counter() - config_start

            print(
                f"  Data config {config}: "
                f"runtime {elapsed:.3f} s"
            )

        for path in output_files:
            if not path.exists():
                print(
                    "ERROR: expected output file was not created: "
                    f"{path}"
                )
                return False

            data = np.loadtxt(path)

            if data.shape[0] != points + 1:
                print(
                    f"ERROR: {path} should contain "
                    f"{points + 1} rows, got {data.shape[0]}"
                )
                return False

        print("Data-output GPU smoke test passed")
        return True

    finally:
        for path in output_files:
            try:
                path.unlink()
            except FileNotFoundError:
                pass


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def main() -> int:
    print("GPU Fortran build and smoke test")
    print(f"Python executable: {sys.executable}")

    if not check_environment():
        return 1

    if not compile_shared_library(SIM_SOURCE, SIM_LIBRARY):
        return 1

    if not compile_shared_library(DATA_SOURCE, DATA_LIBRARY):
        return 1

    print("\n--- Testing simulation shared library ---")

    if not test_module():
        return 1

    print("\n--- Testing data-output shared library ---")

    if not test_data_output():
        return 1

    print("\nDone.")
    print("Run the optimizer inside WSL with:")
    print("    python3 CG_77.py")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
