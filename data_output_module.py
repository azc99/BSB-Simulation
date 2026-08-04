import ctypes
import os
from pathlib import Path

import numpy as np
from numpy.ctypeslib import ndpointer


# This wrapper is for Linux/WSL.
if os.name == "nt":
    raise RuntimeError(
        "data_output_module.py must be run using Linux Python inside WSL. "
        "The GPU library is a Linux .so file, not a Windows .dll."
    )


_module_dir = Path(__file__).resolve().parent
_lib_path = _module_dir / "libdata_output_model.so"

if not _lib_path.exists():
    raise FileNotFoundError(
        f"Fortran shared library not found: {_lib_path}\n"
        "Compile output_66.f with nvfortran before importing this module."
    )


try:
    _lib = ctypes.CDLL(
        str(_lib_path),
        mode=ctypes.RTLD_LOCAL,
    )
except OSError as exc:
    raise OSError(
        f"Could not load {_lib_path}\n"
        f"Original error: {exc}\n\n"
        "Check missing dependencies with:\n"
        f"    ldd {_lib_path}"
    ) from exc


try:
    _run_model = _lib.run_model_
except AttributeError as exc:
    raise AttributeError(
        f"The symbol run_model_ was not found in {_lib_path}.\n"
        "Check exported symbols with:\n"
        f"    nm -D {_lib_path} | grep -i run_model"
    ) from exc


_run_model.argtypes = [
    ctypes.POINTER(ctypes.c_int),                              # j
    ctypes.POINTER(ctypes.c_int),                              # k

    ndpointer(
        dtype=np.float64,
        ndim=2,
        flags="F_CONTIGUOUS",
    ),                                                         # LD(j,k)

    ctypes.POINTER(ctypes.c_int),                              # config
    ctypes.POINTER(ctypes.c_int),                              # points
    ctypes.POINTER(ctypes.c_double),                           # tau

    ndpointer(
        dtype=np.float64,
        ndim=1,
        flags="F_CONTIGUOUS",
    ),                                                         # w_access

    ndpointer(
        dtype=np.float64,
        ndim=1,
        flags="F_CONTIGUOUS",
    ),                                                         # pi_times

    ndpointer(
        dtype=np.float64,
        ndim=2,
        flags="F_CONTIGUOUS",
    ),                                                         # detune
]

_run_model.restype = None


def run_model(
    j,
    k,
    ld,
    config,
    points,
    tau,
    w_access,
    pi_times,
    detune,
):
    """Generate output data for one GPU-accelerated configuration."""

    j = int(j)
    k = int(k)
    config = int(config)
    points = int(points)
    tau = float(tau)

    j_c = ctypes.c_int(j)
    k_c = ctypes.c_int(k)
    config_c = ctypes.c_int(config)
    points_c = ctypes.c_int(points)
    tau_c = ctypes.c_double(tau)

    ld_f = np.asfortranarray(ld, dtype=np.float64)
    w_access_f = np.asfortranarray(w_access, dtype=np.float64)
    pi_times_f = np.asfortranarray(pi_times, dtype=np.float64)
    detune_f = np.asfortranarray(detune, dtype=np.float64)

    if ld_f.shape != (j, k):
        raise ValueError(
            f"LD must have shape ({j}, {k}), got {ld_f.shape}"
        )

    if w_access_f.shape != (7,):
        raise ValueError(
            f"w_access must have shape (7,), got {w_access_f.shape}"
        )

    if pi_times_f.shape != (7,):
        raise ValueError(
            f"pi_times must have shape (7,), got {pi_times_f.shape}"
        )

    if detune_f.shape != (7, 7):
        raise ValueError(
            f"detune must have shape (7, 7), got {detune_f.shape}"
        )

    _run_model(
        ctypes.byref(j_c),
        ctypes.byref(k_c),
        ld_f,
        ctypes.byref(config_c),
        ctypes.byref(points_c),
        ctypes.byref(tau_c),
        w_access_f,
        pi_times_f,
        detune_f,
    )
