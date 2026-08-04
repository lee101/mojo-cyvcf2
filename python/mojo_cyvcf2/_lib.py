"""ctypes loader for the Mojo VCF scanner."""
from __future__ import annotations

import ctypes
import os
import shutil
import subprocess

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
LIB = os.path.join(ROOT, "dist", "libmojo-cyvcf2.so")
I = ctypes.c_int64

_SIGNATURES = {
    "mcv_scan_records": ([I, I, I, I, I], I),
    "mcv_decode_gt": ([I] * 8, None),
}
_loaded = None


def build(force: bool = False) -> str:
    sources = [os.path.join(ROOT, "src", "capi.mojo")]
    if not force and os.path.exists(LIB) and os.path.getmtime(LIB) >= max(map(os.path.getmtime, sources)):
        return LIB
    mojo = shutil.which("mojo")
    if mojo is None:
        pixi = shutil.which("pixi") or os.path.expanduser("~/.pixi/bin/pixi")
        command = [pixi, "run", "--manifest-path", os.path.join(ROOT, "pixi.toml"), "bash", "build/build.sh"]
    else:
        command = ["bash", os.path.join(ROOT, "build", "build.sh")]
    proc = subprocess.run(command, cwd=ROOT, capture_output=True, text=True, timeout=1800)
    if proc.returncode or not os.path.exists(LIB):
        raise RuntimeError((proc.stderr or proc.stdout).strip())
    return LIB


def lib() -> ctypes.CDLL:
    global _loaded
    if _loaded is None:
        _loaded = ctypes.CDLL(build())
        for name, (argtypes, restype) in _SIGNATURES.items():
            fn = getattr(_loaded, name)
            fn.argtypes, fn.restype = argtypes, restype
    return _loaded


def addr(a: np.ndarray) -> int:
    return int(a.ctypes.data)
