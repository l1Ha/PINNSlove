"""Robust loader for the native BKMP2 (H3) potential-energy surface.

The BKMP2 Fortran source lives in ``peslib/Amod_bkmp2.f`` with a small C-binding
wrapper ``peslib/bkmp2_cwrap.f90``.  This module:

* locates a prebuilt ``peslib/libbkmp2.so`` (any platform suffix), or
* **automatically compiles it** with gfortran when a compiler is available, or
* raises a single clear, actionable error otherwise.

It also provides a vectorised grid evaluator with an on-disk cache so that the
expensive per-point ctypes loop is paid only once per grid.

Units: the Fortran routine takes internuclear distances in **Bohr** and returns
the energy in **Hartree**.  Helpers here accept Angstrom and return cm^-1.
"""
from pathlib import Path
import ctypes
import hashlib
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
PESDIR = ROOT / 'peslib'

BOHR = 0.529177210903
HARTREE_TO_CM1 = 219474.6313705

# Fortran source files (order matters: module first, wrapper second)
_SRC_FILES = ['Amod_bkmp2.f', 'bkmp2_cwrap.f90']

_LIB = None            # cached ctypes handle


class BKMP2Unavailable(RuntimeError):
    """Raised when the native BKMP2 library cannot be located or built."""


def _candidate_names():
    names = ['libbkmp2.so', 'libbkmp2.dylib', 'bkmp2.dll']
    return [PESDIR / n for n in names]


def find_library():
    """Return the path of an existing shared library, or None."""
    for p in _candidate_names():
        if p.exists():
            return p
    # also accept any libbkmp2* the user may have dropped in the folder
    hits = sorted(PESDIR.glob('libbkmp2*')) if PESDIR.exists() else []
    return hits[0] if hits else None


def build_library(verbose=True):
    """Compile the shared library with gfortran. Returns its path.

    Raises BKMP2Unavailable with build instructions on failure.
    """
    if not PESDIR.exists():
        raise BKMP2Unavailable(
            f"Potential-energy-surface sources not found: {PESDIR} is missing.\n"
            "Re-extract the accompanying package so that peslib/Amod_bkmp2.f "
            "and peslib/bkmp2_cwrap.f90 are present."
        )
    missing = [f for f in _SRC_FILES if not (PESDIR / f).exists()]
    if missing:
        raise BKMP2Unavailable(
            "BKMP2 Fortran sources are incomplete; missing: "
            + ", ".join(str(PESDIR / m) for m in missing)
        )

    out = _candidate_names()[0]
    cmd = ['gfortran', '-O2', '-shared', '-fPIC', '-std=legacy',
           *[str(PESDIR / f) for f in _SRC_FILES], '-o', str(out)]
    if verbose:
        print('[bkmp2] building native library:\n  ' + ' '.join(cmd), flush=True)
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True)
    except FileNotFoundError as exc:
        raise BKMP2Unavailable(
            "gfortran was not found, and no prebuilt BKMP2 library exists.\n"
            "Install a Fortran compiler and rebuild:\n"
            f"  gfortran -O2 -shared -fPIC -std=legacy "
            f"{PESDIR}/Amod_bkmp2.f {PESDIR}/bkmp2_cwrap.f90 -o {out}\n"
            "(Debian/Ubuntu: sudo apt install gfortran)"
        ) from exc
    if proc.returncode != 0:
        raise BKMP2Unavailable(
            "Compilation of the BKMP2 library failed:\n"
            + (proc.stderr or proc.stdout or '')[-2000:]
        )
    if verbose:
        print(f'[bkmp2] built {out}', flush=True)
    return out


def get_lib(verbose=True, allow_build=True):
    """Return a cached ctypes handle to the BKMP2 library.

    Order: prebuilt file -> automatic gfortran build -> clear failure.
    """
    global _LIB
    if _LIB is not None:
        return _LIB

    path = find_library()
    if path is None:
        if not allow_build:
            raise BKMP2Unavailable(
                "No BKMP2 shared library found in "
                f"{PESDIR} and automatic building is disabled."
            )
        path = build_library(verbose=verbose)

    try:
        lib = ctypes.CDLL(str(path))
        lib.pot3_c.argtypes = [
            np.ctypeslib.ndpointer(np.float64, flags='C_CONTIGUOUS'),
            ctypes.POINTER(ctypes.c_double),
        ]
        lib.pot3_c.restype = None
    except OSError as exc:
        raise BKMP2Unavailable(
            f"Found {path} but it could not be loaded ({exc}).\n"
            "If it was built on another machine, rebuild it here with "
            "scripts/check_environment.py --rebuild."
        ) from exc

    # smoke test: isolated H2 at its equilibrium bond must be ≈ -0.17450 Eh
    # (= -38297.4 cm^-1, the H2 minimum on BKMP2)
    v = ctypes.c_double()
    lib.pot3_c(np.ascontiguousarray([1.401, 30.0, 30.0], dtype=np.float64),
               ctypes.byref(v))
    if not (-0.1755 < v.value < -0.1735):
        raise BKMP2Unavailable(
            f"{path} loaded but returned V = {v.value} Eh for isolated H2 "
            "(expected ≈ -0.17450). The library is likely stale or corrupted; "
            "remove it and rerun scripts/check_environment.py --rebuild."
        )
    _LIB = lib
    return lib


# ---------------------------------------------------------------------------
# potential evaluation
# ---------------------------------------------------------------------------
def V_h3_bohr(r1, r2, r3, lib=None):
    """BKMP2 energy (Hartree) for one geometry, distances in Bohr."""
    lib = lib or get_lib(verbose=False)
    v = ctypes.c_double()
    lib.pot3_c(np.ascontiguousarray([r1, r2, r3], dtype=np.float64),
               ctypes.byref(v))
    return v.value


def V_h3_angstrom(r1, r2, r3, lib=None, cm1=False):
    """BKMP2 energy for one geometry, distances in Angstrom (Eh or cm^-1)."""
    e = V_h3_bohr(r1 / BOHR, r2 / BOHR, r3 / BOHR, lib=lib)
    return e * HARTREE_TO_CM1 if cm1 else e


def potential_grid_rRg(r, R, gamma, cache_dir=None, verbose=True):
    """V(i_R, j_r, k_gamma) in cm^-1 on a (r, R, gamma) product grid.

    The H2 bond is frozen at ``r`` only in the *2D* helper below; here the
    three internuclear distances are reconstructed from the Jacobi geometry
    (atom A = H, diatom BC = H2 with bond length r).

    Shape follows ``np.meshgrid(R, r, gamma, indexing='ij')`` -> (nR, nr, nG).
    Results are cached to ``cache_dir`` keyed by the grid definition.
    """
    r = np.asarray(r, float)
    R = np.asarray(R, float)
    gamma = np.asarray(gamma, float)

    key = hashlib.sha1(
        np.concatenate([r, R, gamma]).tobytes()
        + f'|{r.size}x{R.size}x{gamma.size}|v1'.encode()
    ).hexdigest()[:16]
    cache_dir = Path(cache_dir) if cache_dir else (ROOT / 'results' / '_cache')
    cache_file = cache_dir / f'bkmp2_grid_{key}.npz'
    if cache_file.exists():
        if verbose:
            print(f'[bkmp2] grid cache hit: {cache_file.name}', flush=True)
        return np.load(cache_file)['V']

    lib = get_lib(verbose=verbose)
    rb = r / BOHR
    Rb = R / BOHR
    RR_b, rr_b, GG = np.meshgrid(Rb, rb, gamma, indexing='ij')
    d2 = np.sqrt(RR_b ** 2 + rr_b ** 2 / 4.0 - RR_b * rr_b * np.cos(GG))
    d3 = np.sqrt(RR_b ** 2 + rr_b ** 2 / 4.0 + RR_b * rr_b * np.cos(GG))
    out = np.empty(d2.shape, dtype=float)
    it = np.nditer([rr_b, d2, d3, out], [],
                   [['readonly']] * 3 + [['readwrite']])
    n = 0
    for a, b_, c_, z in it:
        z[...] = V_h3_bohr(float(a), float(b_), float(c_), lib=lib)
        n += 1
    V = out * HARTREE_TO_CM1
    try:
        cache_dir.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(cache_file, V=V)
        if verbose:
            print(f'[bkmp2] grid cache written: {cache_file.name} '
                  f'({n} points)', flush=True)
    except OSError:
        pass
    return V


def potential_2d_rRg(R, gamma, r_fix, cache_dir=None, verbose=True):
    """V(R, gamma) in cm^-1 with the diatomic bond frozen at ``r_fix``."""
    R = np.asarray(R, float)
    gamma = np.asarray(gamma, float)
    key = hashlib.sha1(
        np.concatenate([R, gamma]).tobytes()
        + f'|{R.size}x{gamma.size}|r{r_fix:.6f}|v1'.encode()
    ).hexdigest()[:16]
    cache_dir = Path(cache_dir) if cache_dir else (ROOT / 'results' / '_cache')
    cache_file = cache_dir / f'bkmp2_2d_{key}.npz'
    if cache_file.exists():
        if verbose:
            print(f'[bkmp2] 2D cache hit: {cache_file.name}', flush=True)
        return np.load(cache_file)['V']

    lib = get_lib(verbose=verbose)
    Rb = R / BOHR
    re_b = r_fix / BOHR
    GG = gamma
    RR_b, GGm = np.meshgrid(Rb, GG, indexing='ij')
    d2 = np.sqrt(RR_b ** 2 + re_b ** 2 / 4.0 - RR_b * re_b * np.cos(GGm))
    d3 = np.sqrt(RR_b ** 2 + re_b ** 2 / 4.0 + RR_b * re_b * np.cos(GGm))
    out = np.empty(d2.shape, dtype=float)
    it = np.nditer([d2, d3, out], [], [['readonly']] * 2 + [['readwrite']])
    for b_, c_, z in it:
        z[...] = V_h3_bohr(re_b, float(b_), float(c_), lib=lib)
    V = out * HARTREE_TO_CM1
    try:
        cache_dir.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(cache_file, V=V)
        if verbose:
            print(f'[bkmp2] 2D cache written: {cache_file.name}', flush=True)
    except OSError:
        pass
    return V


def asymptotic_energy(r_fix=0.74144, R_far=30.0):
    """Isolated-diatom asymptote V(H2) in cm^-1 at bond length ``r_fix``."""
    return V_h3_angstrom(r_fix, R_far, R_far, cm1=True)


if __name__ == '__main__':      # pragma: no cover - manual smoke test
    print('library :', find_library())
    get_lib(verbose=True)
    print('V(H2, re)      =', asymptotic_energy(), 'cm^-1')
    print('V(1.401,30,30) =', V_h3_bohr(1.401, 30.0, 30.0), 'Eh')
