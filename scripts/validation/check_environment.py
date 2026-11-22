#!/usr/bin/env python3
"""Environment self-check for the PINNSlove package.

Verifies, in order:
  1. Python / NumPy / SciPy / (optional) PyTorch / Matplotlib versions
  2. NumPy 2.0 ``trapezoid`` compatibility shim
  3. Native BKMP2 library: found, loadable, or auto-buildable (--rebuild to force)
  4. The BKMP2 Hamiltonian probe: V(H2, re) == -38297.4 cm^-1
  5. Classical reference solver smoke test (1D Morse scattering)
  6. Coupled-channel machinery smoke test (hard-sphere phase shift)
  7. Writable output directories (results/, article/fig/)

Exit code 0 when every required check passes, 1 otherwise.

Usage
-----
    python3 scripts/check_environment.py            # diagnose
    python3 scripts/check_environment.py --rebuild  # force-rebuild the BKMP2 lib
"""
import argparse
import platform
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

OK = '  [ok]  '
WARN = ' [warn] '
FAIL = ' [FAIL] '

results = []


def report(tag, msg):
    print(tag + msg, flush=True)
    results.append((tag, msg))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--rebuild', action='store_true',
                    help='delete and rebuild the BKMP2 shared library')
    args = ap.parse_args(argv)

    print('=' * 74)
    print('PINNSlove environment check')
    print('=' * 74)

    # --- 1. interpreters and libraries ----------------------------------
    report(OK, f'Python {platform.python_version()} ({platform.platform()})')
    try:
        import numpy as np
        report(OK, f'NumPy {np.__version__}')
    except Exception as exc:
        report(FAIL, f'NumPy import failed: {exc}')
        return finish()

    try:
        import scipy
        report(OK, f'SciPy {scipy.__version__}')
    except Exception as exc:
        report(FAIL, f'SciPy import failed: {exc}')
        return finish()

    from src.compat import has_torch, trapezoid
    if has_torch():
        import torch
        cuda = torch.cuda.is_available()
        report(OK, f'PyTorch {torch.__version__} (CUDA available: {cuda})')
    else:
        report(WARN, 'PyTorch not installed - the PINN solvers are unavailable, '
                     'but all classical reference / coupled-channel tools work.')

    try:
        import matplotlib
        report(OK, f'Matplotlib {matplotlib.__version__}')
    except Exception:
        report(WARN, 'Matplotlib not installed - figure generation will be skipped.')

    # --- 2. NumPy 2.0 shim ----------------------------------------------
    x = np.linspace(0.0, 1.0, 11)
    y = x ** 2
    val = float(trapezoid(y, x))
    ref = float(0.5 * np.sum((x[1:] - x[:-1]) * (y[1:] + y[:-1])))
    if abs(val - ref) < 1e-14:
        report(OK, f'trapezoid shim works (NumPy 2.0 API on {np.__version__})')
    else:
        report(FAIL, f'trapezoid shim returned {val} (expected {ref})')

    # --- 3. BKMP2 native library ----------------------------------------
    from src import bkmp2
    if args.rebuild:
        for p in bkmp2.PESDIR.glob('libbkmp2*'):
            p.unlink()
            print(f'        removed {p.name}')
    lib = bkmp2.find_library()
    if lib:
        report(OK, f'BKMP2 library found: {lib.name}')
    else:
        report(WARN, 'No BKMP2 shared library yet - attempting to build one')
    try:
        bkmp2.get_lib(verbose=True)
        report(OK, 'BKMP2 library loads and passes the H2 probe')
    except bkmp2.BKMP2Unavailable as exc:
        report(FAIL, 'BKMP2 unavailable: ' + str(exc).splitlines()[0])
        print('        (the 1D/2D model-water calculations still work; '
              'only the H+H2 BKMP2 scripts are affected)')

    # --- 4. potential probe ---------------------------------------------
    try:
        v = bkmp2.asymptotic_energy()
        if abs(v - (-38297.4)) < 1.0:
            report(OK, f'BKMP2 H2 minimum: {v:.2f} cm^-1 (reference -38297.4)')
        else:
            report(FAIL, f'BKMP2 H2 minimum {v:.2f} cm^-1 is off reference')
    except Exception as exc:
        report(FAIL, f'BKMP2 potential probe failed: {exc}')

    # --- 5. classical reference solver ----------------------------------
    try:
        from src.potentials import Morse
        from src.reference import scatter_fd_1d
        mo = Morse(38293.0, 1.94286, 0.74144, 0.50397)
        _, _, R_, _ = scatter_fd_1d(mo, 2000.0, -0.5, 6.0, npts=4000)
        if abs(abs(R_) - 1.0) < 1e-6:
            report(OK, f'1D Morse scattering reference solver (|R| = {abs(R_):.6f})')
        else:
            report(FAIL, f'1D reference solver returned |R| = {abs(R_):.6f}')
    except Exception as exc:
        report(FAIL, f'1D reference solver failed: {exc}')

    # --- 6. coupled-channel machinery -----------------------------------
    try:
        from src.j_frames import JAngularCoupling
        from src.scatter_extract import solve_radial_single
        from scipy.special import spherical_jn, spherical_yn
        c = JAngularCoupling(J=1, j_max=3)
        U, _, bf, _ = c.full_frame_transformation(parity=-1)
        ortho = float(abs(U.T @ U - np.eye(len(bf))).max())
        S, *_ = solve_radial_single(25.0, 0, 25.0, 1.0, 40.0, 8001)
        jl, nl = spherical_jn(0, 1.0), spherical_yn(0, 1.0)
        err = abs(S - (nl + 1j * jl) / (nl - 1j * jl))
        if ortho < 1e-12 and err < 1e-8:
            report(OK, f'coupled-channel machinery (frame orthogonality '
                       f'{ortho:.1e}, hard-sphere error {err:.1e})')
        else:
            report(FAIL, f'coupled-channel machinery off '
                         f'(ortho {ortho:.1e}, error {err:.1e})')
    except Exception as exc:
        report(FAIL, f'coupled-channel machinery failed: {exc}')

    # --- 7. high-dimensional machinery (4D/6D vibrations, 4D diatom-diatom) ---
    try:
        from src.hamiltonians_nd import CoupledAnharmonicND
        from src.pinn_high_dim import HighDimVariationalPINN
        from src.diatom_diatom_4d import DiatomDiatom4DScattering

        # harmonic limit of the QMC variational estimator (bounded-measure weights)
        prob = CoupledAnharmonicND(dim=4)
        qmc = HighDimVariationalPINN(prob, n_samples=2048, max_poly_deg=2)
        rel = qmc.self_check(tol=0.02)

        # validated eigensolver: residual of the dual-scheme benchmark
        _ev, _vc, res = prob.solve_fci(n_states_per_mode=4, k_eigenvalues=3)

        # 4D diatom-diatom coupled-channel solve + flux conservation
        scat = DiatomDiatom4DScattering(Ec=400.0, nR=200, j_max=2)
        out = scat.solve_scattering()
        uni = out['unitarity_residual']
        if rel < 0.02 and res < 1e-8 and uni < 1e-4:
            report(OK, f'high-dimensional machinery (QMC harmonic self-check {rel:.1e}, '
                       f'benchmark residual {res:.1e}, 4D scattering unitarity {uni:.1e})')
        else:
            report(FAIL, f'high-dimensional machinery off (self-check {rel:.1e}, '
                         f'residual {res:.1e}, unitarity {uni:.1e})')
    except Exception as exc:
        report(FAIL, f'high-dimensional machinery failed: {exc}')

    # --- 8. NumPy 2.x generator-removal audit ---------------------------
    try:
        import re
        src_root = Path(__file__).resolve().parents[2]
        bad = []
        for f in list((src_root / 'src').rglob('*.py')) + list((src_root / 'scripts').rglob('*.py')):
            txt = f.read_text()
            for m in re.finditer(r'np\.(sum|mean|max|min|any|all)\(\s*[A-Za-z_][\w\.\[\]\(\)]*\s+for\s+', txt):
                bad.append(f'{f.name}:{txt[:m.start()].count(chr(10)) + 1}')
        if bad:
            report(FAIL, 'np.sum(generator) removed in NumPy 2.x found at: ' + ', '.join(bad))
        else:
            report(OK, 'NumPy 2.x compatibility (no np.sum(generator) patterns)')
    except Exception as exc:
        report(WARN, f'NumPy 2.x audit skipped: {exc}')

    # --- 9. writable output directories ---------------------------------
    for d in [ROOT / 'results', ROOT / 'article' / 'fig', ROOT / 'results' / '_cache']:
        try:
            d.mkdir(parents=True, exist_ok=True)
            probe = d / '.write_probe'
            probe.write_text('ok')
            probe.unlink()
            report(OK, f'writable: {d.relative_to(ROOT)}/')
        except Exception as exc:
            report(WARN, f'not writable: {d} ({exc})')

    return finish()


def finish():
    print('=' * 74)
    n_fail = sum(1 for t, _ in results if t == FAIL)
    n_warn = sum(1 for t, _ in results if t == WARN)
    if n_fail:
        print(f'RESULT: {n_fail} failure(s), {n_warn} warning(s) - '
              'see the messages above.')
    elif n_warn:
        print(f'RESULT: environment usable with {n_warn} warning(s).')
    else:
        print('RESULT: all checks passed.')
    print('=' * 74)
    return 1 if n_fail else 0


if __name__ == '__main__':
    sys.exit(main())
