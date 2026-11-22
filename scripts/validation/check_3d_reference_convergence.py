#!/usr/bin/env python3
"""Convergence study of the FULL-DIMENSIONAL 3D A+BC reference solution.

The reference for the 3D H + H2(v=0) / BKMP2 scattering problem is the
unconstrained Kronecker finite-difference / Gauss-Legendre-DVR solve of

  [-Br d2/dr2 - BR d2/dR2 + (Br/r^2 + BR/R^2) j^2(gamma) + V(r,R,gamma) - E] chi = 0

with Dirichlet at R_min and both r-ends and a multichannel incoming Robin
condition at R_max.  No analytic solution exists for this surface, so the
reference is itself numerical; this script quantifies ITS OWN discretization
error by

  * comparing the consistent 4th-order operator against the legacy operator
    (whose wall-adjacent rows are inconsistent and do not converge),
  * grid refinement in nR, nr, nG,
  * monitoring a physical observable - the elastic asymptotic channel
    amplitude |A(R)| = |<v=0,j=0 | chi(R)>| at R = 8 A.

Usage
-----
    python3 scripts/check_3d_reference_convergence.py --quick
    python3 scripts/check_3d_reference_convergence.py --energies 400 800
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import scipy.linalg
from scipy.interpolate import RegularGridInterpolator

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))

from run_scatter_3d import solve_3d_scattering, Br, BR          # noqa: E402
from src.reference import LegendreAngularDVR, _fd_second          # noqa: E402
from src.scatter_extract import fd2_4th_dirichlet                 # noqa: E402
from src import bkmp2                                            # noqa: E402


def diatomic_channel(r, gamma, R_probe, legacy_fd=False):
    """Ground diatomic channel vector (flattened (nr, nG)) at R = R_probe.

    With the consistent operator the diatomic eigenproblem is solved on the
    interior r-points (Dirichlet walls) and padded with zeros, exactly as in
    :func:`run_scatter_3d.solve_3d_scattering`.
    """
    nr, nG = len(r), len(gamma)
    hr = r[1] - r[0]
    ang = LegendreAngularDVR(nG)
    V_slice = bkmp2.potential_grid_rRg(r, np.array([R_probe]), gamma).reshape(nr, nG)
    if legacy_fd:
        D2 = _fd_second(nr, hr, 4).toarray()
        H_d = (np.kron(-Br * D2, np.eye(nG))
               + np.kron(np.diag(Br / r ** 2), ang.j2)
               + np.diag(V_slice.ravel()))
        evals, evecs = np.linalg.eigh(H_d)
        idx = np.argsort(evals)
        return evecs[:, idx[0]], evals[idx[0]]
    D2 = fd2_4th_dirichlet(nr, hr)
    ri = np.arange(1, nr - 1)
    D2_ii = np.asarray(D2[1:nr - 1, 1:nr - 1].todense())
    H_red = (np.kron(-Br * D2_ii, np.eye(nG))
             + np.kron(np.diag(Br / r[ri] ** 2), ang.j2)
             + np.diag(V_slice[ri, :].ravel()))
    evals, evecs = scipy.linalg.eig(H_red)          # non-symmetric operator
    idx = np.argsort(evals.real)
    v0 = np.zeros(nr * nG)
    for ir in range(1, nr - 1):
        v0[ir * nG:(ir + 1) * nG] = evecs[(ir - 1) * nG:ir * nG, idx[0]].real
    return v0, evals.real[idx[0]]


def amplitude(chi, R, R_probe, v0):
    """|<v0|chi(R_probe)>| at the radial grid point nearest R_probe."""
    iR = int(np.argmin(np.abs(R - R_probe)))
    return float(abs(np.vdot(v0, chi[iR, :, :].ravel())))


def interp_to_target(chi_c, r_c, R_c, g_c, r_t, R_t, g_t):
    """Interpolate a complex 3D field onto another product grid.

    Real and imaginary parts are interpolated separately; the output array is
    zero-initialized (np.empty here leaves NaN garbage in the accumulations).
    """
    out = np.zeros((len(R_t), len(r_t), len(g_t)), dtype=complex)
    pts = np.stack(np.meshgrid(R_t, r_t, g_t, indexing='ij'), axis=-1)
    for factor, comp in ((1.0, chi_c.real), (1j, chi_c.imag)):
        itp = RegularGridInterpolator((R_c, r_c, g_c), comp, method='linear',
                                      bounds_error=False, fill_value=0.0)
        out += factor * itp(pts)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--energies', type=float, nargs='+', default=[400.0])
    ap.add_argument('--quick', action='store_true',
                    help='skip the expensive refinements (legacy/fixed only)')
    ap.add_argument('--output', default=None)
    args = ap.parse_args(argv)

    OUT = ROOT / 'results'
    OUT.mkdir(exist_ok=True)
    records = []

    for Ec in args.energies:
        print('=' * 92)
        print(f'3D reference convergence study,  Ec = {Ec:.1f} cm^-1')
        print('=' * 92)

        grids = [('legacy 28x120x16', dict(nR=120), True),
                 ('fixed  28x120x16', dict(nR=120), False)]
        if not args.quick:
            grids += [('fixed  28x240x16', dict(nR=240), False),
                      ('fixed  28x120x20', dict(nG=20), False),
                      ('fixed  40x120x16', dict(nr=40), False)]

        sols = {}
        for label, kw, legacy in grids:
            print(f'\n--- {label}  (legacy_fd={legacy}) ---', flush=True)
            r, R, gamma, chi, _, K0, l2, el = solve_3d_scattering(
                Ec, legacy_fd=legacy, **kw)
            v0, eps0 = diatomic_channel(r, gamma, R[-1], legacy_fd=legacy)
            amp8 = amplitude(chi, R, 8.0, v0)
            sols[label] = dict(r=r, R=R, gamma=gamma, chi=chi, K0=K0,
                               amp8=amp8, maxchi=float(np.max(np.abs(chi))),
                               sec=el, l2_prod=l2, legacy=legacy)
            print(f'    {np.prod(chi.shape)} unknowns, {el:.0f} s | '
                  f'max|chi| = {np.max(np.abs(chi)):.6f} | '
                  f'|A(8A)| = {amp8:.8f} | '
                  f'channel-projected rel. L2 = {l2*100:.4f}%', flush=True)

        base = sols['fixed  28x120x16']
        row = {'E_c': Ec, 'grids': {}}
        for label, s in sols.items():
            row['grids'][label] = dict(
                unknowns=int(np.prod(s['chi'].shape)), seconds=float(s['sec']),
                max_chi=s['maxchi'], amp8=s['amp8'],
                l2_channel_projected=float(s['l2_prod']))
        # legacy vs fixed, on the common (28x120x16) grid
        d_leg = sols['legacy 28x120x16']
        dl = np.linalg.norm(d_leg['chi'] - base['chi']) / np.linalg.norm(base['chi'])
        row['legacy_vs_fixed_L2'] = float(dl)
        print(f'\nlegacy vs fixed (same grid, field difference): rel. L2 = {dl:.3e}')
        print(f'  |A(8A)| legacy = {d_leg["amp8"]:.8f}   fixed = {base["amp8"]:.8f}'
              f'   (difference {abs(d_leg["amp8"]-base["amp8"])/base["amp8"]*100:.4f}%)')

        # refinements, interpolated onto the baseline grid
        row['refinement_vs_baseline_L2'] = {}
        for label, s in sols.items():
            if label in ('legacy 28x120x16', 'fixed  28x120x16'):
                continue
            fine_on_base = interp_to_target(s['chi'], s['r'], s['R'], s['gamma'],
                                            base['r'], base['R'], base['gamma'])
            d = np.linalg.norm(fine_on_base - base['chi']) / np.linalg.norm(base['chi'])
            row['refinement_vs_baseline_L2'][label] = float(d)
            print(f'{label} -> baseline grid: rel. L2 = {d:.3e} | '
                  f'|A(8A)| diff = {abs(s["amp8"]-base["amp8"])/base["amp8"]*100:.4f}%')
        records.append(row)

    dest = Path(args.output) if args.output else OUT / 'check_3d_reference_convergence.json'
    with open(dest, 'w') as f:
        json.dump(records, f, indent=2)
    print(f'\nSaved {dest}')

    legacy_err = max(r['legacy_vs_fixed_L2'] for r in records)
    refin = max([v for r in records for v in r['refinement_vs_baseline_L2'].values()],
                default=0.0)
    print(f'\nworst legacy-vs-consistent operator difference : {legacy_err:.3e}')
    print(f'worst refinement residual (consistent operator): {refin:.3e}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
