#!/usr/bin/env python3
"""Server-side convergence study of the full-dimensional 3D A+BC reference.

The 3D H + H2(v=0) / BKMP2 reference is the unconstrained Kronecker
finite-difference / Legendre-DVR solve of the (r, R, gamma) Schrodinger
equation.  It is itself numerical, so this script measures ITS OWN error by
grid refinement and by monitoring the elastic asymptotic channel amplitude

    |A(R)| = |<v=0,j=0 | chi(R)>|   at R = 8 A.

Configurations (all with the consistent 4th-order operator, exact Dirichlet at
R_min and both r-ends, 4th-order multichannel Robin at R_max):

  baseline   nr=40, nR=120, nG=16   r in [0.40, 1.9]   R in [0.7, 8.5]
  nR240      in-plane radial refinement
  nG20       angular DVR refinement
  nr60       vibrational-grid refinement
  Rbox10     R_max = 10 A (same hR)
  rbox22     r_max = 2.2 A (same hr)

Usage
-----
    python3 scripts/server_3d_study.py --launch --energy 400      # parallel
    python3 scripts/server_3d_study.py --aggregate --energy 400   # table
"""
import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))

STORE = ROOT / 'results' / '3d_study'

CONFIGS = {
    'baseline': dict(nr=40, nR=120, nG=16, r_min=0.40, r_max=1.9, R_max=8.5),
    'nR240':    dict(nr=40, nR=240, nG=16, r_min=0.40, r_max=1.9, R_max=8.5),
    'nG20':     dict(nr=40, nR=120, nG=20, r_min=0.40, r_max=1.9, R_max=8.5),
    'nr60':     dict(nr=60, nR=120, nG=16, r_min=0.40, r_max=1.9, R_max=8.5),
    'Rbox10':   dict(nr=40, nR=156, nG=16, r_min=0.40, r_max=1.9, R_max=10.0),
    'rbox22':   dict(nr=52, nR=120, nG=16, r_min=0.40, r_max=2.2, R_max=8.5),
}


def run_one(label, energy):
    from run_scatter_3d import solve_3d_scattering
    from check_3d_reference_convergence import diatomic_channel, amplitude

    kw = dict(CONFIGS[label])
    kw['R_min'] = 0.7
    print(f'[study] {label}: {kw}  E_c = {energy}', flush=True)
    r, R, gamma, chi_ref, chi_pred, K0, l2, el = solve_3d_scattering(energy, **kw)
    v0, eps0 = diatomic_channel(r, gamma, R[-1], legacy_fd=False)
    amp8 = amplitude(chi_ref, R, 8.0, v0)
    STORE.mkdir(parents=True, exist_ok=True)
    dest = STORE / f'{label}_{int(energy)}.npz'
    np.savez_compressed(
        dest, r=r, R=R, gamma=gamma, chi=chi_ref,
        K0=K0, l2_prod=l2, elapsed=el, eps0=eps0, amp8=amp8,
        max_chi=float(np.max(np.abs(chi_ref))))
    print(f'[study] {label}: done in {el:.0f} s, l2_prod = {l2*100:.4f}%, '
          f'eps0 = {eps0:.3f} cm^-1, |A(8A)| = {amp8:.8f} -> {dest.name}', flush=True)


def launch(energy, labels):
    STORE.mkdir(parents=True, exist_ok=True)
    (STORE / 'logs').mkdir(exist_ok=True)
    procs = []
    for label in labels:
        cfg = CONFIGS[label]
        # skip if already present
        if (STORE / f'{label}_{int(energy)}.npz').exists():
            print(f'[study] {label}: already present, skipping')
            continue
        log = open(STORE / 'logs' / f'{label}_{int(energy)}.log', 'w')
        p = subprocess.Popen([sys.executable, str(Path(__file__).resolve()),
                              '--run', label, '--energy', str(energy)],
                             stdout=log, stderr=subprocess.STDOUT, cwd=str(ROOT))
        procs.append((label, p, log))
        print(f'[study] launched {label} (pid {p.pid})')
    return procs


def aggregate(energy):
    from check_3d_reference_convergence import interp_to_target

    data = {}
    for label in CONFIGS:
        f = STORE / f'{label}_{int(energy)}.npz'
        if f.exists():
            data[label] = dict(np.load(f))
    if 'baseline' not in data:
        print('baseline result missing - run --launch first')
        return 1
    base = data['baseline']
    print('=' * 96)
    print(f'3D reference convergence study  (E_c = {energy:.1f} cm^-1)')
    print('=' * 96)
    print(f'{"config":10s} {"unknowns":>9s} {"t (s)":>7s} {"max|chi|":>9s} '
          f'{"eps0 (cm^-1)":>13s} {"|A(8A)|":>11s} {"relL2 vs base":>14s} {"l2_prod":>9s}')
    chi_b, r_b, R_b, g_b = base['chi'], base['r'], base['R'], base['gamma']
    norm_b = np.linalg.norm(chi_b)
    rows = {}
    for label, d in data.items():
        if label == 'baseline':
            dl2 = 0.0
        else:
            on_base = interp_to_target(d['chi'], d['r'], d['R'], d['gamma'],
                                       r_b, R_b, g_b)
            dl2 = float(np.linalg.norm(on_base - chi_b) / norm_b)
        amp = float(d['amp8'])
        rows[label] = dict(unknowns=int(np.prod(d['chi'].shape)),
                           seconds=float(d['elapsed']), max_chi=float(d['max_chi']),
                           eps0=float(d['eps0']), amp8=amp,
                           relL2_vs_baseline=dl2, l2_prod=float(d['l2_prod']),
                           amp8_diff_vs_baseline=float(abs(amp - base['amp8']) / base['amp8']))
        print(f'{label:10s} {rows[label]["unknowns"]:9d} {rows[label]["seconds"]:7.0f} '
              f'{rows[label]["max_chi"]:9.6f} {rows[label]["eps0"]:13.3f} '
              f'{amp:11.7f} {dl2:14.3e} {rows[label]["l2_prod"]*100:8.4f}%')
    ref_err = max((v['relL2_vs_baseline'] for k, v in rows.items() if k != 'baseline'),
                  default=0.0)
    amp_err = max((v['amp8_diff_vs_baseline'] for k, v in rows.items() if k != 'baseline'),
                  default=0.0)
    print(f'\nreference self-convergence: worst field change under refinement '
          f'= {ref_err*100:.4f}%   worst |A(8A)| change = {amp_err*100:.4f}%')
    print(f'channel-projected vs full (production grid): '
          f'{rows["baseline"]["l2_prod"]*100:.4f}%')
    out = {'E_c': energy, 'rows': rows, 'reference_self_convergence': ref_err,
           'reference_amp_convergence': amp_err, 'configs': {k: v for k, v in CONFIGS.items()}}
    dest = ROOT / 'results' / f'check_3d_reference_convergence_{int(energy)}.json'
    with open(dest, 'w') as f:
        json.dump(out, f, indent=2)
    print(f'saved {dest}')
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--energy', type=float, default=400.0)
    ap.add_argument('--launch', action='store_true')
    ap.add_argument('--aggregate', action='store_true')
    ap.add_argument('--run', default=None, help=argparse.SUPPRESS)
    ap.add_argument('--labels', nargs='+', default=None)
    args = ap.parse_args(argv)

    if args.run:
        run_one(args.run, args.energy)
        return 0
    if args.launch:
        launch(args.energy, args.labels or list(CONFIGS))
        return 0
    if args.aggregate:
        return aggregate(args.energy)
    ap.print_help()
    return 0


if __name__ == '__main__':
    sys.exit(main())
