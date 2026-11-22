#!/usr/bin/env python3
"""Production-grade full-dimensional 3D A+BC scattering driver (H + H2 on BKMP2).

Runs each collision energy as an independent single-threaded process (well suited
to a cluster: `--launch` fans out over the 192-core node), then merges everything
into the standard summary JSON and regenerates all publication figures.

Production configuration (much finer than the exploratory baseline 40x120x16):
    nr = 64   vibrational points on r in [0.40, 2.20] A
    nR = 256  collision points   on R in [0.70, 10.00] A
    nG = 20   Gauss-Legendre angular DVR points
    -> 64*256*20 = 327,680 coupled complex unknowns solved by sparse direct LU

Usage
-----
    python3 scripts/run_3d_production.py --launch            # fan out all energies
    python3 scripts/run_3d_production.py --run --energy 400  # single energy worker
    python3 scripts/run_3d_production.py --aggregate         # merge + figures
"""
import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))

STORE = ROOT / 'results' / '3d_prod'
ENERGIES = [100.0, 200.0, 300.0, 400.0, 600.0, 800.0, 1000.0, 1200.0]

CONFIGS = {
    # production grid: converged reference, wide boxes, fine angular DVR
    'production': dict(nr=64, nR=256, nG=20, r_max=2.20, R_max=10.0),
    # smoke test used for pre-flight checks
    'smoke':      dict(nr=64, nR=64, nG=20, r_max=2.20, R_max=10.0),
    'baseline':   dict(nr=40, nR=120, nG=16, r_max=1.90, R_max=8.5),
}


def _worker(energy, cfg_name='production', with_figure=True):
    from run_scatter_3d import solve_3d_scattering, extract_s00_3d, plot_3d_slices

    cfg = CONFIGS[cfg_name]
    t0 = time.time()
    (r, R, gamma, chi_ref, chi_pred, K0,
     l2_err, elapsed) = solve_3d_scattering(
        energy, nr=cfg['nr'], nR=cfg['nR'], nG=cfg['nG'],
        r_min=0.40, r_max=cfg['r_max'], R_min=0.7, R_max=cfg['R_max'])

    S00_ref = extract_s00_3d(chi_ref, R, K0, r, gamma, cfg['nr'], cfg['nG'])
    S00_pred = extract_s00_3d(chi_pred, R, K0, r, gamma, cfg['nr'], cfg['nG'])
    asym_ref = float(np.mean(np.abs(chi_ref[R > 0.75 * cfg['R_max']])))
    asym_pred = float(np.mean(np.abs(chi_pred[R > 0.75 * cfg['R_max']])))

    rec = {
        'E_c_cm1': energy,
        'K0_inv_A': float(K0),
        'l2_relative_error_3d': float(l2_err),
        'asym_amplitude_fd': asym_ref,
        'asym_amplitude_pred': asym_pred,
        'S00_ref_abs': float(abs(S00_ref)),
        'S00_pred_abs': float(abs(S00_pred)),
        'S00_ref_arg_deg': float(np.degrees(np.angle(S00_ref))),
        'S00_pred_arg_deg': float(np.degrees(np.angle(S00_pred))),
        'P_elastic': float(abs(S00_ref) ** 2),
        'P_reaction_inelastic': float(max(1.0 - abs(S00_ref) ** 2, 0.0)),
        'grid': {'nr': cfg['nr'], 'nR': cfg['nR'], 'nG': cfg['nG'],
                 'r_max': cfg['r_max'], 'R_max': cfg['R_max']},
        'unknowns': int(cfg['nr'] * cfg['nR'] * cfg['nG']),
        'elapsed_sec': float(elapsed),
    }

    STORE.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        STORE / f'prod_{int(energy)}.npz',
        r=r, R=R, gamma=gamma, chi_ref=chi_ref, chi_pred=chi_pred,
        metrics=json.dumps(rec))

    if with_figure:
        from run_scatter_3d import plot_3d_slices
        fig_path = ROOT / 'article' / 'fig' / f'fig_scatter_3d_bkmp2_{int(energy)}.png'
        plot_3d_slices(r, R, gamma, chi_ref, chi_pred, energy, K0, l2_err, fig_path)

    print(f"[prod] E_c={energy:6.1f} cm^-1 done in {elapsed:6.1f}s "
          f"({rec['unknowns']} unknowns) | relL2 = {l2_err*100:.4f}% | "
          f"|S00| ref = {abs(S00_ref):.6f} pred = {abs(S00_pred):.6f} | "
          f"arg ref = {rec['S00_ref_arg_deg']:.4f} deg", flush=True)
    return rec


def launch(energies, cfg_name='production'):
    STORE.mkdir(parents=True, exist_ok=True)
    (STORE / 'logs').mkdir(exist_ok=True)
    for E in energies:
        out = STORE / f'prod_{int(E)}.npz'
        if out.exists():
            print(f'[prod] E_c={E}: already present, skipping')
            continue
        log = open(STORE / 'logs' / f'prod_{int(E)}.log', 'w')
        subprocess.Popen(
            [sys.executable, str(Path(__file__).resolve()),
             '--run', '--energy', str(E), '--config', cfg_name],
            stdout=log, stderr=subprocess.STDOUT, cwd=str(ROOT))


def aggregate(make_figures=True):
    recs = []
    chis = {}
    for f in sorted(STORE.glob('prod_*.npz')):
        d = np.load(f)
        recs.append(json.loads(str(d['metrics'])))
        chis[int(round(recs[-1]['E_c_cm1']))] = d
    if not recs:
        print('[prod] nothing to aggregate - run --launch first')
        return 1
    recs.sort(key=lambda r: r['E_c_cm1'])

    OUT = ROOT / 'results'
    # keep the exploratory baseline around
    base = OUT / 'scatter_3d_bkmp2_summary.json'
    if base.exists():
        (OUT / 'scatter_3d_bkmp2_summary_baseline.json').write_text(base.read_text())
    with open(base, 'w') as f:
        json.dump(recs, f, indent=2)

    err = np.array([r['l2_relative_error_3d'] for r in recs]) * 100.0
    print('=' * 92)
    print('PRODUCTION full-dimensional 3D scan (nr=64, nR=256, nG=20, r<=2.2 A, R<=10 A)')
    print('=' * 92)
    print(f"{'E_c (cm^-1)':>12s} {'unknowns':>10s} {'t (s)':>7s} {'relL2 (%)':>10s} "
          f"{'|S00| ref':>10s} {'|S00| pred':>11s} {'argS00 (deg)':>13s}")
    for r in recs:
        print(f"{r['E_c_cm1']:12.1f} {r.get('unknowns', 0):10d} {r['elapsed_sec']:7.0f} "
              f"{r['l2_relative_error_3d']*100:10.4f} {r['S00_ref_abs']:10.6f} "
              f"{r['S00_pred_abs']:11.6f} {r['S00_ref_arg_deg']:13.4f}")
    print(f"\nmean relL2 = {err.mean():.4f}%   max relL2 = {err.max():.4f}%   "
          f"all < 1%: {bool((err < 1.0).all())}")

    if make_figures:
        sys.path.insert(0, str(ROOT / 'scripts'))
        from plot_3d_energy_scan import main as plot_scan
        plot_scan()
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--launch', action='store_true')
    ap.add_argument('--run', action='store_true')
    ap.add_argument('--aggregate', action='store_true')
    ap.add_argument('--energy', type=float, default=None)
    ap.add_argument('--energies', type=float, nargs='+', default=ENERGIES)
    ap.add_argument('--config', default='production', choices=sorted(CONFIGS))
    ap.add_argument('--no-figures', action='store_true')
    args = ap.parse_args(argv)

    if args.run:
        assert args.energy is not None
        rec = _worker(args.energy, args.config, with_figure=not args.no_figures)
        print(json.dumps(rec))
        return 0
    if args.launch:
        launch(args.energies, args.config)
        return 0
    if args.aggregate:
        return aggregate(make_figures=not args.no_figures)
    ap.print_help()
    return 0


if __name__ == '__main__':
    sys.exit(main())
