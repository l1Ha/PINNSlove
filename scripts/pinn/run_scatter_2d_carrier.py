"""Benchmark: adiabatic Langer-Airy carrier (+ optional smooth correction)
against the converged CAP FD reference.

Reports, honestly and without any post-hoc scaling:
  * validation PDE residual (the only quantity used for model selection);
  * interior L2 error against the reference (400x60, eta* from the plateau
    scan) alongside the reference's own uncertainty;
  * the honest gate outcome: whether the smooth correction was accepted on
    the held-out PDE residual (if not, the pure carrier is reported).
"""
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from src.potentials import WaterPES
from src.pinn.pinn_scatter_2d_carrier import CarrierPINN2D
from run_scatter_2d_cap_benchmark import compare, fd

ETA_DEFAULT = {200.0: 100.0, 500.0: 200.0}


def main():
    pes = WaterPES()
    root = Path(__file__).resolve().parents[2]
    eta_file = root / 'results' / 'scatter2d_cap_benchmark.json'
    eta_star = {}
    if eta_file.exists():
        for rec in json.load(open(eta_file))['energies']:
            eta_star[rec['E_c']] = rec['eta_star']
    out = {'energies': []}
    for E_c in (200.0, 500.0):
        eta = float(eta_star.get(E_c, ETA_DEFAULT[E_c]))
        chi_ref, R, g, K, E_tot = fd(pes, E_c, 400, 60, eta)
        chi_lo, R2, g2, *_ = fd(pes, E_c, 250, 50, eta)
        ref_unc = compare(chi_ref, R, g, chi_lo, R2, g2)
        rec = {'E_c': E_c, 'eta_star': eta, 'reference_uncertainty_l2': ref_unc,
               'runs': []}
        print(f'--- E_c={E_c:.0f} eta*={eta} ref_unc={ref_unc:.4f}', flush=True)
        for n_ch, n_corr in [(8, 0), (8, 300), (8, 1500)]:
            t0 = time.time()
            sol = CarrierPINN2D(pes, E_c, n_R=400, n_gamma=60,
                                cap_start=5.0, cap_eta=eta, cap_width=2.0,
                                n_channels=n_ch, n_corr=n_corr,
                                sigma_R=(8.0, 20.0), sigma_g=6.0, seed=0)
            sv = sol.solve(verbose=False)
            chi = sol.evaluate()
            l2_49 = compare(chi, sol.R, sol.gam, chi_ref, R, g, 4.9)
            l2_69 = compare(chi, sol.R, sol.gam, chi_ref, R, g, 6.9)
            run = dict(n_channels=n_ch, n_corr=n_corr,
                       carrier=sv['carrier'], correction=sv['correction'],
                       l2_interior_4p9=l2_49, l2_interior_6p9=l2_69,
                       wall_s=round(time.time() - t0, 1))
            rec['runs'].append(run)
            print(f'  n_ch={n_ch} n_corr={n_corr}: '
                  f'carrier_val={sv["carrier"]["val_residual"]:.3e} '
                  f'corr_acc={sv["correction"]["accepted"]} '
                  f'corr_val={sv["correction"]["val_residual"]:.3e} '
                  f'L2(4.9)={l2_49:.4f} L2(6.9)={l2_69:.4f} '
                  f'({run["wall_s"]}s)', flush=True)
            np.savez(f'results/_scatter2d_carrier_{int(E_c)}_{n_ch}_{n_corr}.npz',
                     R=sol.R, gamma=sol.gam, chi_ref=chi_ref, chi_carrier=chi)
        out['energies'].append(rec)
        json.dump(out, open('results/scatter2d_carrier.json', 'w'), indent=2)
    print('carrier 基准完成 -> results/scatter2d_carrier.json', flush=True)


if __name__ == '__main__':
    main()
