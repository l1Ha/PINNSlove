"""Diabatic vs adiabatic carrier benchmark against the converged CAP reference.

All configurations share the converged FD reference (400x60, eta* from the
plateau scan) and the same honest protocol:
  * model selection (ridge lambda, correction acceptance) uses ONLY the PDE
    validation residual on held-out nodes;
  * reported errors are exterior L2 against the reference, alongside the
    reference's own grid uncertainty;
  * the diabatic-transform diagnostics (antisymmetry of A, orthonormality of
    U, off-diagonal ratio of Veff) are printed and saved.

Runs: adiabatic carrier + 300 correction features (previous best), diabatic
carrier alone, diabatic carrier + 300 correction features.
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
from src.pinn.pinn_scatter_2d_diabatic import DiabaticCarrierPINN2D
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
        rec = {'E_c': E_c, 'eta_star': eta,
               'reference_uncertainty_l2': ref_unc, 'runs': []}
        print(f'--- E_c={E_c:.0f} eta*={eta} ref_unc={ref_unc:.4f}', flush=True)
        configs = [('adiabatic+corr300',
                    lambda: CarrierPINN2D(pes, E_c, n_R=400, n_gamma=60,
                                          n_channels=8, n_corr=300,
                                          cap_start=5.0, cap_eta=eta,
                                          cap_width=2.0)),
                   ('diabatic',
                    lambda: DiabaticCarrierPINN2D(pes, E_c, n_R=400,
                                                  n_gamma=60, n_channels=8,
                                                  n_corr=0, cap_start=5.0,
                                                  cap_eta=eta, cap_width=2.0)),
                   ('diabatic+corr300',
                    lambda: DiabaticCarrierPINN2D(pes, E_c, n_R=400,
                                                  n_gamma=60, n_channels=8,
                                                  n_corr=300, cap_start=5.0,
                                                  cap_eta=eta, cap_width=2.0))]
        for label, build in configs:
            t0 = time.time()
            sol = build()
            sv = sol.solve(verbose=False)
            chi = sol.evaluate()
            l2_49 = compare(chi, sol.R, sol.gam, chi_ref, R, g, 4.9)
            l2_69 = compare(chi, sol.R, sol.gam, chi_ref, R, g, 6.9)
            run = dict(label=label,
                       carrier=sv['carrier'], correction=sv['correction'],
                       l2_interior_4p9=l2_49, l2_interior_6p9=l2_69,
                       wall_s=round(time.time() - t0, 1))
            if hasattr(sol, 'dia'):
                run['dia_asymmetry'] = sol.dia.asym
                run['dia_orth_err'] = sol.dia.orth_err
                run['dia_offdiag_ratio'] = sol.dia.offdiag_ratio
                run['column_residuals'] = sol.column_residuals()
                cres = run['column_residuals']
                run['col_res_max'] = max(cres.values())
                run['col_res_open_mean'] = float(np.mean(
                    [v for k, v in cres.items() if k.startswith('dU+_0')
                     or k.startswith('dU+_1') or k.startswith('dU+_2')]))
            rec['runs'].append(run)
            print(f'  {label:>18s}: val={sv["carrier"]["val_residual"]:.3e} '
                  f'L2(4.9)={l2_49:.4f} L2(6.9)={l2_69:.4f} '
                  f'({run["wall_s"]}s)', flush=True)
            if 'col_res_max' in run:
                print(f'      diabatic diagnostics: asym={run["dia_asymmetry"]:.1e} '
                      f'orth={run["dia_orth_err"]:.1e} '
                      f'offdiag_ratio={run["dia_offdiag_ratio"]:.3f} '
                      f'col_res(open mean)={run["col_res_open_mean"]:.3f} '
                      f'max={run["col_res_max"]:.3f}', flush=True)
            np.savez(root / 'results' / f'_scatter2d_diabatic_{int(E_c)}_{label}.npz',
                     R=sol.R, gamma=sol.gam, chi_ref=chi_ref, chi_run=chi)
        out['energies'].append(rec)
        json.dump(out, open(root / 'results' / 'scatter2d_diabatic.json', 'w'),
                  indent=2)
    print('diabatic 基准完成 -> results/scatter2d_diabatic.json', flush=True)


if __name__ == '__main__':
    main()
