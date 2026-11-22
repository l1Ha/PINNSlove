"""ELM (random-Fourier-feature linear collocation) benchmark vs the converged
CAP FD reference, with an analytic-operator self-test.

Per docs/theory_scattering_2d.md: the strong-form PINN fails structurally
(122-156%); linear collocation in a fixed feature span is the theory-directed
remedy.  This script:
  1. self-tests the analytic L-phi formula against finite differences
     (guards against silent formula bugs -- a wrong j2 sign would otherwise
     poison every number below);
  2. solves the scattered-wave BVP with ScatteringPINN2DELM for several
     feature counts / frequency scales (ridge chosen by cross-validated PDE
     residual only);
  3. reports honest L2 errors against the converged CAP reference
     (400x60, eta* from the plateau scan), alongside the reference
     uncertainty, and extracts the multichannel Hankel coefficients.
"""
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from src.potentials import WaterPES
from src.pinn.pinn_scatter_2d_elm import ScatteringPINN2DELM
from run_scatter_2d_cap_benchmark import compare, fd

ETA_DEFAULT = {200.0: 100.0, 500.0: 200.0}   # only if the scan json is absent


def self_test_operator(pes, E_c, seed=7):
    """Analytic L[w*phi] vs finite differences on random features."""
    elm = ScatteringPINN2DELM(pes, E_c, n_quad=(50, 10), n_feat=3, seed=seed,
                              cap_start=5.0, cap_eta=100.0, cap_width=2.0)
    rng = np.random.default_rng(seed)
    R = rng.uniform(0.8, 11.0, 7)
    # avoid nodes too close to gamma in {0, pi} (cot singular)
    gam = rng.uniform(0.25, np.pi - 0.25, 7)
    hR, hg = 1e-4, 1e-4

    def phi(Rv, gv):
        th = elm._theta(Rv, gv)
        return np.cos(th) * elm._envelope(Rv)[None, :]   # (M, N)

    d2R_num = (phi(R + hR, gam) - 2 * phi(R, gam) + phi(R - hR, gam)) / hR ** 2
    d2g = (phi(R, gam + hg) - 2 * phi(R, gam) + phi(R, gam - hg)) / hg ** 2
    cot = np.cos(gam) / np.sin(gam)
    dg = (phi(R, gam + hg) - phi(R, gam - hg)) / (2 * hg)
    j2_num = -(d2g + cot * dg)
    V = elm._V_cap(R, gam)
    C = pes.BR / R ** 2 + pes.Br / pes.re ** 2
    L_num = (-pes.BR * d2R_num + C[None, :] * j2_num
             + V[None, :] * phi(R, gam) - elm.E_tot * phi(R, gam))  # (M, N)
    L_ana = elm._apply_L(R, gam)                                     # (N, M)
    rel = float(np.max(np.abs(L_ana.T - L_num)) /
                max(np.max(np.abs(L_num)), 1e-300))
    print(f'operator self-test: max relative deviation = {rel:.2e}')
    assert rel < 5e-5, 'analytic operator formula disagrees with FD'
    return rel


def main():
    pes = WaterPES()
    root = Path(__file__).resolve().parents[2]
    eta_file = root / 'results' / 'scatter2d_cap_benchmark.json'
    eta_star = {}
    if eta_file.exists():
        for rec in json.load(open(eta_file))['energies']:
            eta_star[rec['E_c']] = rec['eta_star']
    out = {'self_test_rel_dev': self_test_operator(pes, 200.0), 'energies': []}
    for E_c in (200.0, 500.0):
        eta = float(eta_star.get(E_c, ETA_DEFAULT[E_c]))
        chi_ref, R, g, K, E_tot = fd(pes, E_c, 400, 60, eta)
        chi_lo, R2, g2, *_ = fd(pes, E_c, 250, 50, eta)
        ref_unc = compare(chi_ref, R, g, chi_lo, R2, g2)
        rec = {'E_c': E_c, 'eta_star': eta, 'reference_uncertainty_l2': ref_unc,
               'runs': []}
        print(f'--- E_c={E_c:.0f}  eta*={eta}  ref uncertainty={ref_unc:.4f}',
              flush=True)
        for n_feat, sigR in [(2000, (8.0, 20.0)), (2000, (8.0, 20.0, 50.0)),
                             (4000, (8.0, 20.0, 50.0, 120.0)),
                             (8000, (20.0, 60.0, 150.0, 300.0))]:
            t0 = time.time()
            elm = ScatteringPINN2DELM(pes, E_c, n_quad=(400, 80),
                                      n_feat=n_feat, sigma_R=sigR, seed=0,
                                      cap_start=5.0, cap_eta=eta,
                                      cap_width=2.0)
            sv = elm.solve(verbose=False)
            ps = elm.evaluate(R, g)
            l2_49 = compare(ps, R, g, chi_ref, R, g, 4.9)
            l2_69 = compare(ps, R, g, chi_ref, R, g, 6.9)
            hk = elm.extract_hankel()
            run = dict(n_feat=n_feat, sigma_R=list(sigR),
                       ridge=sv['best_ridge'],
                       val_pde_residual=sv['best_val_residual'],
                       l2_interior_4p9=l2_49, l2_interior_6p9=l2_69,
                       hankel_flux_ratio=hk['flux_ratio'],
                       hankel_n_open=hk['n_open'],
                       wall_s=round(time.time() - t0, 1))
            rec['runs'].append(run)
            print(f'  M={n_feat} sigma_R={sigR}: L2(4.9)={l2_49:.4f} '
                  f'L2(6.9)={l2_69:.4f} val_res={sv["best_val_residual"]:.2e} '
                  f'flux_ratio={hk["flux_ratio"]:.3f} ({run["wall_s"]}s)',
                  flush=True)
            np.savez(f'results/_scatter2d_elm_{int(E_c)}_{n_feat}.npz',
                     R=R, gamma=g, chi_ref=chi_ref, chi_elm=ps)
        out['energies'].append(rec)
        json.dump(out, open('results/scatter2d_elm.json', 'w'), indent=2)
    print('ELM 基准完成 -> results/scatter2d_elm.json', flush=True)


if __name__ == '__main__':
    main()
