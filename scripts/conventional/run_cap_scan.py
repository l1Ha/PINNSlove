"""CAP strength scan + box/grid convergence study for the 2D FD reference.

The uncapped deep-well box reference is NOT convergent (see
run_scatter2d_convergence.py: grid-to-grid interior L2 ~ 130-160%, box-size
interior change 83-136%).  A quadratic CAP absorber -i*eta*s^2 (R > cap_start)
removes the box resonance.  cap_eta is selected by the PLATEAU criterion: the
interior solution (R < cap_start) must be stationary under both box-size
change (R_hi = 9 -> 12) and grid refinement at fixed cap_start/cap_width.
This is a numerical-regularization convergence study; eta is not fitted to
any target answer.

Output: results/cap_scan.json + stdout table.
"""
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from src.potentials import WaterPES
from src.conventional.scatter_2d import scatter_fd_2d

ETAS = [0, 100, 300, 1000, 3000, 10000]
E_LIST = [200.0, 500.0]
CAP_START, CAP_WIDTH = 5.0, 2.0
EDGE = 0.05


def interior(a, edge=EDGE):
    n = len(a)
    return a[int(edge * n):int((1 - edge) * n)]


def compare(chi_a, R_a, g_a, chi_b, R_b, g_b, R_max_interior=4.9):
    """Interior (R < R_max_interior, 5% edge strips removed) L2 difference,
    both solutions interpolated onto a common grid."""
    from scipy.interpolate import RegularGridInterpolator
    R_lo = max(R_a[0], R_b[0])
    R_hi_cm = min(R_a[-1], R_b[-1], R_max_interior)
    m_a = (R_a >= R_lo) & (R_a <= R_hi_cm)
    R_sel = interior(R_a[m_a])
    g_sel = interior(g_a)
    RR, GG = np.meshgrid(R_sel, g_sel, indexing='ij')
    pts = np.stack([RR.ravel(), GG.ravel()], axis=-1)

    def ip(chi, R, g):
        out = np.empty(RR.shape, dtype=complex)
        for i, part in enumerate((chi.real, chi.imag)):
            f = RegularGridInterpolator((R, g), part, bounds_error=False,
                                        fill_value=None)
            r = f(pts).reshape(RR.shape)
            if i == 0:
                out.real = r
            else:
                out.imag = r
        return out

    a, b = ip(chi_a, R_a, g_a), ip(chi_b, R_b, g_b)
    return float(np.sqrt(np.sum(np.abs(a - b) ** 2))
                 / max(np.sqrt(np.sum(np.abs(a) ** 2)), 1e-300))


def main():
    pes = WaterPES()
    rows = []
    for E_c in E_LIST:
        for eta in ETAS:
            kw = dict(n_R=150, n_gamma=30, cap_start=CAP_START,
                      cap_width=CAP_WIDTH)
            if eta:
                kw['cap_eta'] = eta
            t0 = time.time()
            chi, R, g, K, E_tot = scatter_fd_2d(pes, E_c, **kw)
            chiH, RH, gH, *_ = scatter_fd_2d(pes, E_c, R_hi=12.0, **{
                k: v for k, v in kw.items() if k != 'R_hi'} | {})
            chiG, RG, gG, *_ = scatter_fd_2d(pes, E_c, n_R=250, n_gamma=50,
                                             R_hi=9.0, cap_start=CAP_START,
                                             cap_width=CAP_WIDTH,
                                             cap_eta=eta or None)
            rec = dict(E_c=E_c, eta=eta,
                       norm=float(np.sqrt(np.sum(np.abs(chi) ** 2))),
                       box_l2=compare(chi, R, g, chiH, RH, gH),
                       grid_l2=compare(chi, R, g, chiG, RG, gG),
                       wall_s=round(time.time() - t0, 1))
            rows.append(rec)
            print(f"E_c={E_c:.0f}  eta={eta:>6}: norm={rec['norm']:9.3f}  "
                  f"box L2(9 vs 12)={rec['box_l2']:.4f}  "
                  f"grid L2(150x30 vs 250x50)={rec['grid_l2']:.4f}  "
                  f"({rec['wall_s']}s)", flush=True)
    Path('results').mkdir(exist_ok=True)
    json.dump(rows, open('results/cap_scan.json', 'w'), indent=2)
    print('CAP 扫描完成 -> results/cap_scan.json', flush=True)


if __name__ == '__main__':
    main()
