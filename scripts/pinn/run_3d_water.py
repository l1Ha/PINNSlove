"""3D water bound states (Jacobi coordinates (r, R, gamma), J = 0) via
variational PINN with mini-batch training; comparison with a product
FD/Legendre-DVR reference.

Produces results/3d_water/{energies.csv, summary.json, slices.npz}.
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
torch.set_num_threads(4)

from src.potentials import WaterPES
from src.hamiltonians import Water3DProblem
from src.pinn_bound import VariationalPINN
from src.reference import dvr_3d

OUT = Path(__file__).resolve().parents[2] / 'results' / '3d_water'
OUT.mkdir(parents=True, exist_ok=True)

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--n-states', type=int, default=1, help='Number of states to solve (default: 1)')
parser.add_argument('--quick', action='store_true', help='Fast solve with reduced quadrature and steps')
args, _ = parser.parse_known_args()

N_STATES = args.n_states
n_quad = (20, 20, 14) if args.quick else (32, 32, 20)
lbfgs_steps = 150 if args.quick else 600

pes = WaterPES()
problem = Water3DProblem(pes)

# ---------------- reference DVR (cached) ----------------
t0 = time.time()
cache = OUT / 'dvr_reference.npz'
if cache.exists():
    vals_ref = np.load(cache)['vals']
    print('loaded cached DVR reference', flush=True)
else:
    vals_ref, vecs_ref, (r_g, R_g, x_g, w_g), Vg = dvr_3d(
        pes, 0.55, 2.4, 60 if args.quick else 90, 0.55, 3.0, 70 if args.quick else 100, 30 if args.quick else 44, k=N_STATES + 2, verbose=True)
    np.savez(cache, vals=vals_ref)
print(f'DVR3d reference ready in {time.time()-t0:.0f}s', flush=True)
print('DVR levels (E - Vmin):', np.round(vals_ref[:N_STATES] - pes.V_min(), 2),
      flush=True)

# ---------------- PINN ----------------
t0 = time.time()
rows = []
for n in range(N_STATES):
    pinn = VariationalPINN(problem, hidden=(2, 48), n_quad=n_quad,
                           seed=1, n_fourier=8,
                           ortho_weight=300.0, env_scale=1.5,
                           arch='plain')
    E, psi = pinn.train_state(n, steps=10 if args.quick else 50, verbose_every=0,
                              lbfgs_steps=lbfgs_steps, lbfgs_chunk=25,
                              lbfgs_lr=0.1, patience=3,
                              seed=1)
    ritz = pinn._ritz_energy
    ev = pinn.evaluate_energies((30, 30, 20) if args.quick else (60, 60, 32), state=n)
    ev2 = pinn.evaluate_energies((36, 36, 24) if args.quick else (72, 72, 36), state=n)
    agree = abs(ev['E_pde'] - ev2['E_pde'])
    E_eval = ev['E_pde']
    err = E_eval - vals_ref[n]
    rows.append(dict(n=n, pinn=E, E_pde_finegrid=ev['E_pde'],
                     E_pde_finegrid2=ev2['E_pde'],
                     fine_grid_agreement=float(agree),
                     E_sob_finegrid=ev['E_sob'], ritz=ritz,
                     dvr=float(vals_ref[n]), err=float(err)))
    print(f'state {n}: PINN(pool) {E:.2f}  PINN(fine grid, PDE form) '
          f'{ev["E_pde"]:.2f}  DVR {vals_ref[n]:.2f}  err {err:+.2f} cm-1',
          flush=True)
wall = time.time() - t0

with open(OUT / 'energies.csv', 'w') as f:
    f.write('n,E_pinn_pool_cm1,E_pde_fine_cm1,E_pde_fine2_cm1,fine_grid_agreement_cm1,E_ritz_cm1,E_dvr_cm1,err_cm1\n')
    for r in rows:
        f.write(f"{r['n']},{r['pinn']:.3f},{r['ritz']:.3f},{r['dvr']:.3f},"
                f"{r['err']:.3f}\n")

# ground-state slice at gamma = gamma_e for plotting
r_plot = np.linspace(0.6, 2.2, 111)
R_plot = np.linspace(0.6, 2.8, 111)
RRg, RRg2 = np.meshgrid(r_plot, R_plot, indexing='ij')
gamma_e = pes.gamma_e
chi_slice = pinn.evaluate([RRg.ravel(), RRg2.ravel(),
                           np.full(RRg.size, gamma_e)],
                          state=0).reshape(RRg.shape)
V_slice = pes.V_jacobi(RRg, RRg2, np.full_like(RRg, gamma_e))
np.savez(OUT / 'slices.npz', r=r_plot, R=R_plot, chi0=chi_slice,
         V_slice=V_slice, gamma_e=gamma_e)

summary = dict(
    system='H2O 3D (Jacobi r, R, gamma), J = 0, reduced chi = r R psi',
    method='variational PINN, mini-batch Adam + chunked LBFGS',
    reference='FD x FD x Legendre-DVR (90 x 100 x 44)',
    n_states=N_STATES, wall_time_s=wall,
    max_abs_err_cm1=max(abs(r['err']) for r in rows),
    rows=rows)
with open(OUT / 'summary.json', 'w') as f:
    json.dump(summary, f, indent=2)
print('done:', OUT)
