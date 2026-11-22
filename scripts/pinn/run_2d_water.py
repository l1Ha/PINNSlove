"""2D water bound states (frozen OH bond, coordinates (R, gamma)) via
variational PINN; comparison with a finite-difference/DVR reference.

Produces results/2d_water/{energies.csv, summary.json, wavefunctions.npz}.
"""
import itertools
import json
import sys
import time
from pathlib import Path

import numpy as np
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
torch.set_num_threads(2)

from src.potentials import WaterPES
from src.hamiltonians import Water2DProblem
from src.pinn_bound import VariationalPINN
from src.reference import dvr_2d

OUT = Path(__file__).resolve().parents[2] / 'results' / '2d_water'
OUT.mkdir(parents=True, exist_ok=True)

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--n-states', type=int, default=3, help='Number of bound states to solve (default: 3)')
parser.add_argument('--quick', action='store_true', help='Fast solve with reduced quadrature and steps')
args, _ = parser.parse_known_args()

N_STATES = args.n_states
n_quad = (64, 48) if args.quick else (96, 72)
lbfgs_steps = 300 if args.quick else 800

pes = WaterPES()
problem = Water2DProblem(pes)

# ---------------- reference DVR ----------------
t0 = time.time()
vals_ref, vecs_ref, R_ref, x_ref, w_ref, Vgrid = dvr_2d(
    pes, pes.re, 0.55, 3.2, 350 if args.quick else 500, 60 if args.quick else 80, k=N_STATES + 2)
print(f'DVR reference done in {time.time()-t0:.0f}s', flush=True)
print('DVR levels (E - Vmin):', np.round(vals_ref[:N_STATES] - pes.V_min(), 2),
      flush=True)

# ---------------- PINN ----------------
pinn = VariationalPINN(problem, hidden=(2, 48), n_quad=n_quad, seed=1,
                       n_fourier=16, ortho_weight=500.0, env_scale=1.5,
                       basis_scales=(1.0, 1.8), k_extra=6)

t0 = time.time()
rows = []
for n in range(N_STATES):
    E, psi = pinn.train_state(n, steps=20 if args.quick else 50, verbose_every=0,
                              lbfgs_steps=lbfgs_steps, seed=1,
                              lbfgs_chunk=30, lbfgs_lr=0.1, patience=4)
    ritz = pinn._ritz_energy
    ev = pinn.evaluate_energies(n_quad, state=n)
    ev2 = pinn.evaluate_energies((80, 60) if args.quick else (128, 96), state=n)
    agree = abs(ev['E_pde'] - ev2['E_pde'])
    err = ev['E_pde'] - vals_ref[n]
    rows.append(dict(n=n, pinn=E, E_pde_finegrid=ev['E_pde'],
                     E_pde_finegrid2=ev2['E_pde'], fine_grid_agreement=float(agree),
                     E_sob_finegrid=ev['E_sob'], ritz=ritz,
                     dvr=float(vals_ref[n]), err=float(err)))
    print(f'state {n}: PINN(pool) {E:.2f}  PINN(fine, PDE) {ev["E_pde"]:.2f}  '
          f'Ritz {ritz:.2f}  DVR {vals_ref[n]:.2f}  err {err:+.2f} cm-1',
          flush=True)
wall = time.time() - t0

ovs = [float(abs(np.sum(pinn.w_np * pinn.states[i]['psi']
                        * pinn.states[j]['psi'])))
       for i, j in itertools.combinations(range(N_STATES), 2)]

# wavefunction snapshots on a polar grid for plotting
R_plot = np.linspace(0.6, 3.0, 141)
gam_plot = np.linspace(0.05, np.pi - 0.05, 121)
RR, GG = np.meshgrid(R_plot, gam_plot, indexing='ij')
psi_grid = np.stack([
    pinn.evaluate([RR.ravel(), GG.ravel()], state=n).reshape(RR.shape)
    for n in range(N_STATES)])
V_plot = pes.V_jacobi(pes.re, RR, GG)

with open(OUT / 'energies.csv', 'w') as f:
    f.write('n,E_pinn_pool_cm1,E_pde_fine_cm1,E_pde_fine2_cm1,fine_grid_agreement_cm1,E_ritz_cm1,E_dvr_cm1,err_cm1\n')
    for r in rows:
        f.write(f"{r['n']},{r['pinn']:.3f},{r['E_pde_finegrid']:.3f},"
                f"{r['E_pde_finegrid2']:.3f},{r['fine_grid_agreement']:.3f},"
                f"{r['ritz']:.3f},{r['dvr']:.3f},{r['err']:.3f}\n")

np.savez(OUT / 'wavefunctions.npz', R=R_plot, gamma=gam_plot,
         psi=psi_grid, V=V_plot, R_ref=R_ref, x_ref=x_ref,
         psi_ref=vecs_ref[:, :N_STATES].T.reshape(N_STATES, len(R_ref), len(x_ref)),
         E_ref=vals_ref[:N_STATES])

summary = dict(
    system='H2O 2D frozen-bond model (r = re), Jacobi (R, gamma)',
    method='variational PINN (variance loss + trainable Gaussian envelope)',
    reference='FD x Legendre-DVR (500 x 80)',
    n_states=N_STATES, wall_time_s=wall,
    max_abs_err_cm1=max(abs(r['err']) for r in rows),
    max_state_overlap=max(ovs),
    rows=rows)
with open(OUT / 'summary.json', 'w') as f:
    json.dump(summary, f, indent=2)
print('done:', OUT)
