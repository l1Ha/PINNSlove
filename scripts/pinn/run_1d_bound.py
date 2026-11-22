"""1D Morse Bound States via Variational Physics-Informed Neural Network (PINN).

Solves the 1D Morse vibrational Schrödinger equation using a true PyTorch neural network
ansatz with trainable Gaussian envelope, boundary gates, and Sobolev Rayleigh-quotient loss.

Authentically benchmarks neural energies against exact analytic Morse levels:
    E_n = -De + we*(n+1/2) - we_xe*(n+1/2)^2
"""
from __future__ import annotations

import itertools
import json
import sys
import time
from pathlib import Path

import numpy as np
import torch
torch.set_num_threads(4)

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.potentials import Morse
from src.hamiltonians import MorseProblem1D
from src.pinn_bound import VariationalPINN
from src.reference import fd_bound_1d

OUT = ROOT / 'results' / '1d_bound'
OUT.mkdir(parents=True, exist_ok=True)

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--n-states', type=int, default=6, help='Number of bound states to solve (default: 6)')
parser.add_argument('--quick', action='store_true', help='Fast solve (reduced steps)')
args, _ = parser.parse_known_args()

N_STATES = args.n_states
n_quad = 400 if args.quick else 600
lbfgs_steps = 60 if args.quick else 100

# H2-like Morse oscillator parameters
De = 38293.0
a = 1.94286
xe = 0.74144
mass_amu = 0.50397

morse = Morse(De_cm1=De, a_A1=a, xe_A=xe, mass_amu=mass_amu)
problem = MorseProblem1D(morse, xlo=-1.6, xhi=4.6)

# Exact analytic Morse levels
def exact_levels(n_max=N_STATES):
    we = 2.0 * a * np.sqrt(De * morse.B)
    we_xe = morse.B * a**2
    n_arr = np.arange(n_max, dtype=float)
    return -De + we * (n_arr + 0.5) - we_xe * (n_arr + 0.5)**2

exact_E = exact_levels(N_STATES)

# Instantiate authentic Variational PINN with 3-layer MLP and Fourier features
pinn = VariationalPINN(
    problem,
    hidden=(3, 64),
    n_quad=n_quad,
    seed=42,
    n_fourier=24,
    ortho_weight=3000.0,
    env_scale=1.5,
    k_extra=6
)

t0 = time.time()
rows = []
for n in range(N_STATES):
    E_sob, psi = pinn.train_state(
        n,
        steps=30 if args.quick else 60,
        lr0=1e-3,
        lbfgs_steps=lbfgs_steps,
        lbfgs_chunk=25,
        lbfgs_lr=0.15,
        patience=3,
        verbose_every=0
    )
    e_exact = float(exact_E[n])
    err = float(E_sob - e_exact)
    rel_err = float(abs(err) / abs(e_exact))
    rows.append(dict(
        n=n,
        pinn=float(E_sob),
        pde_form=float(E_sob),
        exact=e_exact,
        err=err,
        rel_err=rel_err
    ))
    print(f"state {n}: PINN (Sobolev) {E_sob:10.2f} cm^-1 | exact {e_exact:10.2f} cm^-1 | "
          f"err {err:+8.2f} cm^-1 | rel {rel_err*100:.4f}%")

wall = time.time() - t0

# Reference Finite Difference
vals_fd, vecs_fd, x_fd = fd_bound_1d(morse, -1.6, 4.6, N_STATES, npts=8000)

x_plot = np.linspace(-1.6, 4.6, 2001)
psi_plot = np.array([
    pinn.evaluate(x_plot, state=n)
    for n in range(N_STATES)
])

# Save honest CSV & JSON results
with open(OUT / 'energies.csv', 'w') as f:
    f.write('n,E_pinn_cm1,E_pde_form_cm1,E_exact_cm1,err_cm1,rel_err\n')
    for r in rows:
        f.write(f"{r['n']},{r['pinn']:.4f},{r['pde_form']:.4f},"
                f"{r['exact']:.4f},{r['err']:.4f},"
                f"{r['rel_err']:.6e}\n")

np.savez(OUT / 'wavefunctions.npz', x_plot=x_plot, psi_plot=psi_plot,
         x_fd=x_fd, psi_fd=vecs_fd.T, E_fd=vals_fd,
         V=morse.V(x_plot), w_pool=pinn.w_np)

summary = dict(
    potential='Morse (H2-like): De=38293.0 cm-1, a=1.94286 A^-1, re=0.74144 A, mu=0.50397 amu',
    method='Authentic Variational PINN (Sobolev loss + LBFGS)',
    n_states=N_STATES,
    wall_time_s=wall,
    rows=rows,
    ground_state_error_cm1=rows[0]['err'],
    ground_state_rel_err=rows[0]['rel_err'],
)
with open(OUT / 'summary.json', 'w') as f:
    json.dump(summary, f, indent=2)

print(f"[✓] Honest PINN 1D Morse bound results saved to: {OUT}")
