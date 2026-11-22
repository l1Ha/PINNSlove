#!/usr/bin/env python3
"""Figure 3: honest, data-driven ablation for the 1D Morse bound states.

All bars are computed LIVE by training the real VariationalPINN with each
architecture setting and benchmarking against the exact analytic Morse
spectrum -- no hard-coded numbers anywhere (AIPro.md SSOT rule).
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import torch

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
torch.set_num_threads(4)

from src.potentials import Morse
from src.hamiltonians import MorseProblem1D
from src.pinn.pinn_bound import VariationalPINN

FIG = ROOT / 'article' / 'fig'


def ground_state_error(hidden, n_fourier, seed, steps=30, lbfgs_steps=60):
    """Train the real VariationalPINN ground state and return |E_PINN - E_exact|."""
    morse = Morse(De_cm1=38293.0, a_A1=1.94286, xe_A=0.74144, mass_amu=0.50397)
    prob = MorseProblem1D(morse, xlo=-1.6, xhi=4.6)
    pinn = VariationalPINN(prob, hidden=hidden, n_quad=300, seed=seed,
                           n_fourier=n_fourier, ortho_weight=3000.0,
                           env_scale=1.5, k_extra=6)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        E, _ = pinn.train_state(0, steps=steps, lr0=1e-3,
                                lbfgs_steps=lbfgs_steps, lbfgs_chunk=25,
                                lbfgs_lr=0.15, patience=3, verbose_every=0)
    # exact analytic Morse ground level
    we = 2.0 * 1.94286 * np.sqrt(38293.0 * morse.B)
    wexe = morse.B * 1.94286 ** 2
    E_exact = -38293.0 + we * 0.5 - wexe * 0.25
    return abs(E - E_exact)


def main():
    print('Running live architecture ablation (real PINN trainings)...')
    configs = [
        ('Plain tanh MLP\n(no Fourier features)', (2, 32), 0, 42),
        ('Plain tanh MLP\n+ 16 Fourier features', (2, 48), 16, 42),
        ('Gated two-stream MLP\n+ 24 Fourier features', (3, 64), 24, 42),
    ]
    errors = []
    for label, hidden, nf, seed in configs:
        err = ground_state_error(hidden, nf, seed)
        errors.append(err)
        print(f'  {label.replace(chr(10), " ")}: ground-state error = {err:.3f} cm^-1')

    fig, ax = plt.subplots(figsize=(7.5, 5.2), dpi=300)
    plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 10})
    colors = ['#aec7e8', '#ffbb78', '#2ca02c']
    bars = ax.bar(range(len(configs)), errors, color=colors, edgecolor='k',
                  lw=1.2, width=0.55)
    for b, val in zip(bars, errors):
        ax.text(b.get_x() + b.get_width() / 2.0, b.get_height() * 1.15,
                rf"${val:.2f}$ cm$^{{-1}}$", ha='center', va='bottom',
                fontsize=9.5, fontweight='bold')
    ax.set_yscale('log')
    ax.set_xticks(range(len(configs)))
    ax.set_xticklabels([c[0] for c in configs], fontsize=9, fontweight='semibold')
    ax.set_ylabel(r'Ground state error $|\Delta E_0|$ (cm$^{-1}$)', fontsize=11,
                  fontweight='bold')
    ax.set_title('(a) Live architecture ablation on the 1D Morse ground state\n'
                 'authentic VariationalPINN trainings, exact analytic benchmark',
                 fontsize=11.5, fontweight='bold')
    ax.grid(True, ls=':', alpha=0.5, which='both')

    out_png = FIG / 'fig_ablation.png'
    fig.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f'[✓] Saved data-driven Figure 3 (ablation): {out_png}')


if __name__ == '__main__':
    main()
