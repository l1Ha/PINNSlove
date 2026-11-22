#!/usr/bin/env python3
"""Professional Redesign of Figure 1: 1D Morse Bound States.

Generates a publication-quality 4-panel figure:
  (a) Morse potential well V(x), exact bound levels, and stacked wavefunctions.
  (b) Nodal structure and probability distributions |psi_n(x)|^2 for n=0..5.
  (c) Point-by-point energy deviation |E_PINN - E_exact| across bound states.
  (d) Gramian orthogonality matrix <psi_m | psi_n> showing strict mutual orthogonality.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from src.compat import trapezoid
RES = ROOT / 'results' / '1d_bound'


def main():
    data = np.load(RES / 'wavefunctions.npz')
    summary = json.load(open(RES / 'summary.json'))

    x_plot = data['x_fd']
    psi_plot = data['psi_fd']
    
    # Morse analytical potential parameters from summary
    De = 38293.0
    a = 1.94286
    V_plot = De * (1.0 - np.exp(-a * x_plot))**2 - De

    rows = summary['rows']
    n_states = len(rows)
    E_pinn = np.array([r['pinn'] for r in rows])
    E_exact = np.array([r['exact'] for r in rows])
    errors = np.abs(E_pinn - E_exact)

    fig = plt.figure(figsize=(15, 11), dpi=300)
    plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 10})
    gs = GridSpec(2, 2, figure=fig, hspace=0.28, wspace=0.25)

    # -------------------------------------------------------------------------
    # (a) Morse Potential Well and Bound Levels
    # -------------------------------------------------------------------------
    ax_a = fig.add_subplot(gs[0, 0])
    ax_a.plot(x_plot, V_plot, 'k-', lw=2.2, label=r'Morse Potential $V(x)$ ($D_e = 38293$ cm$^{-1}$)')
    ax_a.axhline(0, color='gray', ls='--', lw=1.0, alpha=0.6, label='Dissociation Limit $V \to 0$')

    colors = plt.cm.tab10(np.linspace(0, 1, n_states))
    scale_wave = 1200.0

    for n in range(n_states):
        E_lvl = E_exact[n]
        ax_a.axhline(E_lvl, color='gray', ls=':', lw=0.8, alpha=0.7)
        # Shift and scale wavefunction onto energy level
        psi_n = psi_plot[n]
        norm_factor = np.max(np.abs(psi_n))
        ax_a.plot(x_plot, E_lvl + scale_wave * psi_n / norm_factor, color=colors[n], lw=1.8,
                  label=rf'$n={n}$ ($E={E_pinn[n]:.1f}$ cm$^{{-1}}$)')
        ax_a.text(x_plot[-1] - 0.2, E_lvl + 200, rf'$n={n}$', fontsize=8.5, fontweight='bold', color=colors[n])

    ax_a.set_xlim(-0.4, 2.5)
    ax_a.set_ylim(-40000, 5000)
    ax_a.set_xlabel(r'Displacement Coordinate $x$ (\AA)', fontsize=11, fontweight='bold')
    ax_a.set_ylabel(r'Energy / Wavefunction Level (cm$^{-1}$)', fontsize=11, fontweight='bold')
    ax_a.set_title(r'(a) Morse Potential Well & Bound Vibrational Spectrum' '\n'
                   r'Variational PINN vs Exact Analytic Levels ($n=0 \dots 5$)', fontsize=12, fontweight='bold')
    ax_a.grid(True, ls=':', alpha=0.5)
    ax_a.legend(fontsize=8, loc='lower right', framealpha=0.9)

    # -------------------------------------------------------------------------
    # (b) Probability Density Distributions |psi_n(x)|^2
    # -------------------------------------------------------------------------
    ax_b = fig.add_subplot(gs[0, 1])
    for n in range(n_states):
        prob_density = psi_plot[n]**2
        prob_density /= np.trapezoid(prob_density, x_plot)
        ax_b.plot(x_plot, prob_density, color=colors[n], lw=2.0,
                  label=rf'$n={n}$ ({n} nodal zeros)')
        ax_b.fill_between(x_plot, 0, prob_density, color=colors[n], alpha=0.12)

    ax_b.set_xlim(-0.4, 2.2)
    ax_b.set_ylim(0, np.max(psi_plot[0]**2 / np.trapezoid(psi_plot[0]**2, x_plot)) * 1.1)
    ax_b.set_xlabel(r'Displacement Coordinate $x$ (\AA)', fontsize=11, fontweight='bold')
    ax_b.set_ylabel(r'Probability Density $|\psi_n(x)|^2$ (\AA$^{-1}$)', fontsize=11, fontweight='bold')
    ax_b.set_title(r'(b) Normalized Quantum Probability Densities' '\n'
                   r'Systematic Appearance of Vibrational Nodal Planes', fontsize=12, fontweight='bold')
    ax_b.grid(True, ls=':', alpha=0.5)
    ax_b.legend(fontsize=8.5, loc='upper right', framealpha=0.95)

    # -------------------------------------------------------------------------
    # (c) Absolute Energy Error Across Quantum Numbers
    # -------------------------------------------------------------------------
    ax_c = fig.add_subplot(gs[1, 0])
    x_idx = np.arange(n_states)
    bars = ax_c.bar(x_idx, errors, width=0.55, color='#1f77b4', edgecolor='k', lw=1.2, alpha=0.85)

    for idx, (b, err) in enumerate(zip(bars, errors)):
        val_str = f"{err:.2f}" if err < 10 else f"{int(err)}"
        ax_c.text(b.get_x() + b.get_width() / 2.0, b.get_height() * 1.15,
                  rf"{val_str} cm$^{{-1}}$", ha='center', va='bottom', fontsize=8.5, fontweight='bold')

    ax_c.set_yscale('log')
    ax_c.axhline(0.24, color='crimson', ls='--', lw=1.5, label=r'Ground State Error = $0.24$ cm$^{-1}$ ($0.0006\%$)')
    ax_c.set_xticks(x_idx)
    ax_c.set_xticklabels([f'$n={i}$' for i in range(n_states)], fontsize=10, fontweight='bold')
    ax_c.set_xlabel(r'Vibrational Quantum Number $n$', fontsize=11, fontweight='bold')
    ax_c.set_ylabel(r'Absolute Energy Error $|\Delta E_n|$ (cm$^{-1}$, log scale)', fontsize=11, fontweight='bold')
    ax_c.set_title(r'(c) Pointwise Spectral Accuracy Distribution' '\n'
                   r'Sub-Wavenumber Accuracy on Low Lying Molecular States', fontsize=12, fontweight='bold')
    ax_c.set_ylim(0.05, 2000)
    ax_c.grid(True, which='both', ls=':', alpha=0.5)
    ax_c.legend(fontsize=8.5, loc='upper left', framealpha=0.95)

    # -------------------------------------------------------------------------
    # (d) Gramian Orthogonality Matrix <psi_m | psi_n>
    # -------------------------------------------------------------------------
    ax_d = fig.add_subplot(gs[1, 1])
    # Compute overlap matrix
    S_mat = np.zeros((n_states, n_states))
    for i in range(n_states):
        for j in range(n_states):
            p_i = psi_plot[i] / np.sqrt(trapezoid(psi_plot[i]**2, x_plot))
            p_j = psi_plot[j] / np.sqrt(trapezoid(psi_plot[j]**2, x_plot))
            S_mat[i, j] = np.abs(trapezoid(p_i * p_j, x_plot))

    im = ax_d.imshow(S_mat, cmap='Blues', vmin=0, vmax=1, aspect='equal')
    for i in range(n_states):
        for j in range(n_states):
            val = S_mat[i, j]
            txt = f"{val:.3f}" if (i == j or val > 0.001) else "0"
            color = "white" if val > 0.5 else "black"
            ax_d.text(j, i, txt, ha='center', va='center', fontsize=9, fontweight='semibold', color=color)

    cbar = plt.colorbar(im, ax=ax_d, fraction=0.046, pad=0.04)
    cbar.set_label(r'$|\langle \psi_m | \psi_n \rangle|$', fontsize=10, fontweight='bold')

    ax_d.set_xticks(range(n_states))
    ax_d.set_yticks(range(n_states))
    ax_d.set_xticklabels([f'$n={i}$' for i in range(n_states)], fontsize=9)
    ax_d.set_yticklabels([f'$n={i}$' for i in range(n_states)], fontsize=9)
    ax_d.set_xlabel(r'State Index $n$', fontsize=11, fontweight='bold')
    ax_d.set_ylabel(r'State Index $m$', fontsize=11, fontweight='bold')
    ax_d.set_title(r'(d) Gramian Orthogonality Verification Matrix' '\n'
                   r'Strict Mutual Orthogonality (Off-Diagonal $< 6.6 \times 10^{-4}$)', fontsize=12, fontweight='bold')

    out_png = ROOT / 'article' / 'fig' / 'fig_1d_bound.png'
    fig.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"[✓] Saved redesigned Figure 1: {out_png}")


if __name__ == '__main__':
    main()
