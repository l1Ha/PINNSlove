#!/usr/bin/env python3
"""Professional Redesign of Figure 4: 1D Morse Fixed-Energy Scattering.

Generates a publication-quality multi-panel figure:
  (a) Morse scattering potential V(x), collision energies, and de Broglie wavelength compression.
  (b) Genuine PINN vs Semiclassical Langer vs FD wavefunctions at E = 1000 cm^-1.
  (c) Wavefunction matching at E = 2000 cm^-1 across turning point and asymptotic outgoing region.
  (d) Wavefunction matching at E = 4000 cm^-1.
  (e) Quantitative comparison: Genuine PINN phase-slip behavior vs First-Principles Semiclassical Langer.
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
RES = ROOT / 'results' / '1d_scatter'
FIG = ROOT / 'article' / 'fig'


def main():
    d = np.load(RES / 'wavefunctions.npz')
    x = d['x']
    Es = d['Es']
    psi_pinn = d['psi_pinn']
    psi_sc = d['psi_sc']
    psi_ref = d['psi_ref']
    V = d['V']
    summary_data = json.load(open(RES / 'summary.json'))
    summary = summary_data.get('results', summary_data)

    fig = plt.figure(figsize=(15, 11), dpi=300)
    plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 10})
    gs = GridSpec(3, 2, figure=fig, hspace=0.35, wspace=0.25)

    sub = slice(0, len(x), 4)
    x_sub = x[sub]

    # (a) Potential Well & Physical Regimes
    ax_a = fig.add_subplot(gs[0, :])
    ax_a.plot(x_sub, V[sub], 'k-', lw=2.2, label=r'Morse Potential $V(x)$ ($D_e = 38293$ cm$^{-1}$)')
    ax_a.axhline(0, color='gray', ls='--', lw=1.0, alpha=0.6, label='Asymptotic Continuum Threshold $V \to 0$')

    colors_E = ['#1f77b4', '#2ca02c', '#d62728']
    for idx_e, E_val in enumerate(Es):
        ax_a.axhline(E_val, color=colors_E[idx_e], ls=':', lw=1.5,
                     label=rf'$E = {int(E_val)}$ cm$^{{-1}}$ ($k = {summary[idx_e]["k"]:.2f}$ \AA$^{{-1}}$)')

    ax_a.axvspan(-0.3, 0.4, color='orange', alpha=0.12, label='Attractive Well: $6.3\\times$ de Broglie Compression')
    ax_a.axvspan(0.4, 2.8, color='lightblue', alpha=0.12, label='Asymptotic Scattering Region (Robin BC)')

    ax_a.set_xlim(-0.4, 2.5)
    ax_a.set_ylim(-40000, 8000)
    ax_a.set_xlabel(r'Scattering Coordinate $x$ (\AA)', fontsize=11, fontweight='bold')
    ax_a.set_ylabel(r'Energy $V(x)$ (cm$^{-1}$)', fontsize=11, fontweight='bold')
    ax_a.set_title(r'(a) 1D Morse Scattering Potential and Continuum Energy Levels' '\n'
                   r'Deep Attractive Well Inducing Severe Interior Spatial Wavelength Compression', fontsize=12, fontweight='bold')
    ax_a.grid(True, ls=':', alpha=0.5)
    ax_a.legend(fontsize=8.5, loc='lower right', ncol=3, framealpha=0.95)

    # (b) Wavefunction at E = 1000 cm^-1
    ax_b = fig.add_subplot(gs[1, 0])
    p_ref_0 = psi_ref[0][sub]
    p_lgr_0 = psi_sc[0][sub]
    p_pnn_0 = psi_pinn[0][sub]

    ax_b.plot(x_sub, p_ref_0.real, 'k-', lw=1.8, label=r'Re $\psi_{\mathrm{FD}}$ (Reference)')
    ax_b.plot(x_sub, p_lgr_0.real, 'g--', lw=1.8, label=rf'Re $\psi_{{\mathrm{{carrier}}}}$ ($L^2 = {summary[0]["l2_sc"]*100:.2f}\%$)')
    ax_b.plot(x_sub, p_pnn_0.real, 'r:', lw=1.8, label=rf'Re $\psi_{{\mathrm{{PINN}}}}$ ($L^2 = {summary[0]["l2_pinn"]*100:.1f}\%$)')

    ax_b.set_xlim(-0.3, 2.0)
    ax_b.set_ylabel(r'Wavefunction $\psi(x)$', fontsize=10, fontweight='bold')
    ax_b.set_title(rf'(b) Collision Energy $E = 1000$ cm$^{{-1}}$: First-Principles Semiclassical vs Neural Solution' '\n'
                   rf'Langer $|R| = {summary[0]["absR_sc"]:.4f}$, PINN $|R| = {summary[0]["absR_pinn"]:.4f}$', fontsize=11, fontweight='bold')
    ax_b.grid(True, ls=':', alpha=0.5)
    ax_b.legend(fontsize=8, loc='upper left', framealpha=0.9)

    # (c) Wavefunction at E = 2000 cm^-1
    ax_c = fig.add_subplot(gs[1, 1])
    p_ref_1 = psi_ref[1][sub]
    p_lgr_1 = psi_sc[1][sub]
    p_pnn_1 = psi_pinn[1][sub]

    ax_c.plot(x_sub, p_ref_1.real, 'k-', lw=1.8, label=r'Re $\psi_{\mathrm{FD}}$')
    ax_c.plot(x_sub, p_lgr_1.real, 'g--', lw=1.8, label=rf'Re $\psi_{{\mathrm{{carrier}}}}$ ($L^2 = {summary[1]["l2_sc"]*100:.2f}\%$)')
    ax_c.plot(x_sub, p_pnn_1.real, 'r:', lw=1.8, label=rf'Re $\psi_{{\mathrm{{PINN}}}}$ ($L^2 = {summary[1]["l2_pinn"]*100:.1f}\%$)')

    ax_c.set_xlim(-0.3, 2.0)
    ax_c.set_title(rf'(c) Collision Energy $E = 2000$ cm$^{{-1}}$: Wavefunction Profile Comparison' '\n'
                   rf'Langer $|R| = {summary[1]["absR_sc"]:.4f}$, PINN $|R| = {summary[1]["absR_pinn"]:.4f}$', fontsize=11, fontweight='bold')
    ax_c.grid(True, ls=':', alpha=0.5)
    ax_c.legend(fontsize=8, loc='upper left', framealpha=0.9)

    # (d) Wavefunction at E = 4000 cm^-1
    ax_d = fig.add_subplot(gs[2, 0])
    p_ref_2 = psi_ref[2][sub]
    p_lgr_2 = psi_sc[2][sub]
    p_pnn_2 = psi_pinn[2][sub]

    ax_d.plot(x_sub, p_ref_2.real, 'k-', lw=1.8, label=r'Re $\psi_{\mathrm{FD}}$')
    ax_d.plot(x_sub, p_lgr_2.real, 'g--', lw=1.8, label=rf'Re $\psi_{{\mathrm{{carrier}}}}$ ($L^2 = {summary[2]["l2_sc"]*100:.2f}\%$)')
    ax_d.plot(x_sub, p_pnn_2.real, 'r:', lw=1.8, label=rf'Re $\psi_{{\mathrm{{PINN}}}}$ ($L^2 = {summary[2]["l2_pinn"]*100:.1f}\%$)')

    ax_d.set_xlim(-0.3, 2.0)
    ax_d.set_xlabel(r'Scattering Coordinate $x$ (\AA)', fontsize=10, fontweight='bold')
    ax_d.set_ylabel(r'Wavefunction $\psi(x)$', fontsize=10, fontweight='bold')
    ax_d.set_title(rf'(d) Collision Energy $E = 4000$ cm$^{{-1}}$: High-Frequency Regime' '\n'
                   rf'Langer $|R| = {summary[2]["absR_sc"]:.4f}$, PINN $|R| = {summary[2]["absR_pinn"]:.4f}$', fontsize=11, fontweight='bold')
    ax_d.grid(True, ls=':', alpha=0.5)
    ax_d.legend(fontsize=8, loc='upper left', framealpha=0.9)

    # (e) Quantitative Error Comparison
    ax_e = fig.add_subplot(gs[2, 1])
    x_bar = np.arange(len(Es))
    w_bar = 0.35

    l2_pnn = [s['l2_pinn'] * 100 for s in summary]
    l2_lgr = [s['l2_sc'] * 100 for s in summary]

    b1 = ax_e.bar(x_bar - w_bar/2, l2_pnn, width=w_bar, color='#d62728', edgecolor='k', lw=1.0,
                  label=r'Physics-informed carrier + cross-validated neural correction')
    b2 = ax_e.bar(x_bar + w_bar/2, l2_lgr, width=w_bar, color='#2ca02c', edgecolor='k', lw=1.0,
                  label=r'Pure semiclassical Langer-Airy ($L^2 \approx 1.4-1.8\%$)')

    for b, val in zip(b1, l2_pnn):
        ax_e.text(b.get_x() + b.get_width() / 2.0, b.get_height() + 2.0,
                  rf"${val:.1f}\%$", ha='center', va='bottom', fontsize=8.5, fontweight='bold', color='crimson')
    for b, val in zip(b2, l2_lgr):
        ax_e.text(b.get_x() + b.get_width() / 2.0, b.get_height() + 2.0,
                  rf"${val:.2f}\%$", ha='center', va='bottom', fontsize=8.5, fontweight='bold', color='darkgreen')

    ax_e.axhline(1.0, color='k', ls='--', lw=1.2, label=r'1% Target Error Reference')
    ax_e.set_xticks(x_bar)
    ax_e.set_xticklabels([f'$E = {int(E)}$ cm$^{{-1}}$' for E in Es], fontsize=9.5, fontweight='semibold')
    ax_e.set_xlabel(r'Collision Energy', fontsize=10, fontweight='bold')
    ax_e.set_ylabel(r'Relative $L^2$ Scattering Error ($\%$)', fontsize=10, fontweight='bold')
    ax_e.set_ylim(0, 3)
    ax_e.set_title(r'(e) S-Matrix Scattering Accuracy Benchmark' '\n'
                   r'Carrier + Cross-Validated Correction vs Pure Semiclassical', fontsize=11, fontweight='bold')
    ax_e.grid(True, ls=':', alpha=0.5)
    ax_e.legend(fontsize=8, loc='upper right', framealpha=0.95)

    out_png = FIG / 'fig_1d_scatter.png'
    fig.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"[✓] Saved updated Figure 4: {out_png}")


if __name__ == '__main__':
    main()
