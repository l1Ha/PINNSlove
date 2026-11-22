#!/usr/bin/env python3
"""Professional Redesign of Figure 2: 2D and 3D Water Vibrational Wavefunctions and PES.

Generates:
  1. article/fig/fig_2d_states.png: 2x2 grid of frozen-bond 2D wavefunctions (states 0..3)
     with clear nodal zero contours, energy labels, and professional typography.
  2. article/fig/fig_3d_slice.png: 1x2 side-by-side high-resolution comparison of
     3D ground state chi(r, R) at gamma=gamma_e vs the authentic PES cut V(r, R).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / 'results'
FIG = ROOT / 'article' / 'fig'


def plot_2d_states():
    d = np.load(RES / '2d_water' / 'wavefunctions.npz')
    R = d['R']
    gam_deg = np.rad2deg(d['gamma'])
    psi = d['psi']
    summary = json.load(open(RES / '2d_water' / 'summary.json'))

    rows = summary['rows']
    E_pinn = [r['pinn'] for r in rows]
    E_dvr = [r['dvr'] for r in rows]

    # Create 2x2 panel figure
    fig, axes = plt.subplots(2, 2, figsize=(10, 8.5), dpi=300)
    plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 10})

    state_titles = [
        r'(a) State 0 (Ground State, $\nu_0$)',
        r'(b) State 1 (Bending Mode $\nu_2$, $+1837$ cm$^{-1}$)',
        r'(c) State 2 (Symmetric Stretch $\nu_1$, $+3802$ cm$^{-1}$)',
        r'(d) State 3 (Asymmetric Stretch $\nu_3$, $+4734$ cm$^{-1}$)',
    ]

    for idx, ax in enumerate(axes.flat):
        P = psi[idx]
        vmax = np.max(np.abs(P))
        cf = ax.contourf(gam_deg, R, P, levels=30, cmap='RdBu_r', vmin=-vmax, vmax=vmax)
        # Nodal line
        cs = ax.contour(gam_deg, R, P, levels=[0.0], colors='k', linewidths=1.2, linestyles='--')

        # Equilibrium point (gamma_e ~ 104.5 deg, R_e ~ 1.43 A)
        ax.plot([104.5], [1.43], 'm*', ms=10, mew=1.2, mfc='yellow', label=r'Equilibrium $(\gamma_e, R_e)$' if idx == 0 else "")

        ax.set_title(f"{state_titles[idx]}\n"
                     rf"$E_{{\mathrm{{PINN}}}} = {E_pinn[idx]:.1f}$ cm$^{{-1}}$ (DVR: ${E_dvr[idx]:.1f}$ cm$^{{-1}}$)",
                     fontsize=10, fontweight='bold')
        ax.set_xlim(50, 160)
        ax.set_ylim(0.8, 2.5)
        ax.grid(True, ls=':', alpha=0.4)

        if idx >= 2:
            ax.set_xlabel(r'Bond Angle $\gamma$ (deg)', fontsize=10, fontweight='bold')
        if idx % 2 == 0:
            ax.set_ylabel(r'Jacobi Distance $R$ (\AA)', fontsize=10, fontweight='bold')

        cbar = fig.colorbar(cf, ax=ax, fraction=0.046, pad=0.03)
        cbar.ax.tick_params(labelsize=8)

    axes[0, 0].legend(fontsize=8, loc='upper left', framealpha=0.9)
    plt.tight_layout()

    out_png = FIG / 'fig_2d_states.png'
    fig.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"[✓] Saved redesigned 2D states figure: {out_png}")


def plot_3d_slice():
    d = np.load(RES / '3d_water' / 'slices.npz')
    r = d['r']
    R = d['R']
    chi0 = d['chi0']
    V_slice = d['V_slice']
    gam_e = np.rad2deg(float(d['gamma_e']))

    fig, axes = plt.subplots(1, 2, figsize=(11, 5.0), dpi=300)
    plt.rcParams.update({'font.family': 'sans-serif', 'font.size': 10})

    # Left: Wavefunction chi0(r, R)
    ax0 = axes[0]
    cf0 = ax0.contourf(r, R, chi0.T, levels=30, cmap='viridis')
    ax0.plot([0.958], [1.427], 'r*', ms=11, mew=1.2, mfc='yellow', label=r'Minimum $(r_e, R_e) = (0.96, 1.43)$ \AA')
    ax0.set_title(rf'(a) Ground State $\chi_0(r, R)$ at $\gamma = {gam_e:.1f}^\circ$' '\n'
                  r'Localized Wave Packet ($ZPE = 4642$ cm$^{-1}$)', fontsize=11, fontweight='bold')
    ax0.set_xlabel(r'Bond Length $r_{\mathrm{OH}}$ (\AA)', fontsize=10, fontweight='bold')
    ax0.set_ylabel(r'Jacobi Distance $R$ (\AA)', fontsize=10, fontweight='bold')
    ax0.set_xlim(0.65, 1.65)
    ax0.set_ylim(0.9, 2.3)
    ax0.grid(True, ls=':', alpha=0.4)
    ax0.legend(fontsize=8.5, loc='upper right', framealpha=0.9)
    cb0 = fig.colorbar(cf0, ax=ax0, fraction=0.046, pad=0.04)
    cb0.set_label(r'$\chi_0$ (arb. units)', fontsize=9, fontweight='bold')

    # Right: Potential energy surface V(r, R)
    ax1 = axes[1]
    levels_V = np.linspace(-81000, -35000, 32)
    cf1 = ax1.contourf(r, R, V_slice.T, levels=levels_V, cmap='inferno')
    cs1 = ax1.contour(r, R, V_slice.T, levels=np.linspace(-80000, -50000, 7), colors='white', linewidths=0.7, alpha=0.6)
    ax1.plot([0.958], [1.427], 'r*', ms=11, mew=1.2, mfc='yellow', label=r'Potential Well Bottom')
    ax1.set_title(rf'(b) Water Potential $V(r, R)$ at $\gamma = {gam_e:.1f}^\circ$' '\n'
                  r'Anharmonic Morse-Coupled Well Geometry', fontsize=11, fontweight='bold')
    ax1.set_xlabel(r'Bond Length $r_{\mathrm{OH}}$ (\AA)', fontsize=10, fontweight='bold')
    ax1.set_ylabel(r'Jacobi Distance $R$ (\AA)', fontsize=10, fontweight='bold')
    ax1.set_xlim(0.65, 1.65)
    ax1.set_ylim(0.9, 2.3)
    ax1.grid(True, ls=':', alpha=0.4)
    ax1.legend(fontsize=8.5, loc='upper right', framealpha=0.9)
    cb1 = fig.colorbar(cf1, ax=ax1, fraction=0.046, pad=0.04)
    cb1.set_label(r'$V$ (cm$^{-1}$)', fontsize=9, fontweight='bold')

    plt.tight_layout()
    out_png = FIG / 'fig_3d_slice.png'
    fig.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"[✓] Saved redesigned 3D slice figure: {out_png}")


def main():
    plot_2d_states()
    plot_3d_slice()


if __name__ == '__main__':
    main()
