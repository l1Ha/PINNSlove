"""Re-plot the 2D BKMP2 H + H2 scattering comparison figures from saved NPZ data.

Usage:  python3 scripts/plot_h2_bkmp2.py [--energies 400 800]
"""
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
FIG = ROOT / 'article' / 'fig'

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument('--energies', type=int, nargs='+', default=[400, 800],
                help='collision energies whose NPZ files should be plotted')
args = ap.parse_args()

for Ec in args.energies:
    d = np.load(ROOT / f'results/h2_scatter_{Ec}_bkmp2.npz', allow_pickle=True)
    R = d['R']
    gamma = d['gamma']
    chi_ref = d['chi_ref']
    chi_cc = d['chi_cc']
    m = d['metrics'].item()
    K = m['K_inv_A']
    l2 = m['l2_relative_error']
    
    fig, axes = plt.subplots(2, 3, figsize=(13, 7.5), sharex=True, sharey=True)
    vmax = max(np.abs(chi_ref).max(), np.abs(chi_cc).max())
    
    gam_deg = gamma * 180.0 / np.pi
    
    # Row 0: Classical FD
    im0 = axes[0, 0].pcolormesh(R, gam_deg, np.abs(chi_ref).T, cmap='viridis', shading='auto', vmin=0, vmax=vmax)
    axes[0, 0].set_title(r'Reference FD: $|\chi(R, \gamma)|$', fontsize=11)
    axes[0, 0].set_ylabel(r'$\gamma$ (deg)', fontsize=11)
    fig.colorbar(im0, ax=axes[0, 0], fraction=0.046, pad=0.04)
    
    im1 = axes[0, 1].pcolormesh(R, gam_deg, np.real(chi_ref).T, cmap='coolwarm', shading='auto', vmin=-vmax, vmax=vmax)
    axes[0, 1].set_title(r'Reference FD: $\mathrm{Re}[\chi]$', fontsize=11)
    fig.colorbar(im1, ax=axes[0, 1], fraction=0.046, pad=0.04)
    
    im2 = axes[0, 2].pcolormesh(R, gam_deg, np.imag(chi_ref).T, cmap='coolwarm', shading='auto', vmin=-vmax, vmax=vmax)
    axes[0, 2].set_title(r'Reference FD: $\mathrm{Im}[\chi]$', fontsize=11)
    fig.colorbar(im2, ax=axes[0, 2], fraction=0.046, pad=0.04)
    
    # Row 1: Close-coupling
    im3 = axes[1, 0].pcolormesh(R, gam_deg, np.abs(chi_cc).T, cmap='viridis', shading='auto', vmin=0, vmax=vmax)
    axes[1, 0].set_title(rf'Close-coupling: $|\chi(R, \gamma)|$ (Rel $L^2$: {l2*100:.1f}%)', fontsize=11)
    axes[1, 0].set_xlabel(r'$R$ (Å)', fontsize=11)
    axes[1, 0].set_ylabel(r'$\gamma$ (deg)', fontsize=11)
    fig.colorbar(im3, ax=axes[1, 0], fraction=0.046, pad=0.04)
    
    im4 = axes[1, 1].pcolormesh(R, gam_deg, np.real(chi_cc).T, cmap='coolwarm', shading='auto', vmin=-vmax, vmax=vmax)
    axes[1, 1].set_title(r'Close-coupling: $\mathrm{Re}[\chi]$', fontsize=11)
    axes[1, 1].set_xlabel(r'$R$ (Å)', fontsize=11)
    fig.colorbar(im4, ax=axes[1, 1], fraction=0.046, pad=0.04)
    
    diff = np.abs(chi_cc - chi_ref)
    im5 = axes[1, 2].pcolormesh(R, gam_deg, diff.T, cmap='magma', shading='auto')
    axes[1, 2].set_title(r'Absolute Error $|\chi_{\mathrm{Close-coupling}} - \chi_{\mathrm{FD}}|$', fontsize=11)
    axes[1, 2].set_xlabel(r'$R$ (Å)', fontsize=11)
    fig.colorbar(im5, ax=axes[1, 2], fraction=0.046, pad=0.04)
    
    fig.suptitle(f'H + H2(v=0) Reactive Scattering on BKMP2 PES ($E_c = {Ec}$ ' + r'$\mathrm{cm}^{-1}$' + f', $K = {K:.3f}$ ' + r'$\mathrm{\AA}^{-1}$' + ')', fontsize=13, y=0.98)
    plt.tight_layout()
    plt.savefig(FIG / f'fig_scatter_h2_bkmp2_{Ec}.png', dpi=300)
    plt.close()
    print(f'Successfully replotted fig_scatter_h2_bkmp2_{Ec}.png')
