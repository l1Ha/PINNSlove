"""Generate all result figures from the saved .npz/.csv/.json files.

Usage:  python3 scripts/make_figures.py
Missing inputs are reported and skipped; existing figures are regenerated.
"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))          # allow "import src.*" from any cwd
RES = ROOT / 'results'
plt.rcParams.update({'font.size': 10, 'figure.dpi': 130})


def fig_1d_bound():
    d = np.load(RES / '1d_bound' / 'wavefunctions.npz')
    from src.reference import gauss_legendre_1d as _gl
    x, _ = _gl(-1.6, 4.6, d['psi_pool'].shape[1])
    x = np.sort(x)
    order = np.argsort(_gl(-1.6, 4.6, d['psi_pool'].shape[1])[0])
    V = d['V']
    psi = d['psi_pool'][:, order]
    summary = json.load(open(RES / '1d_bound' / 'summary.json'))
    E = np.array([r['pinn'] for r in summary['rows']])
    E_ex = np.array([r['exact'] for r in summary['rows']])
    fig, ax = plt.subplots(figsize=(7, 5))
    xf = d['x_plot']
    ax.plot(xf, d['V'], 'k-', lw=1.2, label='Morse $V(x)$')
    scale = 2600.0
    for n in range(len(E)):
        off = E_ex[n] + 1300.0
        ax.plot(x, off + scale * psi[n] / np.max(np.abs(psi[n])), 'b-',
                lw=0.9, label='PINN' if n == 0 else None)
        ax.axhline(off, color='gray', lw=0.4, ls=':')
    ax.set_ylim(-38293 - 3000, 8000)
    ax.set_xlabel('x  [Å]')
    ax.set_ylabel('energy / wavefunction  [cm$^{-1}$]')
    ax.set_title('1D Morse bound states: PINN vs exact levels\n'
                 f'max |$\\Delta$E| = {np.max(np.abs(E - E_ex)):.1f} cm$^{{-1}}$'
                 f' (states 0-{len(E)-1})')
    ax.legend(loc='upper right')
    fig.tight_layout()
    fig.savefig(RES / '1d_bound' / 'fig_1d_bound.png')
    plt.close(fig)


def fig_1d_scatter():
    d = np.load(RES / '1d_scatter' / 'wavefunctions.npz')
    x, V = d['x'], d['V']
    Es = d['Es']
    fig, axes = plt.subplots(len(Es), 1, figsize=(7, 2.4 * len(Es)),
                             sharex=True)
    axes = np.atleast_1d(axes)
    for ax, E, ps_p, ps_r in zip(axes, Es, d['psi_pinn'], d['psi_ref']):
        ax.plot(x, ps_r.real, 'k-', lw=1.0, label='Re reference')
        ax.plot(x, ps_r.imag, 'k--', lw=1.0, label='Im reference')
        ax.plot(x, ps_p.real, 'r-', lw=0.9, label='Re PINN')
        ax.plot(x, ps_p.imag, 'r--', lw=0.9, label='Im PINN')
        ax.set_ylabel(f'$\\psi$  (E={float(E):.0f})')
    axes[0].legend(ncol=2, fontsize=8)
    axes[-1].set_xlabel('x  [Å]')
    fig.tight_layout()
    fig.savefig(RES / '1d_scatter' / 'fig_1d_scatter.png')
    plt.close(fig)


def fig_2d_water():
    d = np.load(RES / '2d_water' / 'wavefunctions.npz')
    R, gam, psi = d['R'], d['gamma'], d['psi']
    summary = json.load(open(RES / '2d_water' / 'summary.json'))
    E = np.array([r['pinn'] for r in summary['rows']])
    Er = np.array([r['dvr'] for r in summary['rows']])
    n = min(4, len(E))
    fig, axes = plt.subplots(1, n, figsize=(3.0 * n, 3.0))
    for i, ax in enumerate(axes[:n]):
        P = psi[i]
        lv = np.max(np.abs(P))
        ax.contourf(np.rad2deg(gam), R, P, levels=20, cmap='RdBu_r',
                    vmin=-lv, vmax=lv)
        ax.contour(np.rad2deg(gam), R, P, levels=[0.0], colors='k',
                   linewidths=0.5)
        ax.set_title(f'state {i}: +{E[i]-Er[0]:.0f} cm$^{{-1}}$\n'
                     f'(DVR +{Er[i]-Er[0]:.0f})', fontsize=9)
        ax.set_xlabel('$\\gamma$  [deg]')
    axes[0].set_ylabel('R  [Å]')
    fig.suptitle('2D water (frozen OH): PINN wavefunctions', y=1.02)
    fig.tight_layout()
    fig.savefig(RES / '2d_water' / 'fig_2d_states.png', bbox_inches='tight')
    plt.close(fig)


def fig_3d_water():
    d = np.load(RES / '3d_water' / 'slices.npz')
    r, R, chi0, V = d['r'], d['R'], d['chi0'], d['V_slice']
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.6))
    im0 = axes[0].contourf(r, R, chi0.T, levels=24, cmap='viridis')
    axes[0].set_title('PINN ground state $\\chi(r,R)$ at $\\gamma=\\gamma_e$')
    fig.colorbar(im0, ax=axes[0])
    im1 = axes[1].contourf(r, R, V.T, levels=np.linspace(-81000, -40000, 25),
                           cmap='inferno')
    axes[1].set_title('PES $V(r,R)$ at $\\gamma=\\gamma_e$')
    fig.colorbar(im1, ax=axes[1])
    for ax in axes:
        ax.set_xlabel('r  [Å]')
        ax.set_ylabel('R  [Å]')
    fig.tight_layout()
    fig.savefig(RES / '3d_water' / 'fig_3d_slice.png')
    plt.close(fig)


if __name__ == '__main__':
    for f in [fig_1d_bound, fig_1d_scatter, fig_2d_water, fig_3d_water]:
        try:
            f()
            print(f.__name__, 'ok')
        except Exception as e:
            print(f.__name__, 'SKIPPED:', repr(e))
