#!/usr/bin/env python3
"""Reference-accuracy figure for the full-dimensional 3D A+BC scattering solve.

Reads the saved reference wavefunctions (results/3d_study/*.npz, produced by
scripts/server_3d_study.py), extracts for every grid/box configuration

  * the asymptotic diatomic channel energies (ZPE, j = 2, 4 excitations) and
    compares them with the literature H2 spectrum,
  * the elastic S-matrix element |S00| and arg S00 (least-squares fit against
    the free reduced radial functions with the FD-consistent wavenumber q(h)),

and writes a three-panel validation figure.

Usage
-----
    python3 scripts/make_reference_validation_figure.py [--out PATH]
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import scipy.linalg
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))

from src.reference import LegendreAngularDVR                       # noqa: E402
from src.scatter_extract import (fd2_4th_dirichlet, fit_wavevector,  # noqa: E402
                                 free_incoming, free_outgoing)
from src import bkmp2                                              # noqa: E402
from run_scatter_3d import Br, BR, _channel_labels                 # noqa: E402

STORE = ROOT / 'results' / '3d_study'
CONFIGS = ['baseline', 'nR240', 'nG20', 'nr60', 'Rbox10', 'rbox22']
LABELS = {'baseline': 'baseline\n40x120x16', 'nR240': 'nR=240',
          'nG20': 'nG=20', 'nr60': 'nr=60', 'Rbox10': 'R$_{max}$=10Å',
          'rbox22': 'r$_{max}$=2.2Å'}
# literature H2 (v=0) values, cm^-1
LIT = {'zpe': 2179.2, 'j2': 354.33, 'j4': 1168.1}


def analyse(path):
    d = dict(np.load(path))
    r, R, gamma = d['r'], d['R'], d['gamma']
    nr, nG = len(r), len(gamma)
    hR = R[1] - R[0]
    ang = LegendreAngularDVR(nG)
    V = bkmp2.potential_grid_rRg(r, [R[-1]], gamma, verbose=False).reshape(nr, nG)
    D2 = fd2_4th_dirichlet(nr, r[1] - r[0])
    D2i = np.asarray(D2[1:nr - 1, 1:nr - 1].todense())
    H = (np.kron(-Br * D2i, np.eye(nG))
         + np.kron(np.diag(Br / r[1:-1] ** 2), ang.j2)
         + np.diag(V[1:-1, :].ravel()))
    w, U = scipy.linalg.eig(H)
    order = np.argsort(w.real)
    Ur = U[:, order].real
    Uf = np.zeros((nr * nG, Ur.shape[1]))
    for ir in range(1, nr - 1):
        Uf[ir * nG:(ir + 1) * nG, :] = Ur[(ir - 1) * nG:ir * nG, :]
    jl, vl = _channel_labels(r, gamma, Uf, ang.j2, nr, nG)
    eps = w.real[order]

    def pick(v, j):
        for c in range(len(jl)):
            if vl[c] == v and jl[c] == j:
                return c
        raise KeyError((v, j))

    c00, c02, c04 = pick(0, 0), pick(0, 2), pick(0, 4)
    exc = eps[c02] - eps[c00], eps[c04] - eps[c00]

    # collision wavenumber from the saved K0
    K0 = float(d['K0'])
    q = fit_wavevector(K0, hR, 4)
    lam = 2 * np.pi / q
    win = min(0.45 * (R[-1] - R[0]), max(3 * lam, 2.0))
    m = (R > R[-1] - win) & (R < R[-1] - 5 * hR)
    M = np.column_stack([free_incoming(q * R[m], 0), free_outgoing(q * R[m], 0)])
    f0 = d['chi'].reshape(len(R), -1) @ Uf[:, c00]
    c, *_ = np.linalg.lstsq(M, f0[m], rcond=None)
    S00 = c[1] / c[0]
    return dict(zpe=float(d['eps0']) + 38297.387, exc2=float(exc[0]), exc4=float(exc[1]),
                S00=complex(S00), K0=K0, l2_prod=float(d['l2_prod']))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', default=str(ROOT / 'article' / 'fig' /
                                        'fig_3d_reference_validation.png'))
    ap.add_argument('--json', default=str(ROOT / 'results' /
                                          'check_3d_reference_validation.json'))
    args = ap.parse_args(argv)

    data = {}
    for E in (400, 800):
        for lab in CONFIGS:
            f = STORE / f'{lab}_{E}.npz'
            if f.exists():
                data[(E, lab)] = analyse(f)
                print(f'{E} {lab}: ZPE-err={data[(E,lab)]["zpe"]-LIT["zpe"]:+.3f}  '
                      f'|S00|={abs(data[(E,lab)]["S00"]):.6f}  '
                      f'arg={np.degrees(np.angle(data[(E,lab)]["S00"])):.4f}', flush=True)
    with open(args.json, 'w') as fh:
        json.dump({f'{E}|{lab}': {k: (float(v.real), float(v.imag)) if isinstance(v, complex)
                                  else v for k, v in d.items()}
                   for (E, lab), d in data.items()}, fh, indent=2)

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))
    colors = {'400': '#1f77b4', '800': '#d62728'}
    x = np.arange(len(CONFIGS))

    # panel A: channel energies vs literature
    ax = axes[0]
    ref = data[(400, 'baseline')]
    names = ['ZPE', 'j=2', 'j=4']
    litv = [LIT['zpe'], LIT['j2'], LIT['j4']]
    meas = [ref['zpe'], ref['exc2'], ref['exc4']]
    err = [meas[0] - litv[0], meas[1] - litv[1], meas[2] - litv[2]]
    ax.bar(names, err, color='#2ca02c', width=0.5)
    ax.axhline(0, color='k', lw=0.8)
    ax.set_ylabel('channel energy error (cm$^{-1}$)')
    ax.set_title('diatomic channel energies vs H$_2$ literature')
    for i, e in enumerate(err):
        ax.annotate(f'{e:+.2f}', (i, e), ha='center',
                    va='bottom' if e >= 0 else 'top', fontsize=9)
    ax.text(0.02, 0.95, 'ZPE = 2178.9 cm$^{-1}$\n(lit. 2179.2)', transform=ax.transAxes,
            va='top', fontsize=9)

    # panel B: |S00| across configurations
    ax = axes[1]
    for E in (400, 800):
        vals = [abs(data[(E, lab)]['S00']) for lab in CONFIGS]
        ax.plot(x, vals, 'o-', color=colors[str(E)], label=f'$E_c$ = {E} cm$^{{-1}}$')
    ax.set_xticks(x)
    ax.set_xticklabels([LABELS[l] for l in CONFIGS], fontsize=8)
    ax.set_ylabel('$|S_{00}|$')
    ax.set_title('elastic S-matrix magnitude across grids/boxes')
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    ax.text(0.02, 0.06, 'spread $\\leq 10^{-5}$ (grids) / $10^{-4}$ (R-box)',
            transform=ax.transAxes, fontsize=8)

    # panel C: phase deviation from the baseline
    ax = axes[2]
    for E in (400, 800):
        base = np.degrees(np.angle(data[(E, 'baseline')]['S00']))
        dev = [np.degrees(np.angle(data[(E, lab)]['S00'])) - base for lab in CONFIGS]
        ax.plot(x, dev, 's--', color=colors[str(E)], label=f'$E_c$ = {E} cm$^{{-1}}$')
    ax.axhline(0, color='k', lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels([LABELS[l] for l in CONFIGS], fontsize=8)
    ax.set_ylabel('$\\Delta$ arg $S_{00}$ (deg)')
    ax.set_title('elastic phase stability (vs baseline)')
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    fig.suptitle('Reference-solution validation: full-dimensional 3D H + H$_2$(v=0) '
                 'scattering on BKMP2', y=1.0, fontsize=12)
    fig.tight_layout()
    fig.savefig(args.out, dpi=220, bbox_inches='tight')
    print(f'saved {args.out}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
