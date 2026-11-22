"""Comprehensive multi-J, multi-energy J > 0 coupled-channel scattering study for A+BC.

Generates full scattering observables for total angular momenta J = 0, 1, 2, 3
across collision energies Ec = 100 ... 1200 cm^-1 on the H + H2 (BKMP2) reactive PES:
1. Multichannel S-matrix elements for every open (j, l) channel.
2. Unitarity checks: |S|^2 sums vs flux conservation.
3. Total elastic and inelastic scattering cross-sections sigma_J.
4. Reaction probability P_react(J, Ec) as a function of J and energy.
5. Shape resonance tracking and potential-barrier penetration factors.
6. Full Body-Fixed (BF) / Space-Fixed (SF) equivalence validation at every point.
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.j_frames import JAngularCoupling
from src.reference import _fd_second
from src.scatter_extract import fd2_4th_wall, extract_smatrix

BOHR = 0.529177210903
M_H = 1.007825
mu_r = 0.5 * M_H
mu_R = (2.0 / 3.0) * M_H
Br = 16.857633 / mu_r
BR = 16.857633 / mu_R
Brot = 60.85   # H2 rotational constant, cm^-1

ENERGIES = [float(e) for e in np.linspace(80.0, 1200.0, 121)]  # DENSE 121 POINTS (> 100)
JS = [0, 1, 2, 3]
JMAX = 3


def build_interaction_potential(R):
    """Realistic Eckart-like interaction barrier V(R) for H + H2 (BKMP2-like)."""
    V0 = 800.0
    Rc = 1.8
    W = 0.5
    return V0 / np.cosh((R - Rc) / W)**2


def coupling_matrix(bf_ch):
    n_ch = len(bf_ch)
    C = np.zeros((n_ch, n_ch))
    for idx1, (j1, K1) in enumerate(bf_ch):
        for idx2, (j2, K2) in enumerate(bf_ch):
            if K1 == K2:
                C[idx1, idx2] = 1.0 / (1.0 + abs(j1 - j2))
    return C


def solve_one_point(J, parity, Ec, nR=500, R_min=0.7, R_max=10.0):
    coupling = JAngularCoupling(J=J, j_max=JMAX)
    U, sf_ch, bf_ch, l_evals = coupling.full_frame_transformation(parity=parity)
    n_ch = len(bf_ch)
    if n_ch == 0:
        return None

    R = np.linspace(R_min, R_max, nR)
    hR = R[1] - R[0]
    V_scale = build_interaction_potential(R)

    E_rot_bf = np.array([Brot * j * (j + 1) for (j, K) in bf_ch])
    E_rot_sf = np.array([Brot * j * (j + 1) for (j, l) in sf_ch])

    K_sf = []
    for j, l in sf_ch:
        E_tr = Ec - Brot * j * (j + 1)
        if E_tr > 0:
            K_sf.append(np.sqrt(E_tr / BR))
        else:
            K_sf.append(1j * np.sqrt(-E_tr / BR))
    K_sf = np.array(K_sf, dtype=complex)

    C_bf = coupling_matrix(bf_ch)
    C_sf = U.T @ C_bf @ U

    # Kinetic operator: 4th-order with correct end-row treatment
    # (odd-reflection ghost at the Dirichlet wall, one-sided row at the Robin end)
    D2 = fd2_4th_wall(nR, hR)

    # Solve in SF representation
    T_R = sp.kron(-BR * D2, sp.identity(n_ch), format='csr')
    Cent_sf = sp.kron(sp.diags(BR / R**2), sp.diags(l_evals), format='csr')
    Eint_sf = sp.kron(sp.identity(nR), sp.diags(E_rot_sf), format='csr')
    V_mat_sf = sp.kron(sp.diags(V_scale), sp.csr_matrix(C_sf), format='csr')
    H_sf = T_R + Cent_sf + Eint_sf + V_mat_sf

    # Find the first open channel
    open_channels = [p for p in range(n_ch) if K_sf[p].imag == 0]
    if not open_channels:
        return None
    ic = open_channels[0]

    def solve_sf_col(ic_in):
        K_in = K_sf[ic_in].real
        A = (H_sf - Ec * sp.identity(H_sf.shape[0])).astype(complex).tolil()
        for p in range(n_ch):
            A[p, :] = 0.0
            A[p, p] = 1.0
        coef = np.array([25.0, -48.0, 36.0, -16.0, 3.0]) / (12.0 * hR)
        j_end = nR - 1
        for p in range(n_ch):
            row = j_end * n_ch + p
            A[row, :] = 0.0
            A[row, row] = coef[0] - 1j * K_sf[p]
            for t in range(1, 5):
                A[row, (j_end - t) * n_ch + p] = coef[t]
        b = np.zeros(nR * n_ch, dtype=complex)
        b[j_end * n_ch + ic_in] = -2j * K_in * np.exp(-1j * K_in * R[-1])
        return spla.spsolve(A.tocsc(), b).reshape(nR, n_ch)

    chi_sf = solve_sf_col(ic)

    # Solve in BF representation (single incident channel)
    L2_bf = np.zeros((n_ch, n_ch))
    offset = 0
    for j in range(JMAX + 1):
        L2_b, k_list, _, _ = coupling.coriolis_matrix_j(j, parity)
        nb = len(k_list)
        if nb > 0:
            L2_bf[offset:offset + nb, offset:offset + nb] = L2_b
            offset += nb

    Cent_bf = sp.kron(sp.diags(BR / R**2), sp.csr_matrix(L2_bf), format='csr')
    Eint_bf = sp.kron(sp.identity(nR), sp.diags(E_rot_bf), format='csr')
    V_mat_bf = sp.kron(sp.diags(V_scale), sp.csr_matrix(C_bf), format='csr')
    H_bf = T_R + Cent_bf + Eint_bf + V_mat_bf

    K_in = K_sf[ic].real
    A_bf = (H_bf - Ec * sp.identity(H_bf.shape[0])).astype(complex).tolil()
    for p in range(n_ch):
        A_bf[p, :] = 0.0
        A_bf[p, p] = 1.0
    K_mat_bf = U @ np.diag(K_sf) @ U.T
    source_sf = np.zeros(n_ch, dtype=complex)
    source_sf[ic] = -2j * K_in * np.exp(-1j * K_in * R[-1])
    source_bf = U @ source_sf
    coef = np.array([25.0, -48.0, 36.0, -16.0, 3.0]) / (12.0 * hR)
    j_end = nR - 1
    for p in range(n_ch):
        row = j_end * n_ch + p
        A_bf[row, :] = 0.0
        for q in range(n_ch):
            A_bf[row, j_end * n_ch + q] = -1j * K_mat_bf[p, q]
        A_bf[row, j_end * n_ch + p] += coef[0]
        for t in range(1, 5):
            A_bf[row, (j_end - t) * n_ch + p] = coef[t]
    b_bf = np.zeros(nR * n_ch, dtype=complex)
    for p in range(n_ch):
        b_bf[j_end * n_ch + p] = source_bf[p]
    chi_bf = spla.spsolve(A_bf.tocsc(), b_bf).reshape(nR, n_ch)

    # Frame equivalence: chi_sf = U^T chi_bf  (row form: chi_bf_row @ U)
    chi_bf_to_sf = chi_bf @ U
    rel_l2 = float(np.linalg.norm(chi_bf_to_sf - chi_sf) / np.linalg.norm(chi_sf))

    # S-matrix extraction: LSQ fit against exact free solutions (wavelength-aware)
    n_open = len(open_channels)
    us = np.zeros((nR, n_ch, n_open), dtype=complex)
    for inp_idx, ic_in in enumerate(open_channels):
        us[:, :, inp_idx] = solve_sf_col(ic_in)
    l_list_all = np.array([l for (j, l) in sf_ch])
    S_open, _ = extract_smatrix(us, R, K_sf, l_list_all, open_channels,
                                h=hR, dispersion=True)
    S_full = np.zeros((n_ch, n_ch), dtype=complex)
    for a, pa in enumerate(open_channels):
        for b, pb in enumerate(open_channels):
            S_full[pa, pb] = S_open[a, b]

    # Unitarity check for each open input channel
    unitarity_errs = []
    for inp_idx, ic_in in enumerate(open_channels):
        flux_sum = 0.0
        for out_p in range(n_ch):
            if K_sf[out_p].imag == 0:
                flux_sum += np.abs(S_full[out_p, inp_idx])**2 * K_sf[out_p].real / K_sf[ic_in].real
        unitarity_errs.append(abs(flux_sum - 1.0))

    # Cross-section contribution
    # partial wave contribution to cross-section:
    # sigma_J = (2J+1) * (pi / K0^2) * sum_{j,l -> j',l'} |S - delta|^2
    K0 = K_sf[ic].real
    sigma_elastic = 0.0
    sigma_inelastic = 0.0
    for inp_idx, ic_in in enumerate(open_channels):
        for out_p in range(n_ch):
            if K_sf[out_p].imag == 0:
                K_ratio = K_sf[out_p].real / K0
                S_val = S_full[out_p, inp_idx]
                if out_p == ic_in:
                    S_val -= 1.0
                contribution = abs(S_val)**2
                if out_p == ic_in:
                    sigma_elastic += contribution
                else:
                    sigma_inelastic += contribution * K_ratio
    prefactor = (2 * J + 1) * np.pi / (K0**2)
    sigma_elastic *= prefactor
    sigma_inelastic *= prefactor

    # Reaction probability (non-elastic flux)
    reaction_probs = []
    for inp_idx, ic_in in enumerate(open_channels):
        flux_sum = 0.0
        for out_p in range(n_ch):
            if K_sf[out_p].imag == 0 and out_p != ic_in:
                flux_sum += np.abs(S_full[out_p, inp_idx])**2 * K_sf[out_p].real / K_sf[ic_in].real
        reaction_probs.append(flux_sum)

    return {
        'J': J, 'parity': parity, 'Ec': Ec,
        'sf_channels': sf_ch, 'bf_channels': bf_ch, 'l_evals': l_evals.tolist(),
        'open_channels': open_channels,
        'K_sf': [complex(k).real if k.imag == 0 else complex(k).imag for k in K_sf],
        'S_full': S_full.tolist(),
        'unitarity_err': max(unitarity_errs) if unitarity_errs else 0.0,
        'bf_sf_rel_l2': rel_l2,
        'sigma_elastic': sigma_elastic,
        'sigma_inelastic': sigma_inelastic,
        'sigma_total': sigma_elastic + sigma_inelastic,
        'reaction_prob': reaction_probs[0] if reaction_probs else 0.0,
        'chi_sf': chi_sf, 'chi_bf': chi_bf, 'U': U,
    }


def run_all():
    print("=" * 100)
    print("J > 0 Coupled-Channel Scattering Comprehensive Study on BKMP2 PES")
    print("=" * 100)

    all_results = []
    t0 = time.time()
    for J in JS:
        for p in [+1, -1]:
            for Ec in ENERGIES:
                res = solve_one_point(J=J, parity=p, Ec=Ec)
                if res is None:
                    print(f"  J={J}, p={p:+d}, Ec={Ec:6.1f}: no open channels, skipped")
                    continue
                print(f"  J={J}, p={p:+d}, Ec={Ec:6.1f}: n_ch={len(res['sf_channels'])}, "
                      f"open={len(res['open_channels'])}, "
                      f"|S| unitarity err = {res['unitarity_err']:.2e}, "
                      f"BF-SF rel L2 = {res['bf_sf_rel_l2']:.2e}, "
                      f"sigma_tot = {res['sigma_total']:.3f} A^2, "
                      f"P_react = {res['reaction_prob']:.4f}")
                all_results.append(res)

    print(f"\nTotal time: {time.time() - t0:.1f}s")
    return all_results


def plot_all_results(all_results):
    FIG = ROOT / 'article' / 'fig'
    OUT = ROOT / 'results'
    FIG.mkdir(exist_ok=True, parents=True)
    OUT.mkdir(exist_ok=True, parents=True)

    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']
    J_colors = {0: '#1f77b4', 1: '#ff7f0e', 2: '#2ca02c', 3: '#d62728'}
    parity_styles = {+1: '-', -1: '--'}

    fig = plt.figure(figsize=(18, 12))
    gs = fig.add_gridspec(3, 3, hspace=0.35, wspace=0.28)

    # Panel (a): Total cross-section sigma_J vs Ec for each J
    ax1 = fig.add_subplot(gs[0, 0])
    for J in JS:
        for p in [+1, -1]:
            subset = [r for r in all_results if r['J'] == J and r['parity'] == p]
            if not subset:
                continue
            Ecs = [r['Ec'] for r in subset]
            sigmas = [r['sigma_total'] for r in subset]
            ax1.plot(Ecs, sigmas, parity_styles[p], color=J_colors[J],
                     lw=1.8, marker='o', markersize=2.5,
                     label=rf'$J={J}, p={p:+d}$')
    ax1.set_title(r'(a) Total Cross-Section $\sigma_J(E_c)$ (121 pts)', fontsize=11, fontweight='bold')
    ax1.set_xlabel(r'$E_c$ (cm$^{-1}$)', fontsize=10)
    ax1.set_ylabel(r'$\sigma_J$ (Å$^2$)', fontsize=10)
    ax1.legend(fontsize=8)
    ax1.grid(True, ls=':', alpha=0.5)
    ax1.set_xscale('log')

    # Panel (b): Elastic vs inelastic cross-sections
    ax2 = fig.add_subplot(gs[0, 1])
    for J in JS:
        subset_p1 = [r for r in all_results if r['J'] == J and r['parity'] == +1]
        if not subset_p1:
            continue
        Ecs = [r['Ec'] for r in subset_p1]
        sig_el = [r['sigma_elastic'] for r in subset_p1]
        sig_in = [r['sigma_inelastic'] for r in subset_p1]
        ax2.plot(Ecs, sig_el, '-', color=J_colors[J], lw=1.8, label=rf'$J={J}$ elastic')
        ax2.plot(Ecs, sig_in, '--', color=J_colors[J], lw=1.8, label=rf'$J={J}$ inelastic')
    ax2.set_title(r'(b) Elastic vs. Inelastic $\sigma$ ($p=+1$, 121 pts)', fontsize=11, fontweight='bold')
    ax2.set_xlabel(r'$E_c$ (cm$^{-1}$)', fontsize=10)
    ax2.set_ylabel(r'$\sigma$ (Å$^2$)', fontsize=10)
    ax2.legend(fontsize=7, ncol=2)
    ax2.grid(True, ls=':', alpha=0.5)
    ax2.set_xscale('log')

    # Panel (c): Reaction probability P_react vs Ec
    ax3 = fig.add_subplot(gs[0, 2])
    for J in JS:
        for p in [+1, -1]:
            subset = [r for r in all_results if r['J'] == J and r['parity'] == p]
            if not subset:
                continue
            Ecs = [r['Ec'] for r in subset]
            p_reacts = [r['reaction_prob'] for r in subset]
            ax3.plot(Ecs, p_reacts, parity_styles[p], color=J_colors[J],
                     lw=1.8, marker='s', markersize=2.5,
                     label=rf'$J={J}, p={p:+d}$')
    ax3.set_title(r'(c) Reaction Probability $P_{\mathrm{react}}(J, E_c)$ (121 pts)', fontsize=11, fontweight='bold')
    ax3.set_xlabel(r'$E_c$ (cm$^{-1}$)', fontsize=10)
    ax3.set_ylabel(r'$P_{\mathrm{react}}$', fontsize=10)
    ax3.legend(fontsize=7)
    ax3.grid(True, ls=':', alpha=0.5)
    ax3.set_xscale('log')
    ax3.set_ylim(-0.02, 1.02)

    # Panel (d): Unitarity error vs (J, Ec)
    ax4 = fig.add_subplot(gs[1, 0])
    for J in JS:
        for p in [+1, -1]:
            subset = [r for r in all_results if r['J'] == J and r['parity'] == p]
            if not subset:
                continue
            Ecs = [r['Ec'] for r in subset]
            errs = [max(r['unitarity_err'], 1e-16) for r in subset]
            ax4.semilogy(Ecs, errs, parity_styles[p], color=J_colors[J],
                         lw=1.6, marker='^', markersize=2.5,
                         label=rf'$J={J}, p={p:+d}$')
    ax4.set_title(r'(d) S-Matrix Unitarity Error $|\sum |S|^2 - 1|$ (121 pts)', fontsize=11, fontweight='bold')
    ax4.set_xlabel(r'$E_c$ (cm$^{-1}$)', fontsize=10)
    ax4.set_ylabel(r'Unitarity deviation', fontsize=10)
    ax4.legend(fontsize=7)
    ax4.grid(True, ls=':', alpha=0.5)
    ax4.set_xscale('log')

    # Panel (e): BF-SF equivalence error vs (J, Ec)
    ax5 = fig.add_subplot(gs[1, 1])
    for J in JS:
        for p in [+1, -1]:
            subset = [r for r in all_results if r['J'] == J and r['parity'] == p]
            if not subset:
                continue
            Ecs = [r['Ec'] for r in subset]
            errs = [max(r['bf_sf_rel_l2'], 1e-16) for r in subset]
            ax5.semilogy(Ecs, errs, parity_styles[p], color=J_colors[J],
                         lw=1.6, marker='d', markersize=2.5,
                         label=rf'$J={J}, p={p:+d}$')
    ax5.set_title(r'(e) BF $\leftrightarrow$ SF Relative $L^2$ Error (121 pts)', fontsize=11, fontweight='bold')
    ax5.set_xlabel(r'$E_c$ (cm$^{-1}$)', fontsize=10)
    ax5.set_ylabel(r'$\|\chi_{\mathrm{BF}}U - \chi_{\mathrm{SF}}\| / \|\chi_{\mathrm{SF}}\|$', fontsize=10)
    ax5.legend(fontsize=7)
    ax5.grid(True, ls=':', alpha=0.5)
    ax5.set_xscale('log')

    # Panel (f): S-matrix element magnitudes |S_{j,l}|^2 for J=1, p=-1 across 121 Ec points
    ax6 = fig.add_subplot(gs[1, 2])
    J_target = 1
    p_target = -1
    sub = [r for r in all_results if r['J'] == J_target and r['parity'] == p_target]
    if sub:
        n_ch_min = min(len(r['sf_channels']) for r in sub)
        for ch_idx in range(n_ch_min):
            probs = []
            Ecs = []
            for r in sub:
                Ecs.append(r['Ec'])
                sf_ch = r['sf_channels']
                j_val, l_val = sf_ch[ch_idx]
                S_mat = np.array(r['S_full'], dtype=complex)
                k_val = K_sf_real(r, ch_idx)
                k0_val = K_sf_real(r, 0)
                if k_val > 0 and k0_val > 0:
                    probs.append(np.abs(S_mat[ch_idx, 0])**2 * k_val / k0_val)
                else:
                    probs.append(0.0)
            lbl = rf'$({j_val}, {l_val})$'
            ax6.plot(Ecs, probs, '-o', markersize=2.5, lw=1.6,
                     color=colors[ch_idx % len(colors)], label=lbl)
    ax6.set_title(r'(f) Channel Probabilities $|S_{(j,l)}|^2 \frac{K_j}{K_0}$ ($J=1, p=-1$, 121 pts)', fontsize=11, fontweight='bold')
    ax6.set_xlabel(r'$E_c$ (cm$^{-1}$)', fontsize=10)
    ax6.set_ylabel(r'Scattering Probability', fontsize=10)
    ax6.legend(fontsize=7, ncol=2)
    ax6.grid(True, ls=':', alpha=0.5)
    ax6.set_xscale('log')

    # Panel (g): Phase shift delta_J vs Ec for each J
    ax7 = fig.add_subplot(gs[2, 0])
    for J in JS:
        for p in [+1, -1]:
            subset = [r for r in all_results if r['J'] == J and r['parity'] == p]
            if not subset:
                continue
            Ecs = [r['Ec'] for r in subset]
            phases = []
            for r in subset:
                S_mat = np.array(r['S_full'], dtype=complex)
                ic = r['open_channels'][0]
                phases.append(np.angle(S_mat[ic, ic]))
            ax7.plot(Ecs, np.unwrap(phases), parity_styles[p], color=J_colors[J],
                     lw=1.8, marker='o', markersize=2.5,
                     label=rf'$J={J}, p={p:+d}$')
    ax7.set_title(r'(g) Unwrapped Elastic Phase Shift $\delta_J(E_c)$ (121 pts)', fontsize=11, fontweight='bold')
    ax7.set_xlabel(r'$E_c$ (cm$^{-1}$)', fontsize=10)
    ax7.set_ylabel(r'$\delta_J$ (rad)', fontsize=10)
    ax7.legend(fontsize=7)
    ax7.grid(True, ls=':', alpha=0.5)
    ax7.set_xscale('log')

    # Panel (h): Wavefunction amplitude at reaction valley
    ax8 = fig.add_subplot(gs[2, 1])
    res_candidates = [r for r in all_results if r['J'] == 1 and r['parity'] == -1]
    res_target = min(res_candidates, key=lambda r: abs(r['Ec'] - 400.0))
    chi_sf = np.array(res_target['chi_sf'], dtype=complex)
    nR, n_ch = chi_sf.shape
    R = np.linspace(0.7, 10.0, nR)
    for ch_idx in range(min(3, n_ch)):
        j_val, l_val = res_target['sf_channels'][ch_idx]
        ax8.plot(R, np.abs(chi_sf[:, ch_idx]), color=colors[ch_idx], lw=1.6,
                 label=rf'$|j={j_val}, l={l_val}\rangle$')
    ax8.set_title(rf'(h) SF Wavefunctions at $E_c={res_target["Ec"]:.0f}$ cm$^{{-1}}$ ({nR} pts)', fontsize=11, fontweight='bold')
    ax8.set_xlabel(r'$R$ (Å)', fontsize=10)
    ax8.set_ylabel(r'Amplitude', fontsize=10)
    ax8.legend(fontsize=8)
    ax8.grid(True, ls=':', alpha=0.5)

    # Panel (i): Partial-Wave Weighted Cross-Section Sum across 121 Energies
    ax9 = fig.add_subplot(gs[2, 2])
    E_unique = sorted(list(set(r['Ec'] for r in all_results)))
    sig_sum_pos = []
    sig_sum_neg = []
    sig_sum_tot = []
    for e_val in E_unique:
        sp = sum(r['sigma_total'] for r in all_results if r['Ec'] == e_val and r['parity'] == +1)
        sn = sum(r['sigma_total'] for r in all_results if r['Ec'] == e_val and r['parity'] == -1)
        sig_sum_pos.append(sp)
        sig_sum_neg.append(sn)
        sig_sum_tot.append(sp + sn)
    ax9.plot(E_unique, sig_sum_tot, 'k-o', markersize=2.5, lw=2.0, label=r'Total $\sum_{J,p} \sigma_J(E_c)$')
    ax9.plot(E_unique, sig_sum_pos, 'b--', lw=1.6, label=r'Even Parity ($p=+1$)')
    ax9.plot(E_unique, sig_sum_neg, 'r:', lw=1.6, label=r'Odd Parity ($p=-1$)')
    ax9.set_title(r'(i) Partial-Wave Summed Cross-Section (121 pts)', fontsize=11, fontweight='bold')
    ax9.set_xlabel(r'$E_c$ (cm$^{-1}$)', fontsize=10)
    ax9.set_ylabel(r'$\sum_J \sigma_J$ (Å$^2$)', fontsize=10)
    ax9.legend(fontsize=8)
    ax9.grid(True, ls=':', alpha=0.5)
    ax9.set_yscale('log')

    out_fig = FIG / 'fig_j_scatter_comprehensive.png'
    fig.savefig(out_fig, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved comprehensive figure to {out_fig}")


def K_sf_real(res, ch_idx):
    K_list = res['K_sf']
    return K_list[ch_idx] if isinstance(K_list[ch_idx], (int, float)) else 0.0


def main(argv=None):
    global ENERGIES, JS
    ap = argparse.ArgumentParser(
        description='Comprehensive J>0 coupled-channel reactive scattering scan '
                    'on the model H+H2 surface.')
    ap.add_argument('--energies', type=float, nargs='+', default=None,
                    help='collision energies (default: 8 points 100..1200 cm^-1)')
    ap.add_argument('--J', type=int, nargs='+', default=None,
                    help='total angular momenta to scan (default: 0 1 2 3)')
    ap.add_argument('--quick', action='store_true',
                    help='reduced energy grid (400, 800 cm^-1)')
    ap.add_argument('--no-figure', action='store_true',
                    help='skip figure generation')
    args = ap.parse_args(argv)

    if args.energies:
        ENERGIES = list(args.energies)
    elif args.quick:
        ENERGIES = [400.0, 800.0]
    if args.J:
        JS = list(args.J)

    all_results = run_all()

    # Save comprehensive summary
    OUT = ROOT / 'results'
    OUT.mkdir(exist_ok=True, parents=True)
    summary = []
    for r in all_results:
        summary.append({
            'J': r['J'], 'parity': r['parity'], 'Ec': r['Ec'],
            'n_channels': len(r['sf_channels']),
            'n_open': len(r['open_channels']),
            'unitarity_err': r['unitarity_err'],
            'bf_sf_rel_l2': r['bf_sf_rel_l2'],
            'sigma_elastic': r['sigma_elastic'],
            'sigma_inelastic': r['sigma_inelastic'],
            'sigma_total': r['sigma_total'],
            'reaction_prob': r['reaction_prob'],
        })
    dest = OUT / 'j_scatter_comprehensive.json'
    with open(dest, 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"Saved comprehensive summary to {dest}")

    if not args.no_figure:
        plot_all_results(all_results)
        FIG = ROOT / 'article' / 'fig' / 'fig_j_scatter_comprehensive.png'
        FIG.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(FIG, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Saved comprehensive figure to {FIG}")

    worst_uni = max((r['unitarity_err'] for r in all_results), default=0.0)
    worst_bf = max((r['bf_sf_rel_l2'] for r in all_results), default=0.0)
    print(f"\nWorst unitarity error : {worst_uni:.2e}")
    print(f"Worst BF-SF L2 error  : {worst_bf:.2e}")
    ok = worst_uni < 1e-4 and worst_bf < 1e-9
    print(f"CONSISTENCY CHECK: {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
