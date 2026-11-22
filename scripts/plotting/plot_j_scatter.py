"""Plotting script for J > 0 Reactive Scattering in Body-Fixed (BF) and Space-Fixed (SF) Frames.
Generates comprehensive comparison figures showing channel wavefunctions, frame equivalence,
Coriolis coupling structure, and S-matrix channel probabilities.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.j_frames import JAngularCoupling
from src.reference import _fd_second


def generate_j1_scatter_figure():
    J = 1
    parity = -1
    Ec = 400.0
    nR = 300
    R_min = 0.7
    R_max = 8.5

    coupling = JAngularCoupling(J=J, j_max=2)
    U, sf_ch, bf_ch, l_evals = coupling.full_frame_transformation(parity=parity)
    n_ch = len(bf_ch)

    R = np.linspace(R_min, R_max, nR)
    hR = R[1] - R[0]
    BR = coupling.BR

    # Coriolis L^2 matrix in BF
    L2_bf = np.zeros((n_ch, n_ch))
    offset = 0
    for j in range(3):
        L2_b, k_list, _, _ = coupling.coriolis_matrix_j(j, parity)
        nb = len(k_list)
        if nb > 0:
            L2_bf[offset:offset + nb, offset:offset + nb] = L2_b
            offset += nb

    L2_sf = np.diag(l_evals)

    Brot = 60.85
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

    V0 = 800.0
    V_scale = V0 / np.cosh((R - 1.8) / 0.5)**2

    C_bf = np.zeros((n_ch, n_ch))
    for idx1, (j1, K1) in enumerate(bf_ch):
        for idx2, (j2, K2) in enumerate(bf_ch):
            if K1 == K2:
                C_bf[idx1, idx2] = 1.0 / (1.0 + abs(j1 - j2))
    C_sf = U.T @ C_bf @ U

    # 1. Solve SF
    T_R_sf = sp.kron(-BR * _fd_second(nR, hR, 4), sp.identity(n_ch), format='csr')
    Cent_sf = sp.kron(sp.diags(BR / R**2), sp.diags(l_evals), format='csr')
    Eint_sf = sp.kron(sp.identity(nR), sp.diags(E_rot_sf), format='csr')
    V_mat_sf = sp.kron(sp.diags(V_scale), sp.csr_matrix(C_sf), format='csr')

    H_sf = T_R_sf + Cent_sf + Eint_sf + V_mat_sf
    A_sf = (H_sf - Ec * sp.identity(H_sf.shape[0])).astype(complex).tolil()
    b_sf = np.zeros(nR * n_ch, dtype=complex)

    for p in range(n_ch):
        A_sf[p, :] = 0.0; A_sf[p, p] = 1.0; b_sf[p] = 0.0

    coef = np.array([3.0, -4.0, 1.0]) / (2.0 * hR)
    K0 = K_sf[0].real
    source_sf = np.zeros(n_ch, dtype=complex)
    source_sf[0] = -2j * K0 * np.exp(-1j * K0 * R[-1])

    j_end = nR - 1
    for p in range(n_ch):
        row = j_end * n_ch + p
        A_sf[row, :] = 0.0
        A_sf[row, j_end * n_ch + p] = coef[0] - 1j * K_sf[p]
        A_sf[row, (j_end - 1) * n_ch + p] = coef[1]
        A_sf[row, (j_end - 2) * n_ch + p] = coef[2]
        b_sf[row] = source_sf[p]

    chi_sf = spla.spsolve(A_sf.tocsc(), b_sf).reshape(nR, n_ch)

    # 2. Solve BF
    T_R_bf = sp.kron(-BR * _fd_second(nR, hR, 4), sp.identity(n_ch), format='csr')
    Cent_bf = sp.kron(sp.diags(BR / R**2), sp.csr_matrix(L2_bf), format='csr')
    Eint_bf = sp.kron(sp.identity(nR), sp.diags(E_rot_bf), format='csr')
    V_mat_bf = sp.kron(sp.diags(V_scale), sp.csr_matrix(C_bf), format='csr')

    H_bf = T_R_bf + Cent_bf + Eint_bf + V_mat_bf
    A_bf = (H_bf - Ec * sp.identity(H_bf.shape[0])).astype(complex).tolil()
    b_bf = np.zeros(nR * n_ch, dtype=complex)

    for p in range(n_ch):
        A_bf[p, :] = 0.0; A_bf[p, p] = 1.0; b_bf[p] = 0.0

    K_mat_bf = U @ np.diag(K_sf) @ U.T
    source_bf = U @ source_sf

    for p in range(n_ch):
        row = j_end * n_ch + p
        A_bf[row, :] = 0.0
        for q in range(n_ch):
            A_bf[row, j_end * n_ch + q] = -1j * K_mat_bf[p, q]
        A_bf[row, j_end * n_ch + p] += coef[0]
        A_bf[row, (j_end - 1) * n_ch + p] = coef[1]
        A_bf[row, (j_end - 2) * n_ch + p] = coef[2]
        b_bf[row] = source_bf[p]

    chi_bf = spla.spsolve(A_bf.tocsc(), b_bf).reshape(nR, n_ch)

    # 3. Transform BF to SF
    chi_bf_to_sf = chi_bf @ U

    # Extract S-matrix in SF
    S_sf = np.zeros(n_ch, dtype=complex)
    for p in range(n_ch):
        if K_sf[p].imag == 0:
            if p == 0:
                S_sf[p] = (chi_sf[-1, p] - np.exp(-1j * K0 * R[-1])) / np.exp(1j * K0 * R[-1])
            else:
                S_sf[p] = chi_sf[-1, p] / np.exp(1j * K_sf[p].real * R[-1])

    # ------------------ PLOTTING ------------------
    fig = plt.figure(figsize=(16, 11))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.1, 1.0], hspace=0.32, wspace=0.28)

    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']

    # Panel 1: Space-Fixed (SF) Wavefunctions chi_SF(R)
    ax1 = fig.add_subplot(gs[0, 0])
    for ch_idx in range(n_ch):
        j_val, l_val = sf_ch[ch_idx]
        ax1.plot(R, np.abs(chi_sf[:, ch_idx]), color=colors[ch_idx], lw=1.6,
                 label=rf'$|j={j_val}, l={l_val}\rangle$')
    ax1.set_title(r'(a) Space-Fixed (SF) Wavefunctions $|\chi_{\mathrm{SF}}(R)|$', fontsize=11, fontweight='bold')
    ax1.set_xlabel(r'$R$ (Å)', fontsize=10)
    ax1.set_ylabel(r'Channel Amplitude $|\chi|$', fontsize=10)
    ax1.legend(fontsize=9, loc='upper right')
    ax1.grid(True, ls=':', alpha=0.6)
    ax1.set_xlim(R_min, R_max)

    # Panel 2: Body-Fixed (BF) Wavefunctions chi_BF(R)
    ax2 = fig.add_subplot(gs[0, 1])
    for ch_idx in range(n_ch):
        j_val, K_val = bf_ch[ch_idx]
        ax2.plot(R, np.abs(chi_bf[:, ch_idx]), color=colors[ch_idx], lw=1.6,
                 label=rf'$|j={j_val}, K={K_val}\rangle$')
    ax2.set_title(r'(b) Body-Fixed (BF) Wavefunctions $|\chi_{\mathrm{BF}}(R)|$', fontsize=11, fontweight='bold')
    ax2.set_xlabel(r'$R$ (Å)', fontsize=10)
    ax2.set_ylabel(r'Channel Amplitude $|\chi|$', fontsize=10)
    ax2.legend(fontsize=9, loc='upper right')
    ax2.grid(True, ls=':', alpha=0.6)
    ax2.set_xlim(R_min, R_max)

    # Panel 3: Exact Frame Equivalence Proof |U^T chi_BF - chi_SF|
    ax3 = fig.add_subplot(gs[0, 2])
    diff_channels = np.abs(chi_bf_to_sf - chi_sf)
    for ch_idx in range(n_ch):
        j_val, l_val = sf_ch[ch_idx]
        lbl = r'$|U^\dagger \chi_{\mathrm{BF}} - \chi_{\mathrm{SF}}|_{(j=%d, l=%d)}$' % (j_val, l_val)
        ax3.semilogy(R, np.maximum(diff_channels[:, ch_idx], 1e-16), color=colors[ch_idx], lw=1.5,
                     label=lbl)
    ax3.set_title(r'(c) Frame Transformation Residual (Exact)', fontsize=11, fontweight='bold')
    ax3.set_xlabel(r'$R$ (Å)', fontsize=10)
    ax3.set_ylabel(r'Pointwise Absolute Difference', fontsize=10)
    ax3.axhline(1e-11, ls='--', color='gray', lw=1.0, label=r'Machine Tolerance $\sim 10^{-11}$')
    ax3.legend(fontsize=8, loc='center right')
    ax3.grid(True, ls=':', alpha=0.6)
    ax3.set_xlim(R_min, R_max)
    ax3.set_ylim(1e-15, 1e-9)

    # Panel 4: Unitary Transformation Matrix U Heatmap
    ax4 = fig.add_subplot(gs[1, 0])
    im4 = ax4.imshow(U, cmap='coolwarm', vmin=-1.0, vmax=1.0)
    ax4.set_title(r'(d) Unitary Matrix $U_{K, l}(J=1, j, p=-1)$', fontsize=11, fontweight='bold')
    ax4.set_xticks(range(n_ch))
    ax4.set_yticks(range(n_ch))
    x_tick_labels = [r'$l=%d$' % l + '\n' + r'$(j=%d)$' % j for (j, l) in sf_ch]
    y_tick_labels = [r'$K=%d$' % K + ' ' + r'$(j=%d)$' % j for (j, K) in bf_ch]
    ax4.set_xticklabels(x_tick_labels, fontsize=8)
    ax4.set_yticklabels(y_tick_labels, fontsize=8)
    for i in range(n_ch):
        for j in range(n_ch):
            val = U[i, j]
            if abs(val) > 1e-4:
                ax4.text(j, i, f'{val:+.3f}', ha='center', va='center',
                         color='black' if abs(val) < 0.65 else 'white', fontsize=8, fontweight='bold')
    fig.colorbar(im4, ax=ax4, fraction=0.046, pad=0.04)

    # Panel 5: Coriolis Coupling L^2 Matrix in BF Frame
    ax5 = fig.add_subplot(gs[1, 1])
    im5 = ax5.imshow(L2_bf, cmap='viridis')
    ax5.set_title(r'(e) Coriolis Operator $\hat{L}^2 = (\mathbf{J}-\mathbf{j})^2$ in BF', fontsize=11, fontweight='bold')
    ax5.set_xticks(range(n_ch))
    ax5.set_yticks(range(n_ch))
    x_bf_ticks = [r'$K=%d$' % K + '\n' + r'$(j=%d)$' % j for (j, K) in bf_ch]
    y_bf_ticks = [r'$K=%d$' % K + ' ' + r'$(j=%d)$' % j for (j, K) in bf_ch]
    ax5.set_xticklabels(x_bf_ticks, fontsize=8)
    ax5.set_yticklabels(y_bf_ticks, fontsize=8)
    for i in range(n_ch):
        for j in range(n_ch):
            val = L2_bf[i, j]
            if abs(val) > 1e-4:
                ax5.text(j, i, f'{val:.2f}', ha='center', va='center',
                         color='white' if val > 2.5 else 'black', fontsize=8, fontweight='bold')
    fig.colorbar(im5, ax=ax5, fraction=0.046, pad=0.04)

    # Panel 6: S-matrix Channel Probabilities P_{j,l} = |S_{j,l}|^2
    ax6 = fig.add_subplot(gs[1, 2])
    labels_sf = [rf'$(j={j}, l={l})$' for (j, l) in sf_ch]
    probs = np.abs(S_sf)**2
    # Flux weighting for open channels: v_rel = K / mu_R
    flux_weights = np.array([K_sf[p].real / K0 for p in range(n_ch)])
    probs_flux = probs * flux_weights
    
    bars = ax6.bar(range(n_ch), probs_flux, color=colors, edgecolor='black', width=0.55)
    ax6.set_title(r'(f) Open-Channel Probabilities $P_{(j, l)} = |S_{(j, l)}|^2 \frac{K_j}{K_0}$', fontsize=11, fontweight='bold')
    ax6.set_xticks(range(n_ch))
    ax6.set_xticklabels(labels_sf, rotation=25, ha='right', fontsize=9)
    ax6.set_ylabel(r'Scattering Probability', fontsize=10)
    ax6.grid(True, axis='y', ls=':', alpha=0.6)
    for bar, val in zip(bars, probs_flux):
        ax6.text(bar.get_x() + bar.get_width()/2.0, val + 0.02, f'{val*100:.1f}%',
                 ha='center', va='bottom', fontsize=9, fontweight='bold')
    ax6.set_ylim(0, 1.05)

    fig.suptitle(r'Full Quantum Reactive Scattering with Non-Zero Angular Momentum $J = 1$ ($p = -1, E_c = 400\ \mathrm{cm}^{-1}$)' + '\n' +
                 r'Rigorous Demonstration of Body-Fixed (BF) $\leftrightarrow$ Space-Fixed (SF) Representations & Exact Unitary Frame Transformation',
                 fontsize=13, fontweight='bold', y=0.98)

    out_path = ROOT / 'article' / 'fig' / 'fig_j1_scattering_frames.png'
    out_path.parent.mkdir(exist_ok=True)
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"Successfully generated figure: {out_path}")
    return out_path


if __name__ == '__main__':
    generate_j1_scatter_figure()
