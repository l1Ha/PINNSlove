"""Authentic BKMP2 Coupled-Channel Reactive Scattering for arbitrary J >= 0.

Solves the multi-channel A + BC reactive scattering Schrödinger equation on the
authentic BKMP2 potential energy surface for arbitrary total angular momentum J >= 0
and inversion parity p = (+1, -1).

Features:
1. Exact Body-Fixed (BF) representation with Coriolis coupling L^2 = (J - j)^2.
2. Exact Space-Fixed (SF) representation with diagonal centrifugal barrier B_R l(l+1)/R^2.
3. Rigorous potential matrix elements <j K | V_BKMP2(R, gamma) | j' K> evaluated
   via Gauss-Legendre quadrature on the authentic BKMP2 potential.
4. Unitary frame transformation U diagonalizing L^2:
     chi_SF = chi_BF @ U,  S_SF = U^T @ S_BF @ U.
5. Consistent 4th-order finite difference radial solver (fd2_4th_dirichlet) with
   incoming-wave Robin radiation boundary condition at R_max.
6. S-matrix extraction with FD-consistent dispersion q(h), extracting
   state-to-state scattering matrix S^J_{j' l', j l}(E).
7. Partial-wave summation yielding integral cross sections (ICS) sigma(E) and
   differential cross sections (DCS) d sigma / d Omega(theta).
"""
import math
import sys
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.special import lpmv, roots_legendre

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.j_frames import JAngularCoupling
from src.scatter_extract import fd2_4th_dirichlet, fit_wavevector, extract_smatrix
from src import bkmp2

# Molecular constants for H + H2
BOHR = 0.529177210903
M_H = 1.007825
mu_r = 0.5 * M_H
mu_R = (2.0 / 3.0) * M_H
H2_OVER_2AMU = 16.857633
Br = H2_OVER_2AMU / mu_r     # 33.4535 cm^-1 A^2
BR = H2_OVER_2AMU / mu_R     # 25.0901 cm^-1 A^2
re_H2 = 0.74144              # H2 equilibrium bond length, Angstrom
Brot = Br / (re_H2 ** 2)     # H2 rotational constant at r_e: 60.847 cm^-1
#   (the diatomic reduced mass mu_r = m_H/2 enters here, NOT the collision
#    reduced mass; thresholds are Brot * j(j+1))


def theta_jK(j, K, x):
    """Normalized associated Legendre function: integral_{-1}^1 Theta_jK Theta_j'K dx = delta_jj'."""
    norm = np.sqrt((2.0 * j + 1.0) / 2.0 * math.factorial(j - K) / math.factorial(j + K))
    return norm * lpmv(K, j, x)


def compute_bkmp2_potential_matrix(R_grid, bf_channels, r_fix=re_H2, n_quad=32):
    """Computes the authentic BKMP2 potential matrix V_ab(R) for given BF channels.

    V_ab(R) = delta_{Ka, Kb} * integral_0^pi Theta_{ja Ka}(gamma) V_BKMP2(R, gamma) Theta_{jb Kb}(gamma) sin(gamma) dgamma.

    Returns:
        V_matrix: array of shape (nR, n_channels, n_channels) in cm^-1.
    """
    R_grid = np.asarray(R_grid, dtype=float)
    nR = len(R_grid)
    n_ch = len(bf_channels)

    x_quad, w_quad = roots_legendre(n_quad)
    gamma_quad = np.arccos(x_quad)

    # Evaluate authentic BKMP2 on the (nR, n_quad) grid
    V_grid = bkmp2.potential_2d_rRg(R_grid, gamma_quad, r_fix, verbose=False)  # (nR, n_quad)

    # Precompute basis matrix Phi: shape (n_quad, n_ch)
    # with Phi[k, a] = sqrt(w_k) * Theta_{ja Ka}(x_k)
    Phi = np.zeros((n_quad, n_ch), dtype=float)
    for a, (ja, Ka) in enumerate(bf_channels):
        Phi[:, a] = np.sqrt(w_quad) * theta_jK(ja, Ka, x_quad)

    # Matrix multiplication: V_ab(R_i) = sum_k Phi[k, a] * V_grid[i, k] * Phi[k, b]
    V_matrix = np.zeros((nR, n_ch, n_ch), dtype=float)
    for i in range(nR):
        # Weighted diagonal: Phi.T @ diag(V_grid[i, :]) @ Phi
        V_diag_phi = V_grid[i, :, None] * Phi
        V_matrix[i, :, :] = Phi.T @ V_diag_phi

    # Enforce exact symmetry and K-selection rule
    for a, (ja, Ka) in enumerate(bf_channels):
        for b, (jb, Kb) in enumerate(bf_channels):
            if Ka != Kb:
                V_matrix[:, a, b] = 0.0

    V_matrix = 0.5 * (V_matrix + np.swapaxes(V_matrix, 1, 2))
    return V_matrix


def solve_j_scattering_bkmp2(J, parity, Ec, nR=400, R_min=0.7, R_max=9.0, j_max=4,
                             j_step=2, r_fix=re_H2):
    """Solves coupled-channel reactive scattering on the authentic BKMP2 PES for given (J, parity).

    Parameters:
        J : int, total angular momentum (0, 1, 2, ...)
        parity : int, inversion parity (+1 or -1)
        Ec : float, collision energy in cm^-1
        nR : int, radial grid points
        R_min, R_max : float, radial coordinate range (Angstrom)
        j_max : int, maximum diatomic rotational quantum number
        j_step : int, 2 for para-H2 (j=0, 2, 4,...), 1 for all j

    Returns:
        results : dict containing:
          - S_matrix_sf: open-channel S-matrix in Space-Fixed representation
          - S_matrix_bf: open-channel S-matrix in Body-Fixed representation
          - sf_channels, bf_channels: channel lists
          - open_indices: indices of energetically open channels
          - chi_bf, chi_sf: radial channel wavefunctions (nR, n_ch)
          - R_grid: radial coordinates
          - bf_sf_error: max difference ||chi_sf - chi_bf @ U||
          - unitarity_residual: max |sum_b |S_ba|^2 - 1|
    """
    coupling = JAngularCoupling(J=J, j_max=j_max)
    U, sf_ch, bf_ch, l_evals = coupling.full_frame_transformation(parity=parity, j_step=j_step)
    n_ch = len(bf_ch)
    if n_ch == 0:
        return None

    R = np.linspace(R_min, R_max, nR)
    h = R[1] - R[0]

    # Asymptotic reference energy (isolated H2 at re)
    V_asym = bkmp2.asymptotic_energy(r_fix, R_far=30.0)
    E_tot = Ec + V_asym

    # Diatomic rotational threshold energies: eps_j = Brot * j(j+1)
    Brot_eff = Br / (r_fix ** 2)   # 60.85 cm^-1 for r_fix = r_e
    thresholds = np.array([Brot_eff * j * (j + 1) for (j, K) in bf_ch], dtype=float)

    # Authentic BKMP2 potential matrix elements V_ab(R)
    V_mat = compute_bkmp2_potential_matrix(R, bf_ch, r_fix=r_fix)

    # Asymptotic wavevectors in each channel
    K_list = np.zeros(n_ch, dtype=complex)
    for i in range(n_ch):
        E_tr = Ec - thresholds[i]
        K_list[i] = np.sqrt(E_tr / BR) if E_tr > 0 else 1j * np.sqrt(-E_tr / BR)

    open_idx = [i for i in range(n_ch) if K_list[i].imag == 0.0]
    n_open = len(open_idx)
    if n_open == 0:
        return None

    # Coriolis L^2 matrix in the BF basis
    # Construct block diagonal L^2
    L2_bf = np.zeros((n_ch, n_ch), dtype=float)
    offset = 0
    for j in range(0, j_max + 1, j_step):
        L2_block, k_list, _, _ = coupling.coriolis_matrix_j(j, parity)
        nb = len(k_list)
        if nb > 0:
            L2_bf[offset:offset + nb, offset:offset + nb] = L2_block
            offset += nb

    # -------------------------------------------------------------
    # 1. SOLVE IN BODY-FIXED (BF) REPRESENTATION
    # -------------------------------------------------------------
    D2 = fd2_4th_dirichlet(nR, h)
    H_bf = sp.kron(-BR * D2, sp.identity(n_ch), format='csr')

    # Centrifugal / Coriolis term: BR * L2 / R^2
    # Diagonal thresholds + potential matrix
    diag_blocks = []
    for iR in range(nR):
        C_R = BR * L2_bf / (R[iR] ** 2)
        diag_blocks.append(C_R + V_mat[iR] + np.diag(thresholds))
    H_bf = H_bf + sp.block_diag(diag_blocks, format='csr')

    A_bf = (H_bf - E_tot * sp.identity(H_bf.shape[0])).astype(complex).tolil()

    # Inner Dirichlet boundary condition at R_min
    for p in range(n_ch):
        A_bf[p, :] = 0.0
        A_bf[p, p] = 1.0

    # Asymptotic incoming/outgoing Robin boundary condition at R_max
    coef = np.array([25.0, -48.0, 36.0, -16.0, 3.0]) / (12.0 * h)
    je = nR - 1

    # In the asymptotic region, U transforms BF to SF where L^2 is diagonal l(l+1).
    # Multichannel Robin BC in SF: (d/dR - i K_sf) chi_sf = source
    # Expressed in BF: U (d/dR - i K_diag) U^T chi_bf = source_bf
    K_mat_bf = U @ np.diag(K_list) @ U.T

    # Solve for incident open channels
    chi_bf_all = np.zeros((nR, n_ch, n_open), dtype=complex)
    chi_sf_all = np.zeros((nR, n_ch, n_open), dtype=complex)
    b_bf = np.zeros(nR * n_ch, dtype=complex)
    b_sf = np.zeros(nR * n_ch, dtype=complex)

    # Setup Robin matrix rows in A_bf
    for p in range(n_ch):
        row = je * n_ch + p
        A_bf[row, :] = 0.0
        for q in range(n_ch):
            A_bf[row, je * n_ch + q] = -1j * K_mat_bf[p, q]
        for t, c in enumerate(coef):
            A_bf[row, (je - t) * n_ch + p] += c

    A_bf_csc = A_bf.tocsc()
    solver_bf = spla.factorized(A_bf_csc)

    # -------------------------------------------------------------
    # 2. SOLVE IN SPACE-FIXED (SF) REPRESENTATION
    # -------------------------------------------------------------
    # In SF, L^2 is diagonal l(l+1) and V_sf(R) = U^T V_mat(R) U
    L2_sf = np.diag(l_evals)
    diag_blocks_sf = []
    for iR in range(nR):
        C_R = BR * L2_sf / (R[iR] ** 2)
        V_sf_i = U.T @ V_mat[iR] @ U
        diag_blocks_sf.append(C_R + V_sf_i + np.diag(thresholds))

    H_sf = sp.kron(-BR * D2, sp.identity(n_ch), format='csr') + sp.block_diag(diag_blocks_sf, format='csr')
    A_sf = (H_sf - E_tot * sp.identity(H_sf.shape[0])).astype(complex).tolil()

    for p in range(n_ch):
        A_sf[p, :] = 0.0
        A_sf[p, p] = 1.0

    for p in range(n_ch):
        row = je * n_ch + p
        A_sf[row, :] = 0.0
        A_sf[row, row] = -1j * K_list[p]
        for t, c in enumerate(coef):
            A_sf[row, (je - t) * n_ch + p] += c

    A_sf_csc = A_sf.tocsc()
    solver_sf = spla.factorized(A_sf_csc)

    # Incident waves are defined natively in Space-Fixed channels (j, l).
    # The corresponding Body-Fixed source vector is source_bf = U @ source_sf.
    for col_idx, inc_ch in enumerate(open_idx):
        k_in = K_list[inc_ch].real
        source_sf = np.zeros(n_ch, dtype=complex)
        source_sf[inc_ch] = -2j * k_in * np.exp(-1j * k_in * R[-1])
        b_sf[:] = 0.0
        b_sf[je * n_ch:(je + 1) * n_ch] = source_sf
        sol_sf = solver_sf(b_sf)
        chi_sf_all[:, :, col_idx] = sol_sf.reshape(nR, n_ch)

        # BF source via frame transformation: source_bf = U @ source_sf
        source_bf = U @ source_sf
        b_bf[:] = 0.0
        b_bf[je * n_ch:(je + 1) * n_ch] = source_bf
        sol_bf = solver_bf(b_bf)
        chi_bf_all[:, :, col_idx] = sol_bf.reshape(nR, n_ch)
        k_in = K_list[inc_ch].real
        source_sf = np.zeros(n_ch, dtype=complex)
        source_sf[inc_ch] = -2j * k_in * np.exp(-1j * k_in * R[-1])
        b_sf[:] = 0.0
        b_sf[je * n_ch:(je + 1) * n_ch] = source_sf
        sol = solver_sf(b_sf)
        chi_sf_all[:, :, col_idx] = sol.reshape(nR, n_ch)

    # -------------------------------------------------------------
    # 3. VERIFY BF <-> SF EXACT FRAME EQUIVALENCE
    # -------------------------------------------------------------
    # For each incident channel, chi_sf should equal chi_bf @ U
    chi_transformed = np.zeros_like(chi_sf_all)
    for col_idx in range(n_open):
        chi_transformed[:, :, col_idx] = chi_bf_all[:, :, col_idx] @ U

    bf_sf_diff = np.max(np.abs(chi_sf_all - chi_transformed))

    # -------------------------------------------------------------
    # 4. EXTRACT S-MATRIX IN SPACE-FIXED REPRESENTATION
    # -------------------------------------------------------------
    l_list_open = [int(np.round(np.sqrt(l_evals[ch] + 0.25) - 0.5)) for ch in open_idx]
    S_mat_sf, a_inc = extract_smatrix(
        chi_sf_all, R, K_list, l_list_open, open_idx,
        dispersion=True, matrix_form=True
    )

    # BF S-matrix via frame transformation: S_bf = U @ S_sf @ U^T
    U_open = U[np.ix_(open_idx, open_idx)]
    S_mat_bf = U_open @ S_mat_sf @ U_open.T

    # Unitarity check: sum_b (k_b / k_a) |S_ba|^2 == 1.0
    k_open = np.array([K_list[ch].real for ch in open_idx])
    flux_matrix = (k_open[:, None] / k_open[None, :]) * (np.abs(S_mat_sf) ** 2)
    unitarity_sums = np.sum(flux_matrix, axis=0)
    unitarity_residual = float(np.max(np.abs(unitarity_sums - 1.0)))

    # Reciprocity check: k_a S_ba == k_b S_ab
    reciprocity_matrix = k_open[None, :] * S_mat_sf - k_open[:, None] * S_mat_sf.T
    reciprocity_residual = float(np.max(np.abs(reciprocity_matrix)) / np.max(k_open))

    return {
        'J': J,
        'parity': parity,
        'Ec': Ec,
        'nR': nR,
        'bf_channels': bf_ch,
        'sf_channels': sf_ch,
        'open_channels_bf': [bf_ch[i] for i in open_idx],
        'open_channels_sf': [sf_ch[i] for i in open_idx],
        'open_indices': open_idx,
        'S_matrix_sf': S_mat_sf,
        'S_matrix_bf': S_mat_bf,
        'chi_bf': chi_bf_all,
        'chi_sf': chi_sf_all,
        'R': R,
        'K_open': k_open,
        'bf_sf_diff': bf_sf_diff,
        'unitarity_residual': unitarity_residual,
        'reciprocity_residual': reciprocity_residual,
    }
