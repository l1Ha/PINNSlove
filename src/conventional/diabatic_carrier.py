"""Diabatic (nonadiabatic) channel carrier machinery for the 2D A+BC wave.

Theory: docs/theory_scattering_2d.md section 12.

The adiabatic eigensystem Phi_j(R, gamma) carries derivative couplings
A_jk = <Phi_j | d_R Phi_k> that dominate the carrier residual (theory
sections 8 and 10).  The parallel-transport (diabatic) transformation

    dU/dR = -A(R) U,   U(R_lo) = I,   Phi_diag = Phi U

annihilates the derivative coupling exactly (A_tilde = U^T A U + U^T d_R U = 0)
and turns the coupled radial equation into a purely multiplicative one

    -B_R w'' + Veff(R) w = E w,
    Veff = U^T diag(eps) U - B_R * U^T A^T A U,                     (Eq. 30)

where the second (O(A^2)) term is the exact discretisation of the
second-derivative projection (it is usually dropped in the literature; here
it is kept and its norm is reported).

This module provides:
  * derivative_coupling(Phi, h)          -- Eq. 31 with an antisymmetry check
  * integrate_diabatic_transform(A, h)   -- Eq. 32, 2nd-order Magnus (expm),
                                            with an orthonormality check
  * DiabaticBasis                        -- Phi_diag, Veff curves, diagnostics
  * build_langer_carriers(curve, ...)    -- the uniform Langer-Airy carrier of
                                            theory section 7.3 with the
                                            discrete-dispersion phase of
                                            section 9.2, applied to ANY
                                            channel curve (adiabatic eps_j or
                                            diabatic Veff_jj).  NOTE: the
                                            formulas duplicate those inside
                                            AdiabaticCarrier2D (same theory
                                            sections); consolidating the two
                                            implementations is a pending
                                            refactor tracked in the docs.
"""
from __future__ import annotations

import numpy as np
from scipy.linalg import expm
from scipy.special import airy


# ---------------------------------------------------------------------------
def derivative_coupling(Phi, h):
    """A_jk(R_i) = <Phi_j | (Phi_k(R_{i+1}) - Phi_k(R_{i-1})) / (2h)>  (Eq. 31).

    Returns (A, asym) with A (n_R, n_gamma, n_gamma) and asym =
    max|A + A^T| (antisymmetry check; should be at the level of the central
    difference truncation error times ||Phi||).
    """
    n_R, n_gamma, _ = Phi.shape
    dPhi = np.empty_like(Phi)
    dPhi[1:-1] = (Phi[2:] - Phi[:-2]) / (2.0 * h)
    dPhi[0] = (Phi[1] - Phi[0]) / h
    dPhi[-1] = (Phi[-1] - Phi[-2]) / h
    A = np.empty((n_R, n_gamma, n_gamma))
    for i in range(n_R):
        A[i] = Phi[i].T @ dPhi[i]
    asym = float(np.max(np.abs(A + np.transpose(A, (0, 2, 1)))))
    return A, asym


# ---------------------------------------------------------------------------
def derivative_coupling_hf(Phi, eps, H_gamma, h, gap_floor=1.0):
    """Hellmann-Feynman-form derivative coupling (fix for the mode-mixing
    failure of the central-difference A, see theory section 13):

        A_jk(R_i) = <Phi_j | d_R H_gamma | Phi_k> / (eps_k - eps_j),  j != k
        A_jj = 0 (real basis; theory Eq. 24)

    with d_R H_gamma obtained by central differences of the node matrices
    (H_gamma is a smooth function of R).  ANTI-SYMMETRY IS EXACT BY
    CONSTRUCTION (H_gamma symmetric, the 1/(eps_k - eps_j) factor flips
    sign).  Pairs with |eps_k - eps_j| < gap_floor (near-degenerate, where
    the adiabatic basis itself is ill-defined) are set to zero and counted.
    Returns (A, n_zeroed, max_gap_ratio).
    """
    n_R, n, _ = Phi.shape
    W = np.empty_like(H_gamma)
    W[1:-1] = (H_gamma[2:] - H_gamma[:-2]) / (2.0 * h)
    W[0] = (H_gamma[1] - H_gamma[0]) / h
    W[-1] = (H_gamma[-1] - H_gamma[-2]) / h
    A = np.zeros((n_R, n, n))
    n_zeroed = 0
    worst = 0.0
    for i in range(n_R):
        Wij = Phi[i].T @ W[i] @ Phi[i]          # symmetric by construction
        gap = eps[i][None, :] - eps[i][:, None]
        small = np.abs(gap) < gap_floor      # includes the diagonal (gap=0)
        n_zeroed += int(np.triu(small, 1).sum())
        denom = np.where(small, 1.0, gap)
        M = Wij / denom
        M[small] = 0.0                       # diagonal forced to zero (Eq. 24)
        A[i] = M
        off = np.abs(np.triu(M, 1))
        if off.size:
            worst = max(worst, float(off.max()))
    return A, n_zeroed, worst


# ---------------------------------------------------------------------------
def integrate_diabatic_transform(A, h):
    """U_{i+1} = expm(-h A(R_{i+1/2})) U_i, U_0 = I   (Eq. 32, 2nd-order Magnus).

    Returns (U, orth_err) with U (n_R, n_gamma, n_gamma) and
    orth_err = max_i ||U_i^T U_i - I||_max.
    """
    n_R, n, _ = A.shape
    U = np.empty((n_R, n, n))
    U[0] = np.eye(n)
    for i in range(n_R - 1):
        A_mid = 0.5 * (A[i] + A[i + 1])
        U[i + 1] = expm(-h * A_mid) @ U[i]
    orth = float(max(np.max(np.abs(U[i].T @ U[i] - np.eye(n)))
                     for i in range(n_R)))
    return U, orth


# ---------------------------------------------------------------------------
class DiabaticBasis:
    """Diabatic channel functions, potential matrix and diagnostics."""

    def __init__(self, Phi, eps, h, B_R, H_gamma=None, gap_floor=1.0):
        self.h = float(h)
        self.B_R = float(B_R)
        self.A_fd, self.asym = derivative_coupling(Phi, h)   # diagnostic only
        if H_gamma is not None:
            # Hellmann-Feynman coupling: antisymmetry exact by construction;
            # the central-difference A fails its own antisymmetry check in
            # the well (mode mixing among near-degenerate eigenvectors)
            self.A, self.n_gap_zeroed, self.max_offdiag_A = \
                derivative_coupling_hf(Phi, eps, H_gamma, h,
                                       gap_floor=gap_floor)
        else:
            self.A, self.n_gap_zeroed, self.max_offdiag_A = self.A_fd, 0, None
        self.U, self.orth_err = integrate_diabatic_transform(self.A, h)
        n_R, n_gamma, _ = Phi.shape
        self.Phi_dia = np.empty_like(Phi)
        self.Vdia = np.empty((n_R, n_gamma, n_gamma))
        self.Veff = np.empty((n_R, n_gamma, n_gamma))
        for i in range(n_R):
            self.Phi_dia[i] = Phi[i] @ self.U[i]
            self.Vdia[i] = self.U[i].T @ (eps[i][:, None] * self.U[i])
            G = self.U[i].T @ (self.A[i].T @ self.A[i]) @ self.U[i]
            self.Veff[i] = self.Vdia[i] - self.B_R * G
        scale = float(np.max(np.abs(self.Vdia)))
        self.offdiag_ratio = float(np.max(np.abs(self.Veff
                                                 - np.diagonal(self.Veff,
                                                               axis1=1,
                                                               axis2=2)
                                                 [:, :, None]))) / scale

    def curves(self):
        """Diagonal diabatic curves Veff_jj(R): (n_R, n_gamma)."""
        return np.diagonal(self.Veff, axis1=1, axis2=2).copy()


# ---------------------------------------------------------------------------
def dispersion_q(k2, h):
    """Discrete-dispersion wavenumber (theory Eq. 21):

        30 - 32 cos(qh) + 2 cos(2qh) = 12 h^2 k2,

    unique root in (0, pi/h) since f'(x) = 4 sin x (8 - cos x) > 0.
    Only meaningful for k2 > 0; callers use the continuous decay constant on
    the forbidden side."""
    tgt = 12.0 * h ** 2 * np.asarray(k2, float)
    x = np.empty_like(tgt)
    for idx, tv in np.ndenumerate(tgt):
        lo, hi = 1e-14, np.pi - 1e-14
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            f = 30.0 - 32.0 * np.cos(mid) + 2.0 * np.cos(2.0 * mid)
            if f < tv:
                lo = mid
            else:
                hi = mid
        x[idx] = 0.5 * (lo + hi)
    return x / h


# ---------------------------------------------------------------------------
def build_langer_carriers(curve, R, E_tot, B_R, norm_window_frac=0.25):
    """Uniform Langer-Airy carrier (theory section 7.3) on ONE channel curve,
    with the discrete-dispersion phase (section 9.2) and unit asymptotic
    e^{+/- i k R} normalisation (chi = R psi has no 1/sqrt(R) factor).

    curve: channel potential V_j(R); turning point = last entry into the
    classically allowed region.  Returns dict with U_plus / U_minus /
    dU_plus / dU_minus / i_c, or None if the curve is forbidden everywhere.
    """
    R = np.asarray(R, float)
    n = len(R)
    h = R[1] - R[0]
    k2 = (E_tot - np.asarray(curve, float)) / B_R
    allow = k2 > 0.0
    if not allow.any():
        return None
    signchg = np.where(np.diff(allow.astype(int)) != 0)[0]
    i_c = 0
    for isc in signchg:
        if allow[isc + 1]:
            i_c = isc + 1
    if allow[0]:
        i_c = 0
    # phase step: discrete-dispersion q in the allowed region, continuous
    # decay constant in the forbidden region
    g = np.where(allow, dispersion_q(np.abs(k2), h) * h,
                 np.sqrt(np.abs(k2)) * h)
    S = np.empty(n)
    S[i_c] = 0.0
    acc = 0.0
    for i in range(i_c + 1, n):
        acc += 0.5 * (g[i] + g[i - 1])
        S[i] = acc
    acc = 0.0
    for i in range(i_c - 1, -1, -1):
        acc += 0.5 * (g[i] + g[i + 1])
        S[i] = -acc
    zeta = np.sign(R - R[i_c]) * (1.5 * np.abs(S)) ** (2.0 / 3.0)
    zeta[np.abs(R - R[i_c]) < 1e-12] = 0.0
    amp = np.empty(n)
    deps = np.gradient(curve, h)
    for i in range(n):
        if i == i_c:
            amp[i] = (B_R / max(abs(deps[i]), 1e-30)) ** (1.0 / 6.0)
        else:
            amp[i] = (np.abs(zeta[i]) / max(abs(k2[i]), 1e-300)) ** 0.25
    Ai, _, Bi, _ = airy(-zeta)
    U_p = amp * (1j * Ai + Bi)          # pure incoming branch
    U_m = amp * (-1j * Ai + Bi)         # pure outgoing branch
    k_asym = np.sqrt(max((E_tot - curve[-1]) / B_R, 0.0))
    if k_asym > 0:
        win = R > (R[-1] - norm_window_frac * (R[-1] - R[0]))
        M = np.column_stack([np.exp(1j * k_asym * R[win]),
                             np.exp(-1j * k_asym * R[win])])
        cp = np.linalg.lstsq(M, U_p[win], rcond=None)[0]
        cm = np.linalg.lstsq(M, U_m[win], rcond=None)[0]
        if abs(cp[0]) > 1e-12:
            U_p = U_p / cp[0]
        if abs(cm[1]) > 1e-12:
            U_m = U_m / cm[1]
    d = lambda f: np.gradient(f, h, edge_order=2)
    return dict(U_plus=U_p, U_minus=U_m, dU_plus=d(U_p), dU_minus=d(U_m),
                i_c=int(i_c), curve=np.asarray(curve, float))
