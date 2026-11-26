"""Adiabatic-channel uniform Langer-Airy carrier for the 2D A+BC scattered wave.

Theory: docs/theory_scattering_2d.md section 7.

    H = -B_R d2_R + C(R) j2 + V(R,gamma),  C = B_R/R^2 + B_r/r_e^2

1. Adiabatic channel basis: at each R diagonalize H_gamma(R) = C(R) j2 +
   diag(V(R,gamma_g)) on the Legendre-DVR nodes;
       H_gamma(R) Phi_j(R,.) = eps_j(R) Phi_j(R,.),
   tracking channels across R by maximum eigenvector overlap (no sorting
   swaps).
2. Per-channel eikonal: k_j(R) = sqrt((E - eps_j(R))/B_R); turning point
   R_c,j from eps_j(R_c) = E; action S_j(R) = sign(R-R_c) int_{R_c}^{R} k
   (allowance side) / |k| (forbidden side) integrated to the turning point.
3. Uniform Langer-Airy carrier (zero-free on the real axis, pure
   incoming/outgoing asymptotics):
       zeta_j = sign(R-R_c,j)(3|S_j|/2)^{2/3},
       A_j    = |zeta_j / (B_R k_j^2)|^{1/4}  (regularized at R_c with
                (B_R/|d eps_j/dR|)^{1/6}),
       U^pm_j = A_j [ +/- i Ai(-zeta_j) + Bi(-zeta_j) ].

Zero fitting parameters: everything is computed from the PES by the above
formulas.  Channels with no classically allowed region in [R_lo, R_hi] are
dropped (they are closed everywhere; amplitude < e^{-40} at R_hi).
"""
from __future__ import annotations

import numpy as np


class AdiabaticCarrier2D:
    """Adiabatic-channel Langer-Airy carrier on a given (R, gamma_DVR) grid."""

    def __init__(self, pes, E_c, R, xg, j2, sig_n=1e-3):
        self.pes = pes
        self.E_c = float(E_c)
        self.E_tot = E_c + float(pes.V_valence(pes.re, pes.re, pes.theta_e)) + pes.De
        self.R = np.asarray(R, float)
        self.xg = np.asarray(xg, float)
        self.gam = np.arccos(self.xg)
        self.j2 = np.asarray(j2)
        self.n_R, self.n_gamma = len(self.R), len(self.gam)
        self.h = self.R[1] - self.R[0]

        # ---- adiabatic channel curves -----------------------------------
        RR, GG = np.meshgrid(self.R, self.gam, indexing='ij')
        Vg = pes.V_jacobi(pes.re, RR, GG)                     # (n_R, n_gamma)
        C = pes.BR / self.R ** 2 + pes.Br / pes.re ** 2       # (n_R,)
        from scipy.optimize import linear_sum_assignment
        self.eps = np.empty((self.n_R, self.n_gamma))
        self.Phi = np.empty((self.n_R, self.n_gamma, self.n_gamma))
        prev = None
        for i in range(self.n_R):
            M = C[i] * self.j2 + np.diag(Vg[i])
            ev, evec = np.linalg.eigh(M)
            if prev is not None:
                # adiabatic continuity: assign each new eigenvector to the
                # previous one with maximal |overlap| (optimal matching)
                Q = np.abs(prev.T @ evec)                     # (old, new)
                row, col = linear_sum_assignment(-Q)
                evec = evec[:, col]
                ev = ev[col]
            sign = np.sign(evec[0, :])
            sign[sign == 0] = 1.0
            evec = evec * sign[None, :]
            self.eps[i] = ev
            self.Phi[i] = evec
            prev = evec

        # ---- per-channel carrier -----------------------------------------
        self.U = {}      # j -> dict(U_plus, U_minus, dU_plus, dU_minus)

    # ------------------------------------------------------------------
    def _q_disp(self, k2):
        """Discrete-dispersion wavenumber solving (theory Eq. 21)

            30 - 32 cos(qh) + 2 cos(2qh) = 12 h^2 k2,

        i.e. the phase for which the 4th-order FD stencil (the reference
        operator) sees exactly the local wavenumber k.  Unique root in
        (0, pi/h) because f'(x) = 4 sin x (8 - cos x) > 0.  Only meaningful
        for k2 > 0 (allowed region); the caller handles the forbidden side
        with the continuous decay constant."""
        tgt = 12.0 * self.h ** 2 * np.asarray(k2, float)
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
        return x / self.h

    # ------------------------------------------------------------------
    def build_channels(self, n_channels=12, mask_eps=None):
        """Construct carriers for the first n_channels adiabatic curves.

        The channel order follows the energy ordering at R_hi (asymptotic
        channel index j), so U^pm_j map onto the physical j channels.
        Returns dict {j: {U+, U-, dU+, dU-}} on the R grid.
        """
        from scipy.special import airy
        idx_hi = self.n_R - 1
        order = np.argsort(self.eps[idx_hi])[:n_channels]
        out = {}
        for jj in order:
            eps_j = self.eps[:, jj]
            # radial equation uses the TOTAL energy E_tot = E_c - D_e
            # (eps_j includes V(R->inf) = -D_e); using E_c here would make
            # every channel "allowed" by ~40 000 cm^-1
            k2 = (self.E_tot - eps_j) / self.pes.BR
            # classically allowed indices
            allow = k2 > 0.0
            if not allow.any():
                continue
            # rightmost sign change gives the asymptotic-side turning point
            signchg = np.where(np.diff(allow.astype(int)) != 0)[0]
            i_c = 0
            for isc in signchg:
                if allow[isc + 1]:      # becoming allowed going outward
                    i_c = isc + 1
            if i_c == 0 and allow[0]:
                i_c = 0
            # discrete-dispersion phase step: q_disp*h in the allowed
            # region (removes the O((qh)^4 * L) phase mismatch with the
            # reference's 4th-order stencil, theory section 9.2), continuous
            # decay constant in the forbidden region
            g = np.where(k2 > 0.0,
                         self._q_disp(np.abs(k2)) * self.h,
                         np.sqrt(np.abs(k2)) * self.h)
            # accumulate |k| from the turning point outwards
            S = np.empty(self.n_R)
            S[i_c] = 0.0
            acc = 0.0
            for i in range(i_c + 1, self.n_R):
                acc += 0.5 * (g[i] + g[i - 1])
                S[i] = acc
                if not allow[i]:
                    # forbidden: keep integrating |k| on the (negative) side
                    pass
            acc = 0.0
            for i in range(i_c - 1, -1, -1):
                acc += 0.5 * (g[i] + g[i + 1])
                S[i] = -acc
            eps_signed = np.sign(self.R - self.R[i_c]) * (1.5 * np.abs(S)) ** (2.0 / 3.0)
            zeta = np.where(np.abs(self.R - self.R[i_c]) < 1e-12, 0.0, eps_signed)
            amp = np.empty(self.n_R)
            deps = np.gradient(eps_j, self.h)
            for i in range(self.n_R):
                if i == i_c:
                    vp = abs(deps[i])
                    amp[i] = (self.pes.BR / max(vp, 1e-30)) ** (1.0 / 6.0)
                else:
                    den = abs(k2[i])
                    amp[i] = (np.abs(zeta[i]) / max(den, 1e-300)) ** 0.25
            # airy(z) -> (Ai, Ai', Bi, Bi')
            Ai, _, Bi, _ = airy(-zeta)
            U_p = amp * (1j * Ai + Bi)       # pure incoming branch
            U_m = amp * (-1j * Ai + Bi)      # pure outgoing branch
            # asymptotic normalization to unit e^{+/- i k R} amplitude:
            # chi = R psi has no 1/sqrt(R) factor, so in the flat region the
            # exact radial factors are e^{+/-i k_j R}.  The projection makes
            # the LS coefficients physically O(1) (c+ ~ 1 for a reflection-
            # free solution) and removes the arbitrary Langer amplitude scale.
            k_asym = np.sqrt(max((self.E_tot - eps_j[-1]) / self.pes.BR, 0.0))
            if k_asym > 0:
                win = self.R > (self.R[-1] - 0.25 * (self.R[-1] - self.R[0]))
                Mp = np.column_stack([np.exp(1j * k_asym * self.R[win]),
                                      np.exp(-1j * k_asym * self.R[win])])
                cp = np.linalg.lstsq(Mp, U_p[win], rcond=None)[0]
                cm = np.linalg.lstsq(Mp, U_m[win], rcond=None)[0]
                if abs(cp[0]) > 1e-12:
                    U_p = U_p / cp[0]
                if abs(cm[1]) > 1e-12:
                    U_m = U_m / cm[1]
            dU_p = self._d1(U_p)
            dU_m = self._d1(U_m)
            out[jj] = dict(U_plus=U_p, U_minus=U_m,
                           dU_plus=dU_p, dU_minus=dU_m,
                           eps=eps_j, i_c=i_c)
        self.U = out
        return out

    def _d1(self, f):
        """First derivative by np.gradient (2nd order); used for diagnostics
        only -- the residual uses the same discrete operator as the reference."""
        return np.gradient(f, self.h, edge_order=2)

    # ------------------------------------------------------------------
    def basis_fields(self, channels):
        """Grid functions of the carrier basis: list of (name, field) with
        field (n_R, n_gamma) = U^pm_j(R) Phi_j(R,gamma)."""
        cols = []
        for jj in channels:
            d = self.U.get(jj)
            if d is None:
                continue
            Phi = self.Phi[:, :, jj]                      # (n_R, n_gamma)
            cols.append((f'U+_{jj}', d['U_plus'][:, None] * Phi))
            cols.append((f'U-_{jj}', d['U_minus'][:, None] * Phi))
        return cols
