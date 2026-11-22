"""Carrier-based physics-informed linear collocation for the 2D scattered wave.

Theory: docs/theory_scattering_2d.md section 7.  The ansatz is

    chi_sc = sum_j [ c+_j U^+_j(R) + c-_j U^-_j(R) ] Phi_j(R,gamma)   (carrier)
             + sum_m d_m w(R) cos(w_m . x + b_m)                       (correction)

where U^pm_j are the adiabatic-channel uniform Langer-Airy functions (they
carry the e^{i k_j(R) R} oscillation analytically) and Phi_j are the adiabatic
angular channel functions.  Because the carrier absorbs the high-frequency
phase, the remaining correction field is smooth (envelope wavenumber ~
O(V'/k) << k), so a few low-frequency features suffice -- restoring the
best-approximation argument of theory section 3.

The residual is formed with the SAME discrete operator as the converged FD
reference (src/conventional/scatter_2d.assemble_2d_scatter) and solved as a
cross-validated ridge least-squares problem; the ridge is chosen on held-out
PDE residual, never on the reference solution.  If the correction features
do not improve the validation residual, the solver falls back to the pure
carrier (honest gate, same protocol as the 1D ScatteringPINN1D).
"""
from __future__ import annotations

import numpy as np

from src.conventional.semiclassical_2d import AdiabaticCarrier2D


class CarrierPINN2D:
    def __init__(self, pes, E_c, R_lo=0.6, R_hi=12.0, n_R=400, n_gamma=80,
                 cap_start=5.0, cap_eta=100.0, cap_width=2.0,
                 n_channels=10, n_corr=0, sigma_R=(8.0, 20.0), sigma_g=6.0,
                 seed=0, val_frac=0.2, wall_window=0.6, wall_weight=1e3):
        from src.conventional.reference import LegendreAngularDVR
        from src.conventional.scatter_2d import assemble_2d_scatter
        self.pes = pes
        self.E_c = float(E_c)
        self.cap = dict(cap_start=cap_start, cap_eta=cap_eta, cap_width=cap_width)
        self.n_corr = int(n_corr)
        self.val_frac = float(val_frac)
        # wall-region Dirichlet window: the scattered wave is exponentially
        # small for R < R_c (the carrier's Bi branch GROWS in the forbidden
        # region, so without this the wall rows would dominate the residual).
        self.wall_window = float(wall_window)
        self.wall_weight = float(wall_weight)

        # discrete operator of the reference (same stencils / BC rows)
        A, b, R, xg, K, E_tot = assemble_2d_scatter(
            pes, E_c, R_lo=R_lo, R_hi=R_hi, n_R=n_R, n_gamma=n_gamma,
            cap_start=cap_start, cap_eta=cap_eta, cap_width=cap_width)
        self.A, self.b, self.R, self.xg, self.K, self.E_tot = A, b, R, xg, K, E_tot
        self.gam = np.arccos(xg)
        self.n_R, self.n_gamma = len(R), len(xg)
        ang = LegendreAngularDVR(self.n_gamma)
        self.carrier = AdiabaticCarrier2D(pes, E_c, R, xg, np.asarray(ang.j2))
        self.channels = self.carrier.build_channels(n_channels=n_channels)
        self.wall_idx = np.where(np.repeat(self.R, self.n_gamma)
                                 <= R_lo + self.wall_window)[0]
        self.names = []
        self.B = self._carrier_columns()
        self.corr_names, self.Bc = self._correction_columns(sigma_R, sigma_g,
                                                           seed) if n_corr else ([], np.zeros((0, 0), complex))

    # ------------------------------------------------------------------
    def _carrier_columns(self):
        cols, names = [], []
        for jj, d in self.channels.items():
            Phi = self.carrier.Phi[:, :, jj]
            for tag, u in (('+', d['U_plus']), ('-', d['U_minus'])):
                cols.append((u[:, None] * Phi).ravel())
                names.append(f'U{tag}_{jj}')
        self.names = names
        return np.column_stack(cols) if cols else np.zeros((self.n_R * self.n_gamma, 0), complex)

    def _correction_columns(self, sigma_R, sigma_g, seed):
        rng = np.random.default_rng(seed)
        n = self.n_R * self.n_gamma
        M = self.n_corr
        cols, names = [], []
        n_per = [M // len(sigma_R)] * len(sigma_R)
        n_per[-1] += M - sum(n_per)
        dR = self.R[-1] - self.R[0]
        a_R, b_Rn = 2.0 / dR, -(self.R[-1] + self.R[0]) / dR
        ag = 1.0 / np.pi
        RR, GG = np.meshgrid(self.R, self.gam, indexing='ij')
        Rh = (a_R * RR + b_Rn).ravel()
        gh = (ag * GG).ravel()
        env = ((RR - self.R[0]) * (self.R[-1] - RR)).ravel()
        for sR, m in zip(sigma_R, n_per):
            wr = rng.normal(0.0, sR, m)[:, None]
            wg = rng.normal(0.0, sigma_g, m)[:, None]
            bb = rng.uniform(0.0, 2.0 * np.pi, m)[:, None]
            cols.append((np.cos(wr * Rh[None] + wg * gh[None] + bb) * env[None]).T)
            names += [f'cor_{i}_{sR}' for i in range(m)]
        return names, np.column_stack(cols)

    # ------------------------------------------------------------------
    def _residual_matrix(self, B):
        L = self.A @ B                              # (n_nodes, n_basis)
        if self.wall_weight > 0 and len(self.wall_idx):
            # rows requiring chi = 0 over the wall window (independent of A's
            # single Dirichlet row at R_lo)
            L = np.vstack([L, self.wall_weight * B[self.wall_idx]])
        return L

    def _rhs(self):
        b = self.b.copy()
        if self.wall_weight > 0 and len(self.wall_idx):
            b = np.concatenate([b, np.zeros(len(self.wall_idx), dtype=complex)])
        return b

    def solve(self, ridge_logspace=np.arange(-14.0, 0.0, 1.0), verbose=True):
        rng = np.random.default_rng(12345)
        out = {}
        # ---- stage 1: pure carrier -----------------------------------
        bb = self._rhs()
        n_tot = bb.shape[0]
        perm = rng.permutation(n_tot)
        n_fit = int(round((1.0 - self.val_frac) * n_tot))
        fit, val = perm[:n_fit], perm[n_fit:]
        Ab = self._residual_matrix(self.B)
        cv = self._cv_ridge(Ab, bb, fit, val, ridge_logspace, verbose)
        c_car, lam_car, r_car = cv['c'], cv['ridge'], cv['val']
        out['carrier'] = dict(ridge=lam_car, val_residual=r_car,
                              n_basis=self.B.shape[1])
        # ---- stage 2: carrier + correction (honest gate) --------------
        acc = {'accepted': False, 'val_residual': r_car, 'ridge': None}
        if self.n_corr:
            Ball = np.hstack([self.B, self.Bc])
            Aall = self._residual_matrix(Ball)
            cv2 = self._cv_ridge(Aall, bb, fit, val, ridge_logspace, verbose)
            if cv2['val'] < r_car * (1.0 - 1e-6):
                acc = dict(accepted=True, val_residual=cv2['val'],
                           ridge=cv2['ridge'])
                self.c_ = cv2['c']
                self.B_used = Ball
            else:
                self.c_ = c_car
                self.B_used = self.B
        else:
            self.c_ = c_car
            self.B_used = self.B
        out['correction'] = acc
        self.val_residual_ = float(acc['val_residual'] if acc['accepted'] else r_car)
        if verbose:
            print(f"  carrier-only: ridge={lam_car:.1e} val_res={r_car:.4e}")
            if self.n_corr:
                print(f"  correction: accepted={acc['accepted']} "
                      f"val_res={acc['val_residual']:.4e}")
        return out

    def _cv_ridge(self, Aop, b, fit, val, ridge_logspace, verbose):
        # column equilibration: carrier columns carry the operator scale
        # (~1e5) while correction columns are O(1), so an un-equilibrated
        # ridge would be dominated by the largest columns and effectively
        # ignore the smooth correction (this was observed as a val residual
        # stuck at ||b||-level before equilibration).
        coln = np.linalg.norm(Aop[fit], axis=0)
        coln[coln == 0] = 1.0
        A = Aop[fit] / coln[None, :]
        bf = b[fit]
        G = A.conj().T @ A
        rhs = A.conj().T @ bf
        evals, evecs = np.linalg.eigh(G)
        base = evecs.conj().T @ rhs
        Av, bv = Aop[val] / coln[None, :], b[val]
        best = (np.inf, None, None)
        for lg in ridge_logspace:
            lam = 10.0 ** lg
            c_s = evecs @ (base / (evals + lam))
            rv = float(np.linalg.norm(Av @ c_s - bv) / np.linalg.norm(bv))
            if verbose:
                print(f'    ridge 1e{lg:+.0f}: val = {rv:.4e}', flush=True)
            if rv < best[0]:
                best = (rv, lam, c_s)
        c = best[2] / coln
        return dict(val=best[0], ridge=best[1], c=c)

    # ------------------------------------------------------------------
    def evaluate(self):
        """Solution on the native (R, gamma) grid: (n_R, n_gamma)."""
        chi = (self.B_used @ self.c_).reshape(self.n_R, self.n_gamma)
        return chi
