"""Physics-informed linear collocation (random-Fourier-feature ELM) for the
2D A+BC scattered-wave BVP.

Theory (docs/theory_scattering_2d.md): the strong-form nonlinear residual
training of ScatteringPINN2D plateaus at a scaled residual ~7e-2, which by
the error-residual conversion delta_chi ~ eps*S/E_c predicts an O(1) solution
error -- exactly the observed 122-156% L2.  The failure is structural, not a
capacity issue.  Restricting the ansatz to the linear span of fixed random
Fourier features turns the problem into a single ridge least-squares solve
whose error is controlled by the best-approximation error of the span, not
by the inverse of the second-derivative operator.

Formulation identical to the FD reference (src/conventional/scatter_2d.py):
    (H + CAP - E_tot) chi_sc = -(V + D_e) e^{+iKR},   chi_sc(R_lo,R_hi) = 0,
with the quadratic CAP absorber.  The ansatz
    chi_sc = w(R) * sum_n c_n cos(w_n . x + b_n),   w(R) = (R-R_lo)(R_hi-R)
satisfies both Dirichlet conditions exactly; L phi_n is computed ANALYTICALLY
(no autograd, no FD on features):

For phi = cos(th), th = wR*Rhat + wg*gammahat + b:
    d2_R phi   = -(2 wR / dR)^2 phi
    j2 phi = -(d2_g + cot(g) d_g) phi = wg^2 phi + wg cot(g) sin(th).

Honesty: the ridge parameter is chosen by cross-validated PDE residual on a
held-out node split; the FD reference is never touched during selection.
"""
import numpy as np




class ScatteringPINN2DELM:
    """Random-Fourier-feature linear collocation for chi_sc (Eq. 1)."""

    def __init__(self, pes, E_c, R_lo=0.6, R_hi=12.0, n_quad=(400, 80),
                 n_feat=2000, sigma_R=(8.0, 20.0), sigma_g=6.0, seed=0,
                 cap_start=None, cap_eta=None, cap_width=2.0, val_frac=0.2):
        from src.conventional.hamiltonians import Water2DProblem
        self.pes = pes
        self.E_c = float(E_c)
        self.E_tot = E_c + float(pes.V_valence(pes.re, pes.re, pes.theta_e)) + pes.De
        self.R_lo, self.R_hi = R_lo, R_hi
        self.K = float(np.sqrt(E_c / pes.BR))
        self.problem = Water2DProblem(pes, R_lo=R_lo, R_hi=R_hi, gate_width=0.3)
        self.cap_start, self.cap_eta, self.cap_width = cap_start, cap_eta, cap_width
        self.val_frac = float(val_frac)
        self._n_gamma = int(n_quad[1])

        rng = np.random.default_rng(seed)
        self.n_feat = int(n_feat)
        # split features across R-frequency scales (envelope vs oscillation)
        n_per = [self.n_feat // len(sigma_R)] * len(sigma_R)
        n_per[-1] += self.n_feat - sum(n_per)
        wR, wg, b = [], [], []
        for sR, n_m in zip(sigma_R, n_per):
            wR.append(rng.normal(0.0, sR, n_m))
            wg.append(rng.normal(0.0, sigma_g, n_m))
            b.append(rng.uniform(0.0, 2.0 * np.pi, n_m))
        self.wR = np.concatenate(wR)[:, None]          # (M, 1)
        self.wg = np.concatenate(wg)[:, None]          # (M, 1)
        self.b = np.concatenate(b)[:, None]            # (M, 1)

        # collocation nodes (Gauss-Legendre; sin(gamma) quadrature weight
        # folded into row weights so the LS functional is the discretized L2
        # norm on the physical measure sin(gamma) dgamma dR)
        pts, w = self.problem.quadrature(n_quad)
        self.R = pts[0]
        self.gam = pts[1]
        self.wrow = np.sqrt(np.asarray(w, float))

        dR = self.R_hi - self.R_lo
        self.a_R, self.b_Rn = 2.0 / dR, -(self.R_hi + self.R_lo) / dR
        self.a_g = 1.0 / np.pi

    # ---- feature internals ---------------------------------------------
    def _theta(self, R, gam):
        """th = wR*Rhat + wg*gammahat + b; returns shape (M, N)."""
        Rh = self.a_R * np.asarray(R, float) + self.b_Rn
        gh = self.a_g * np.asarray(gam, float)
        return self.wR * Rh[None, :] + self.wg * gh[None, :] + self.b

    def _envelope(self, R):
        return (R - self.R_lo) * (self.R_hi - R)

    def _V_cap(self, R, gam):
        V = self.pes.V_jacobi(self.pes.re, R, gam)
        if self.cap_eta:
            s = np.clip((R - self.cap_start) / self.cap_width, 0.0, None)
            V = V - 1j * self.cap_eta * s ** 2
        return V

    def _apply_L(self, R, gam=None):
        """Rows of L phi_n at nodes (R, gam): (N, M) complex.  gam defaults
        to the training node array; self-tests pass their own node pairs.
        The radial second derivative uses the full product rule -- w is
        quadratic, so L[w phi] != w L[phi]:
            d2_R(w phi) = w phi'' + 2 w' phi' + w'' phi,
            w = (R-R_lo)(R_hi-R), w' = R_hi+R_lo-2R, w'' = -2.
        The multiplicative terms (C j2, V+CAP, -E) are local: they pick up a
        single factor of w."""
        if gam is None:
            gam = self.gam
        th = self._theta(R, gam)                         # (M, N)
        phi = np.cos(th)
        w = self._envelope(R)                            # (N,)
        wp = self.R_hi + self.R_lo - 2.0 * R
        wpp = -2.0 * np.ones_like(R)
        dR_phi = -self.wR * self.a_R * np.sin(th)        # (M, N)
        d2R_phi = -(self.wR * self.a_R) ** 2 * phi       # (M, N)
        cot = np.cos(gam) / np.sin(gam)
        # chain rule through gammahat = a_g * gamma (a_g = 1/pi):
        # j2 phi = -(d2g + cot dg) phi = (wg a_g)^2 phi + wg a_g cot sin(th)
        wg_eff = self.wg * self.a_g
        j2_phi = wg_eff ** 2 * phi + wg_eff * cot[None, :] * np.sin(th)
        C = self.pes.BR / R ** 2 + self.pes.Br / self.pes.re ** 2
        radial = -self.pes.BR * (w[None, :] * d2R_phi
                                 + 2.0 * wp[None, :] * dR_phi
                                 + wpp[None, :] * phi)   # (M, N)
        local = (C * j2_phi + self._V_cap(R, gam)[None, :] * phi
                 - self.E_tot * phi) * w[None, :]        # (M, N)
        return (radial + local).T                        # (N, M)

    def _source(self, R):
        """-(H - E) chi_in = -(V + D_e) e^{+iKR} at the nodes (physical V)."""
        V = self.pes.V_jacobi(self.pes.re, R, self.gam)
        return -(V + self.pes.De) * np.exp(1j * self.K * R)
    # NOTE: _source is only called with the full training node array
    # (self.R, self.gam); _V_cap is parameterized for reuse in self-tests.

    # ---- solve ----------------------------------------------------------
    def solve(self, ridge_logspace=np.arange(-14.0, 0.0, 1.0), verbose=True):
        """Cross-validated ridge least-squares; selection never sees the
        reference solution.  M <= 5000: exact hermitian-eigendecomposition
        path (all lambdas from one eigh).  Larger M: scipy.sparse.linalg
        LSQR per lambda (memory-safe).  Returns dict with the lambda sweep."""
        import scipy.sparse.linalg as spla
        A = self._apply_L(self.R) * self.wrow[:, None]      # (N, M)
        bvec = self._source(self.R) * self.wrow             # (N,)
        n = len(bvec)
        rng = np.random.default_rng(12345)
        perm = rng.permutation(n)
        n_fit = int(round((1.0 - self.val_frac) * n))
        fit_idx, val_idx = perm[:n_fit], perm[n_fit:]
        Av, bv = A[val_idx], bvec[val_idx]
        out = {'ridge': [], 'val_residual': []}
        best = (np.inf, None, None)

        if self.n_feat <= 5000:
            Af, bf = A[fit_idx], bvec[fit_idx]
            G = Af.conj().T @ Af                            # (M, M) hermitian
            rhs = Af.conj().T @ bf
            evals, evecs = np.linalg.eigh(G)
            base = evecs.conj().T @ rhs                     # (M,)
            for lg in ridge_logspace:
                lam = 10.0 ** lg
                c = evecs @ (base / (evals + lam))
                rv = float(np.linalg.norm(Av @ c - bv) / np.linalg.norm(bv))
                out['ridge'].append(float(lam))
                out['val_residual'].append(rv)
                if verbose:
                    print(f'  ridge 1e{lg:+.0f}: val PDE residual = {rv:.4e}',
                          flush=True)
                if rv < best[0]:
                    best = (rv, lg, c)
        else:
            # column scaling for conditioning; ridge rows appended to the LSQR
            coln = np.linalg.norm(A[fit_idx], axis=0)
            coln[coln == 0] = 1.0
            As = A / coln[None, :]
            for lg in ridge_logspace:
                lam = 10.0 ** lg
                A_aug = np.vstack([As[fit_idx],
                                   lam * np.eye(self.n_feat, dtype=complex)])
                b_aug = np.concatenate([bvec[fit_idx],
                                        np.zeros(self.n_feat, dtype=complex)])
                c = spla.lsqr(A_aug, b_aug, iter_lim=600, atol=0, btol=0)[0]
                # As @ c_s == A @ (c_s / coln): residual in original units
                rv = float(np.linalg.norm(As[val_idx] @ c - bv)
                           / np.linalg.norm(bv))
                out['ridge'].append(float(lam))
                out['val_residual'].append(rv)
                if verbose:
                    print(f'  ridge 1e{lg:+.0f} (LSQR): val PDE residual = '
                          f'{rv:.4e}', flush=True)
                if rv < best[0]:
                    best = (rv, lg, c)
        self.ridge_ = 10.0 ** best[1]
        self.c_ = None if best[2] is None else (
            best[2] / coln if self.n_feat > 5000 else best[2])
        self.val_residual_ = best[0]
        out['best_ridge'] = self.ridge_
        out['best_val_residual'] = best[0]
        if verbose:
            print(f'  selected ridge = {self.ridge_:.2e}, '
                  f'val residual = {best[0]:.4e}', flush=True)
        return out

    # ---- evaluation & asymptotic extraction -----------------------------
    def evaluate(self, R_grid, gamma_grid):
        RR, GG = np.meshgrid(R_grid, gamma_grid, indexing='ij')
        th = self._theta(RR.ravel(), GG.ravel())            # (M, Np)
        phi = np.cos(th)
        env = self._envelope(RR.ravel())
        chi = (phi.T @ self.c_) * env
        return chi.reshape(RR.shape)

    def extract_hankel(self, R_match=(4.2, 5.0), n_R=40, kappa_cut=40.0):
        """Least-squares projection of the solution onto the EXACT multichannel
        asymptotic basis (theory Eq. 5, zero adjustable parameters): open
        channels -> outgoing Hankel H^(2)_nu(k_j R); weakly closed channels
        (kappa*R_match < kappa_cut) -> the decaying modified Bessel K_nu,
        computed via the scaled form kve*exp(-kap R) to avoid overflow
        (K_nu(nu~6000) overflows float64 otherwise); deeply closed channels
        are excluded -- their amplitude at the matching region is below
        e^{-kappa_cut}, standard close-coupling practice.  Angular structure
        = DVR eigenvectors of j^2.  Returns coeffs, kept channel indices,
        wavenumbers and the open-channel flux ratio sum |c_j|^2 k_j / k_0
        (a physics self-check: <= 1 up to small CAP absorption)."""
        from scipy.special import hankel2, kve
        from src.conventional.reference import LegendreAngularDVR
        ng = self._n_gamma
        ang = LegendreAngularDVR(ng)
        evals, evecs = np.linalg.eigh(np.asarray(ang.j2))
        B_rot = self.pes.Br / self.pes.re ** 2
        Rm = np.linspace(R_match[0], R_match[1], n_R)
        gam = np.arccos(ang.x)          # gamma nodes matching the DVR x grid
        chi = self.evaluate(Rm, gam)                        # (n_R, ng)
        F_cols, kept, k_list, is_open = [], [], [], []
        for ch in range(ng):
            eps_ch = max(evals[ch], 0.0)
            E_ch = self.E_c - B_rot * eps_ch
            if E_ch > 0:
                k = np.sqrt(E_ch / self.pes.BR)
                F_cols.append(hankel2(eps_ch, k * Rm))
                kept.append(ch)
                k_list.append(k)
                is_open.append(True)
            else:
                kap = np.sqrt(-E_ch / self.pes.BR)
                if kap * R_match[0] < kappa_cut:
                    # K_nu(kap R) = kve(nu, kap R) * exp(-kap R) (scaled form)
                    F_cols.append(kve(eps_ch, kap * Rm) * np.exp(-kap * Rm))
                    kept.append(ch)
                    k_list.append(1j * kap)
                    is_open.append(False)
        F = np.array(F_cols)                                # (n_kept, n_R)
        evecs_k = evecs[:, kept]
        # A[iR*ng + ig, ch] = F_ch(R_iR) * Phi_ch(gamma_ig)
        A = np.einsum('rc,gc->rgc', F.T, evecs_k).reshape(len(Rm) * ng,
                                                          len(kept))
        c, *_ = np.linalg.lstsq(A, chi.ravel(), rcond=None)
        is_open = np.array(is_open)
        k_open = np.array([k.real for k, o in zip(k_list, is_open) if o])
        flux = float(np.sum(np.abs(c)[is_open] ** 2 * k_open) / self.K)
        return dict(coeffs=c, kept_channels=np.array(kept),
                    eps=evals[kept], k=np.array(k_list), flux_ratio=flux,
                    n_open=int(is_open.sum()), R_match=R_match)
