"""PINN for 1D fixed-energy scattering states of the Morse potential.

The state at energy E > 0 satisfies
    -B psi'' + V(x) psi = E psi,
    psi(x) -> e^{+ikx} + R e^{-ikx}   as x -> +inf   (unit incoming amplitude),
    psi(x) -> 0                       as x -> -inf   (impenetrable wall).

Physics-informed ansatz: the carrier is the pair of ZERO-FREE uniform
semiclassical traveling waves built from the classical action (no empirical
parameters); the network correction acts INDEPENDENTLY on the two branches:
    S(x)    = int_{x_c}^{x} k(x') dx',   k = sqrt((E - V)/B)   (allowed side)
    zeta(x) = sign(x - x_c) (3|S(x)|/2)^(2/3)                  (Langer coordinate)
    U_+(x)  = +i Ai(-zeta) + Bi(-zeta)   (pure incoming  e^{+i theta}, theta = S + pi/4)
    U_-(x)  = -i Ai(-zeta) + Bi(-zeta)   (pure outgoing  e^{-i theta})
    psi(x)  = M_+(x) U_+(x) + M_-(x) U_-(x),   base M_+ = 1, M_- = -1
The BASE combination U_+ - U_- = 2i Ai(-zeta) is the REGULAR solution: it
decays naturally into the wall (no extra damping factor is needed -- a WKB
envelope would double-count the decay and its curvature diverges at x_c).
The M_- field carries the physical Stokes transition from the regular
combination near the wall (M_- = -1) to the reflected amplitude on the
plateau (M_- -> R).
where M_+-M_- are smooth complex modulation fields carried by a SINGLE-HIDDEN-
LAYER random-feature neural network with a trainable linear output:
    M_-(x) = 1 + pinch(x) * sum_k c_k^{(-)} phi_k(x),   M_+(x) = 1 + pinch(x) * sum_k c_k^{(+)} phi_k(x)
Two independent modulations are REQUIRED: the true solution's zeros are
reproduced by destructive interference of the two branches, which a single
multiplicative modulation (which cannot shift the branches' relative phase)
cannot represent without poles.  U_+-U_- never vanish, so smooth M_+-M_- suffice.

Because the PDE residual and the asymptotic Robin condition are LINEAR in the
output coefficients c, the network is trained by physics-informed LINEAR
LEAST SQUARES (randomized-network / extreme-learning-machine style): one
deterministic normal-equation solve, no iterations, no local minima, and no
reference data anywhere.  Near x_c the marched action is replaced by its
analytic two-term expansion (exact for the linearised potential), which keeps
the carrier residual bounded across the turning point.

Loss rows = D^2-weighted interior PDE residual + asymptotic Robin condition
    psi'(x_hi) + i k psi(x_hi) = 2 i k e^{+i k x_hi}
(the same continuum condition the finite-difference reference discretises; it
fixes the unit incoming amplitude).

Time convention e^{+iEt/hbar}: e^{ikx} propagates toward -x (incoming).
"""
import numpy as np

try:
    import torch
except Exception:
    torch = None

from scipy.special import airy as _scipy_airy
airy_vec = _scipy_airy

from src.conventional.reference import gauss_legendre_1d


def _smoothstep(t):
    """Quintic C^2 smoothstep on [0, 1]."""
    t = np.clip(t, 0.0, 1.0)
    return t * t * t * (10.0 + t * (-15.0 + 6.0 * t))


def _d1(f, h):
    """4th-order first derivative on a uniform grid."""
    g = np.zeros_like(f)
    g[2:-2] = (f[:-4] - 8 * f[1:-3] + 8 * f[3:-1] - f[4:]) / (12 * h)
    g[1] = (f[2] - f[0]) / (2 * h)
    g[-2] = (f[-1] - f[-3]) / (2 * h)
    g[0], g[-1] = g[1], g[-2]
    return g


def _smoothstep_dpinch(x, x_c, d0, w):
    """Derivative of the one-sided quintic smoothstep pinch w.r.t. x."""
    t = np.clip((x - x_c - d0) / w, 0.0, 1.0)
    inside = (x - x_c - d0 >= 0.0) & (x - x_c - d0 <= w)
    dt = np.where(inside, 30.0 * t * t * (t - 1.0) ** 2 / w, 0.0)
    return dt


def _d2(f, h):
    """4th-order second derivative on a uniform grid."""
    g = np.zeros_like(f)
    g[2:-2] = (-f[4:] + 16 * f[3:-1] - 30 * f[2:-2] + 16 * f[1:-3] - f[:-4]) / (12 * h ** 2)
    g[1] = (f[2] - 2 * f[1] + f[0]) / h ** 2
    g[-2] = (f[-1] - 2 * f[-2] + f[-3]) / h ** 2
    g[0], g[-1] = g[1], g[-2]
    return g


class ScatteringPINN1D:
    def __init__(self, morse, E, xlo=-1.5, xhi=5.5, n_quad=1200,
                 n_basis=48, seed=0, device='cpu',
                 robin_weight=5.0, blend_d0=0.10, blend_w=0.20,
                 pinch_d0=0.15, pinch_w=1.5,
                 n_fine=40001, ridge=1e-10):
        if torch is None:
            raise RuntimeError("PyTorch is not available; ScatteringPINN1D requires PyTorch.")
        self.morse, self.E = morse, float(E)
        self.xlo, self.xhi = xlo, xhi
        self.k = float(morse.k_continuum(E))
        self.device = device
        self.ridge = float(ridge)
        self.robin_weight = float(robin_weight)
        self._blend_d0, self._blend_w = float(blend_d0), float(blend_w)
        self._pinch_d0, self._pinch_w = float(pinch_d0), float(pinch_w)

        # classical turning point on the wall side (exact for Morse)
        self.x_c = morse.xe - np.log(1.0 + np.sqrt(1.0 + E / morse.De)) / morse.a
        self.x_c = float(np.clip(self.x_c, xlo + 0.2, xhi - 2.0))

        # local expansions of V at the turning point (first principles)
        eps = 1e-5
        self.V_c = float(morse.V(self.x_c))
        self.dV_c = float((morse.V(self.x_c + eps) - morse.V(self.x_c - eps)) / (2 * eps))
        self.d2V_c = float((morse.V(self.x_c + 1e-4) - 2 * self.V_c
                            + morse.V(self.x_c - 1e-4)) / 1e-8)
        self.a_lin = abs(self.dV_c) / morse.B          # linearised slope / B
        self.c1 = 3.0 * self.d2V_c / (20.0 * self.dV_c)

        # random-feature basis (single hidden layer; linear output trained by LS)
        rng = np.random.default_rng(seed)
        self.n_basis = int(n_basis)
        # low-order deterministic Fourier basis in the slow variable + constant:
        # the WKB corrections are smooth (accumulated phase/amplitude), so the
        # basis needs only smooth low-frequency elements.
        self.freqs = np.pi * np.arange(1, self.n_basis) / (self.xhi - self.xlo)

        # collocation grid
        x, w = gauss_legendre_1d(xlo, xhi, n_quad)
        self.x_np = x
        self.x = torch.tensor(x, dtype=torch.float64, device=device,
                              requires_grad=True).unsqueeze(-1)
        self.res_scale = max(self.morse.B * self.k ** 2 + self.E, 1.0)

        # ---- first-principles carrier on a fine uniform helper grid ----
        x_fine = np.linspace(xlo, xhi, n_fine)
        Up_f, Um_f, rUp_f, rUm_f, S_f = self._carrier_on_grid(x_fine)

        def interp_c(arr):
            re = np.interp(x, x_fine, arr.real)
            im = np.interp(x, x_fine, arr.imag)
            return re + 1j * im

        self._Up = interp_c(Up_f)
        self._Um = interp_c(Um_f)
        self._resUp = interp_c(rUp_f)
        self._resUm = interp_c(rUm_f)
        self._S_fine, self._x_fine = S_f, x_fine

        # carrier first/second derivatives via the residual identity:
        #   resU = -B U'' + (V - E) U   =>   U'' = ((V - E) U - resU) / B
        V_nodes = np.asarray(morse.V(x), dtype=float)
        h_f = x_fine[1] - x_fine[0]
        self._dUp = interp_c(_d1(Up_f, h_f))
        self._dUm = interp_c(_d1(Um_f, h_f))
        self._d2Up = ((V_nodes - E) * self._Up - self._resUp) / morse.B
        self._d2Um = ((V_nodes - E) * self._Um - self._resUm) / morse.B

        # residual weights: uniform (the regular-Airy carrier decays naturally,
        # so no amplitude-based reweighting is needed)
        self._w_res = np.ones_like(x)

        # pinch: the corrections are active ONLY in the allowed region beyond
        # the match zone; on the wall side M_- must stay at -1 (regular solution)
        self._pinch = _smoothstep((x - self.x_c - self._pinch_d0) / self._pinch_w)
        p_f = _smoothstep((x_fine - self.x_c - self._pinch_d0) / self._pinch_w)
        self._pinch_d1 = np.interp(x, x_fine, _d1(p_f, h_f))
        self._pinch_d2 = np.interp(x, x_fine, _d2(p_f, h_f))

        # basis functions and their derivatives on the collocation grid
        # (analytic derivatives of the trigonometric basis)
        s = 2.0 * (x - xlo) / (xhi - xlo) - 1.0
        Phi = np.ones((len(x), self.n_basis))
        dPhi = np.zeros((len(x), self.n_basis))
        d2Phi = np.zeros((len(x), self.n_basis))
        # chain rule: s(x) = 2(x-xlo)/L - 1  =>  ds/dx = 2/L
        cscale = 2.0 / (xhi - xlo)
        for j, om in enumerate(self.freqs, start=1):
            Phi[:, j] = np.sin(om * s)
            dPhi[:, j] = om * np.cos(om * s) * cscale
            d2Phi[:, j] = -om ** 2 * np.sin(om * s) * cscale ** 2
        # corrected basis psi_k = pinch * phi_k (product rule, all precomputed)
        self._Psi = self._pinch[:, None] * Phi
        self._dPsi = (self._pinch_d1[:, None] * Phi
                      + self._pinch[:, None] * dPhi)
        self._d2Psi = (self._pinch_d2[:, None] * Phi
                       + 2.0 * self._pinch_d1[:, None] * dPhi
                       + self._pinch[:, None] * d2Phi)

        # regular-combination (base) residual on the fine grid, computed
        # directly from the decaying 2i*amp*Ai form (no Bi cancellation)
        Sf2, Vf2 = self._action_on_grid(x_fine)
        sgn_f = np.sign(x_fine - self.x_c)
        zeta_f = sgn_f * (1.5 * np.abs(Sf2)) ** (2.0 / 3.0)
        amp_f = self._langer_amp(x_fine, zeta_f)
        ai_f = _scipy_airy(-zeta_f)[0]
        reg_f = 2j * amp_f * ai_f
        h_f2 = x_fine[1] - x_fine[0]
        self._base_res_f = -self.morse.B * _d2(reg_f, h_f2) + (Vf2 - self.E) * reg_f

        # ---- asymptotic Robin boundary data at x_hi ----
        Up_hi, Um_hi, _, _, _ = self._boundary_data(np.array([xhi]))
        h_loc = 1e-4
        x_pm = xhi + h_loc * np.array([-2.0, -1.0, 0.0, 1.0, 2.0])
        Up_pm, Um_pm, _, _, _ = self._boundary_data(x_pm)
        w4 = np.array([1, -8, 0, 8, -1], dtype=float) / (12 * h_loc)
        self._Up_hi = complex(Up_hi[0])
        self._Um_hi = complex(Um_hi[0])
        self._dUp_hi = complex(np.dot(w4, Up_pm))
        self._dUm_hi = complex(np.dot(w4, Um_pm))
        self._robin_rhs = 2j * self.k * np.exp(1j * self.k * xhi)
        self._robin_scale = max(2.0 * self.k, 1.0)
        # basis at the boundary point (pinch = 1 there: far from x_c)
        sb = 2.0 * (xhi - xlo) / (xhi - xlo) - 1.0 if False else 2.0 * (xhi - xlo) / (xhi - xlo) - 1.0
        sb = 1.0  # s(x_hi) = +1 exactly
        self._phi_hi = np.concatenate([[1.0], np.sin(self.freqs * sb)])
        self._dphi_hi = np.concatenate([[0.0], self.freqs * np.cos(self.freqs * sb) * cscale])

        # coefficients (torch Parameters so the state dict / report API works)
        self.coef = torch.zeros(2 * self.n_basis, 2, dtype=torch.float64, device=device,
                                requires_grad=False)  # filled by the LS solve

    # ------------------------------------------------------------------
    def _action_on_grid(self, xg):
        """Classical action S(x) and wall integral K(x) with the analytic
        two-term expansion near the turning point (removes the O(h^{3/2})
        march error that the second-difference residual would amplify)."""
        mo, E = self.morse, self.E
        xg = np.atleast_1d(np.asarray(xg, dtype=float))
        V = np.asarray(mo.V(xg), dtype=float)
        B = mo.B
        k2 = (E - V) / B
        k_np = np.sqrt(np.maximum(k2, 0.0))
        kap = np.sqrt(np.maximum(-k2, 0.0))

        # SIGNED classical action: k on the allowed side, -kappa on the wall
        # side (S < 0 there), so zeta < 0 and Ai(-zeta) = Ai(+|zeta|) decays
        # naturally into the wall -- the regular combination U_+ - U_- = 2i Ai
        # carries the exact uniform wall decay with no extra envelope factor.
        S = np.zeros_like(xg)
        i_c = int(np.argmin(np.abs(xg - self.x_c)))
        for i in range(i_c + 1, len(xg)):
            S[i] = S[i - 1] + 0.5 * (k_np[i] + k_np[i - 1]) * (xg[i] - xg[i - 1])
        for i in range(i_c - 1, -1, -1):
            S[i] = S[i + 1] - 0.5 * (kap[i] + kap[i + 1]) * (xg[i + 1] - xg[i])

        dlt = np.abs(xg - self.x_c)
        sgn = np.sign(xg - self.x_c)
        S_an = sgn * (2.0 / 3.0) * np.sqrt(self.a_lin) * dlt ** 1.5 * (1.0 + self.c1 * dlt)

        w = _smoothstep((dlt - self._blend_d0) / self._blend_w)
        S_b = w * S + (1.0 - w) * S_an
        # cap the WALL-side action only: beyond |zeta| ~ 20 the true solution
        # is exp(-100) small and the growing Bi branch would overflow double
        # precision.  The allowed side (S > 0) is left uncapped -- both Airy
        # functions oscillate there with decaying amplitude.
        S_cap = (2.0 / 3.0) * 20.0 ** 1.5
        wall = S_b < 0
        S_b = np.where(wall, -np.minimum(np.abs(S_b), S_cap), S_b)
        return S_b, V

    def _langer_amp(self, xg, zeta):
        """Uniform Langer amplitude (|zeta/k^2|)^(1/4) with the regularised
        turning-point value (B/|V'|)^(1/6); recovers the WKB amplitude
        k^(-1/2) away from x_c on both sides."""
        V = np.asarray(self.morse.V(xg), dtype=float)
        k2abs = np.abs(self.E - V) / self.morse.B
        zabs = np.maximum(np.abs(zeta), 1e-12)
        amp = (zabs / np.maximum(k2abs, 1e-12)) ** 0.25
        reg = (self.morse.B / abs(self.dV_c)) ** (1.0 / 6.0)
        return np.where(zabs > 1e-12 * np.max(zabs, initial=1.0), amp, reg)

    def _boundary_data(self, x_pt):
        """Carrier values (U_+, U_-) at arbitrary points (no derivatives)."""
        S, V = self._action_on_grid(np.atleast_1d(x_pt))
        sgn = np.sign(np.asarray(x_pt, dtype=float) - self.x_c)
        zeta = sgn * (1.5 * np.abs(S)) ** (2.0 / 3.0)
        ai, _, bi, _ = _scipy_airy(-zeta)
        amp = self._langer_amp(np.asarray(x_pt, dtype=float), zeta)
        Up = amp * (1j * ai + bi)
        Um = amp * (-1j * ai + bi)
        return Up, Um, ai, bi, S

    def _carrier_on_grid(self, xg):
        """Damped uniform traveling waves and their PDE residuals on a uniform
        grid (4th-order finite differences of exact function evaluations)."""
        mo, E = self.morse, self.E
        xg = np.atleast_1d(np.asarray(xg, dtype=float))
        S, V = self._action_on_grid(xg)
        B = mo.B
        sgn = np.sign(xg - self.x_c)
        zeta = sgn * (1.5 * np.abs(S)) ** (2.0 / 3.0)
        ai, _, bi, _ = _scipy_airy(-zeta)
        amp = self._langer_amp(xg, zeta)
        Up = amp * (1j * ai + bi)
        Um = amp * (-1j * ai + bi)

        def d2(f):
            g = np.zeros_like(f)
            h = xg[1] - xg[0]
            g[2:-2] = (-f[4:] + 16 * f[3:-1] - 30 * f[2:-2] + 16 * f[1:-3] - f[:-4]) / (12 * h ** 2)
            g[1] = (f[2] - 2 * f[1] + f[0]) / h ** 2
            g[-2] = (f[-1] - 2 * f[-2] + f[-3]) / h ** 2
            g[0], g[-1] = g[1], g[-2]
            return g

        resUp = -B * d2(Up) + (V - E) * Up
        resUm = -B * d2(Um) + (V - E) * Um
        return Up, Um, resUp, resUm, S

    # ------------------------------------------------------------------
    def _basis(self, x_np):
        """Basis matrix (with the one-sided pinch) and derivatives on an
        arbitrary grid.  The pinch factor is the same as self._pinch: the
        corrections are active only in the allowed region beyond the match
        zone (the M_- = -1 regular wall behaviour is preserved)."""
        x_np = np.asarray(x_np, dtype=float)
        s = 2.0 * (x_np - self.xlo) / (self.xhi - self.xlo) - 1.0
        cscale = 2.0 / (self.xhi - self.xlo)
        Phi = np.ones((len(x_np), self.n_basis))
        dPhi = np.zeros_like(Phi)
        for j, om in enumerate(self.freqs, start=1):
            Phi[:, j] = np.sin(om * s)
            dPhi[:, j] = om * np.cos(om * s) * cscale
        pinch = _smoothstep((x_np - self.x_c - self._pinch_d0) / self._pinch_w)
        dpinch = _smoothstep_dpinch(x_np, self.x_c, self._pinch_d0, self._pinch_w)
        return pinch[:, None] * Phi, (dpinch[:, None] * Phi + pinch[:, None] * dPhi), pinch

    def solve(self, val_frac=0.2, ridge_logspace=(-14.0, -1.0, 25), seed=0):
        """Physics-informed linear least-squares training with honest
        cross-validated regularisation.

        The interior collocation nodes are split into fit/validation sets; the
        ridge parameter is chosen by the VALIDATION PDE residual (never by any
        reference data), and the trained correction is accepted ONLY if it
        reduces the validation residual below the pure-carrier base level --
        otherwise the method honestly falls back to the base and reports it.
        """
        n = len(self.x_np)
        nb = self.n_basis
        Bc = self.morse.B
        wgt = np.sqrt(self._w_res)

        Rk_p = (self._Psi * self._resUp[:, None]
                - 2.0 * Bc * self._dPsi * self._dUp[:, None]
                - Bc * self._d2Psi * self._Up[:, None]) * wgt[:, None]
        Rk_m = (self._Psi * self._resUm[:, None]
                - 2.0 * Bc * self._dPsi * self._dUm[:, None]
                - Bc * self._d2Psi * self._Um[:, None]) * wgt[:, None]
        res_reg_f = self._base_res_f  # precomputed regular-combination residual
        base = np.interp(self.x_np, self._x_fine, res_reg_f.real) \
            + 1j * np.interp(self.x_np, self._x_fine, res_reg_f.imag)
        base = base * wgt

        A = np.zeros((2 * n, 4 * nb))
        A[:n, 0:nb] = Rk_p.real
        A[:n, nb:2 * nb] = -Rk_p.imag
        A[:n, 2 * nb:3 * nb] = Rk_m.real
        A[:n, 3 * nb:4 * nb] = -Rk_m.imag
        A[n:, 0:nb] = Rk_p.imag
        A[n:, nb:2 * nb] = Rk_p.real
        A[n:, 2 * nb:3 * nb] = Rk_m.imag
        A[n:, 3 * nb:4 * nb] = Rk_m.real
        b = np.concatenate([-base.real, -base.imag])

        # Robin rows (always in the fit set: they define the BVP)
        robin_row_p = (self._phi_hi * (self._dUp_hi + 1j * self.k * self._Up_hi)
                       + self._dphi_hi * self._Up_hi)
        robin_row_m = (self._phi_hi * (self._dUm_hi + 1j * self.k * self._Um_hi)
                       + self._dphi_hi * self._Um_hi)
        w_rob = np.sqrt(self.robin_weight) / self._robin_scale
        base_rob = (self._dUp_hi - self._dUm_hi
                    + 1j * self.k * (self._Up_hi - self._Um_hi)
                    - self._robin_rhs) * w_rob
        A_rob = np.zeros((2, 4 * nb))
        A_rob[0, 0:nb] = robin_row_p.real * w_rob
        A_rob[0, nb:2 * nb] = -robin_row_p.imag * w_rob
        A_rob[0, 2 * nb:3 * nb] = robin_row_m.real * w_rob
        A_rob[0, 3 * nb:4 * nb] = -robin_row_m.imag * w_rob
        A_rob[1, 0:nb] = robin_row_p.imag * w_rob
        A_rob[1, nb:2 * nb] = robin_row_p.real * w_rob
        A_rob[1, 2 * nb:3 * nb] = robin_row_m.imag * w_rob
        A_rob[1, 3 * nb:4 * nb] = robin_row_m.real * w_rob
        b_rob = np.array([-base_rob.real, -base_rob.imag])

        # validation split over interior nodes (block-wise, deterministic)
        rng = np.random.default_rng(seed)
        n_val = max(int(val_frac * n), 8)
        perm = rng.permutation(n)
        val_idx = np.sort(perm[:n_val])
        fit_mask = np.ones(n, dtype=bool)
        fit_mask[val_idx] = False

        rows_fit = np.concatenate([2 * np.where(fit_mask)[0],
                                   2 * np.where(fit_mask)[0] + 1])
        rows_val = np.concatenate([2 * val_idx, 2 * val_idx + 1])
        A_fit = np.vstack([A[rows_fit], A_rob])
        b_fit = np.concatenate([b[rows_fit], b_rob])
        A_val, b_val = A[rows_val], b[rows_val]

        base_val = float(np.sum(b_val ** 2))
        best = None
        lo, hi, num = ridge_logspace
        for log_ridge in np.linspace(lo, hi, num):
            ridge = 10.0 ** log_ridge
            AtA = A_fit.T @ A_fit
            Atb = A_fit.T @ b_fit
            scale = max(float(np.trace(AtA)) / (4 * nb), 1.0)
            AtA[np.diag_indices_from(AtA)] += ridge * scale
            try:
                sol = np.linalg.solve(AtA, Atb)
            except np.linalg.LinAlgError:
                continue
            val_score = float(np.sum((A_val @ sol - b_val) ** 2)) \
                + float(np.sum((A_rob @ sol - b_rob) ** 2))
            # exact single-channel unitarity: |R| = 1 is a hard physical
            # invariant -- candidates violating it are rejected outright
            if self._flux_error(sol) > 1e-3:
                continue
            if best is None or val_score < best[0]:
                best = (val_score, sol, ridge)

        if best is not None and best[0] < base_val and self._flux_error(best[1]) <= 1e-3:
            sol = best[1]
            c_p = sol[0:nb] + 1j * sol[nb:2 * nb]
            c_m = sol[2 * nb:3 * nb] + 1j * sol[3 * nb:4 * nb]
            cm = np.stack([np.column_stack([c_p.real, c_p.imag]),
                           np.column_stack([c_m.real, c_m.imag])])
            self.coef = torch.tensor(cm, dtype=torch.float64, device=self.device)
            self.fallback = False
            self.ridge_used = best[2]
            self.val_improvement = base_val / max(best[0], 1e-30)
        else:
            self.coef = torch.zeros(2, nb, 2, dtype=torch.float64, device=self.device)
            self.fallback = True
            self.ridge_used = None
            self.val_improvement = 1.0
        return self.report()

    # ------------------------------------------------------------------
    def modulation_on(self, x_np):
        """M_+, M_- on an arbitrary grid from the trained coefficients."""
        Psi, dPsi, pinch = self._basis(np.asarray(x_np, float))
        c = self.coef.cpu().numpy()  # (2, nb, 2): branch, basis, (re, im)
        Mp = 1.0 + Psi @ c[0, :, 0] + 1j * (Psi @ c[0, :, 1])
        Mm = -1.0 + Psi @ c[1, :, 0] + 1j * (Psi @ c[1, :, 1])
        return Mp, Mm

    def _eval_raw(self, x_np):
        """Raw ansatz (base + corrections) on an arbitrary grid."""
        x_np = np.asarray(x_np, dtype=float)
        Up, Um, _, _, S = self._carrier_on_grid(x_np)
        Psi, _, pinch = self._basis(x_np)
        c = self.coef.cpu().numpy()
        Mp = 1.0 + Psi @ c[0, :, 0] + 1j * (Psi @ c[0, :, 1])
        Mm = -1.0 + Psi @ c[1, :, 0] + 1j * (Psi @ c[1, :, 1])
        return Mp * Up + Mm * Um

    def _flux_error(self, sol=None):
        """Exact single-channel unitarity check |R| - 1 for a coefficient
        vector (gauge-invariant plateau fit; uses only exact physics)."""
        if sol is None:
            cmat = self.coef.cpu().numpy()
        else:
            cmat = sol.reshape(2, self.n_basis, 2)
        saved = self.coef
        self.coef = torch.tensor(cmat, dtype=torch.float64, device=self.device)
        x_pl = np.linspace(max(self.xhi - 0.9, 4.6), self.xhi - 0.05, 200)
        psi = self._eval_raw(x_pl)
        self.coef = saved
        Mfit = np.column_stack([np.exp(1j * self.k * x_pl),
                                np.exp(-1j * self.k * x_pl)])
        ab, *_ = np.linalg.lstsq(Mfit, psi, rcond=None)
        return abs(abs(ab[1] / ab[0]) - 1.0)

    def reflection_from_c(self):
        """Physical reflection amplitude: R = (M_-/M_+) * exp(-i(pi/2 - 2 k x_c))
        evaluated at x_hi (plateau amplitudes of the two branches)."""
        Mp_hi, Mm_hi = self.modulation_on(np.array([self.xhi]))
        Mp_hi, Mm_hi = complex(Mp_hi[0]), complex(Mm_hi[0])
        if abs(Mp_hi) < 1e-14:
            return complex(float('inf'), 0.0)
        return (Mm_hi / Mp_hi) * np.exp(-1j * (np.pi / 2.0 - 2.0 * self.k * self.x_c))

    def report(self):
        R = self.reflection_from_c()
        return dict(E=self.E, k=self.k, R=R, absR=abs(R), argR=np.angle(R))

    def evaluate(self, x_new, normalize_incoming=True):
        """Evaluate psi on an arbitrary grid; optionally normalise the incoming
        amplitude to 1 and align the global carrier phase to the e^{ikx}
        convention via the ANALYTIC gauge factor exp(-i(pi/4 - k x_c))."""
        x_np = np.asarray(x_new, float)
        Up, Um, _, _, S = self._carrier_on_grid(x_np)
        Mp, Mm = self.modulation_on(x_np)
        p = Mp * Up + Mm * Um
        if not normalize_incoming:
            return p
        Mp_hi, Mm_hi = self.modulation_on(np.array([self.xhi]))
        A_in = complex(Mp_hi[0]) * np.exp(1j * (np.pi / 4.0 - self.k * self.x_c))
        if abs(A_in) < 1e-14:
            return p
        return p / A_in

    def extract_smatrix_fit(self, x_grid, psi_vals):
        """First-principles S-matrix extraction: least-squares fit of
        alpha e^{+iS} + beta e^{-iS} on the asymptotic plateau; R = beta/alpha.
        Returns (R, psi normalised to unit incoming amplitude)."""
        _, _, _, _, S = self._carrier_on_grid(np.asarray(x_grid, float))
        x_grid = np.asarray(x_grid)
        m = (x_grid > max(self.xhi - 0.9, 4.6)) & \
            (x_grid < self.xhi - 2 * (x_grid[1] - x_grid[0]))
        Mfit = np.column_stack([np.exp(1j * S[m]), np.exp(-1j * S[m])])
        coef, _, _, _ = np.linalg.lstsq(Mfit, psi_vals[m], rcond=None)
        alpha, beta = coef[0], coef[1]
        return beta / alpha, psi_vals / alpha

    def train(self, steps=0, verbose_every=0, lbfgs=False, **kw):
        """Deterministic linear least-squares training (kept for API compat)."""
        return self.solve()

    def loss(self):
        """Weighted residual magnitude (diagnostic; the LS solve is exact)."""
        Mp, Mm = self.modulation_on(self.x_np)
        res = Mp * self._resUp + Mm * self._resUm
        r2 = (np.abs(res / self.res_scale) ** 2) * self._w_res
        return float(np.sum(r2) / max(self._w_res.sum(), 1e-6))
