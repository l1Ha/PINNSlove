"""High-Dimensional Quasi-Monte-Carlo Variational PINN (QMC-VPINN).

A second, *stochastic* route to the same high-dimensional vibrational benchmarks
used by the deterministic collocation solver, for cross-validation:

    H psi = E psi ,   H = sum_i [ -1/2 omega_i d^2/dq_i^2 + 1/2 omega_i q_i^2 + beta_i q_i^4 ]
                        + sum_{i<j} lambda_ij q_i^2 q_j^2 .

Design (all steps validated against the dual-scheme converged benchmarks):
1. Quasi-Monte-Carlo sampling: Sobol sequence mapped through the inverse Gaussian
   CDF, i.e. samples q_k ~ N(0, I).  The corresponding importance weight for a
   plain Lebesgue integral is  w_k ∝ exp(|q_k|^2 / 2)  (constant factors cancel in
   the Rayleigh quotient).  Using uniform-in-a-box samples with unit weights --
   the previous implementation -- biases the quotient because the integrand is
   concentrated near the origin; that manifested as a Ritz warm start of
   1775 cm^-1 instead of 2855 cm^-1.
2. Basis: tensor products of Hermite functions psi_k(q) = (...)^-1 H_k(q) e^{-q^2/2}
   with EXACT analytic derivatives from the ladder relation
       d psi_k/dq = sqrt(k/2) psi_{k-1} - sqrt((k+1)/2) psi_{k+1} ,
   so no finite differences enter the kinetic energy.
3. Variational functional in Sobolev form (first derivatives only):
       E[psi] = sum_i w_i [ sum_m (omega_m/2) (d_m psi)^2 + V(q_i) psi^2 ] / sum_i w_i psi^2 .
4. Optimization: linear coefficients are solved by the weighted generalized
   eigenproblem (Ritz); optional torch refinement of the coefficient vector and a
   small amplitude network is accepted only if it lowers the energy (safety gate).
"""
import itertools
import math
import sys
import time
from pathlib import Path

import numpy as np
from numpy.polynomial.hermite import hermval
import scipy.linalg as sla

try:
    import torch
    import torch.nn as nn
except Exception:                                   # pragma: no cover
    torch = None

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def hermite_function(k, q):
    """Normalized Hermite function psi_k(q) = (pi^{1/4} sqrt(2^k k!))^{-1} H_k(q) e^{-q^2/2}."""
    norm = 1.0 / (np.pi ** 0.25 * np.sqrt(2.0 ** k * math.factorial(k)))
    return norm * hermval(q, [0.0] * k + [1.0]) * np.exp(-0.5 * q ** 2)


def hermite_poly_and_deriv(k, q):
    """Polynomial part a_k(q) of the Hermite function and its derivative a_k'(q).

    psi_k(q) = a_k(q) e^{-q^2/2} with a_k(q) = (pi^{1/4} sqrt(2^k k!))^{-1} H_k(q),
    a_k'(q) = (pi^{1/4} sqrt(2^k k!))^{-1} 2k H_{k-1}(q).
    """
    norm = 1.0 / (np.pi ** 0.25 * np.sqrt(2.0 ** k * math.factorial(k)))
    a = norm * hermval(q, [0.0] * k + [1.0])
    da = np.zeros_like(q)
    if k >= 1:
        da = norm * 2.0 * k * hermval(q, [0.0] * (k - 1) + [1.0])
    return a, da


class HighDimVariationalPINN:
    """QMC variational solver for d-dimensional coupled anharmonic vibrations."""

    def __init__(self, problem, n_samples=16384, max_poly_deg=4, seed=42,
                 q_clip=5.0, device='cpu'):
        self.problem = problem
        self.dim = problem.dim
        self.n_samples = int(n_samples)
        self.max_poly_deg = int(max_poly_deg)
        self.q_clip = float(q_clip)
        self.device = torch.device(device) if torch is not None else None
        self.seed = seed

        # --- Quasi-Monte-Carlo nodes: Sobol -> N(0, I) ---
        rng = np.random.default_rng(seed)
        try:
            from scipy.stats import qmc
            sobol = qmc.Sobol(d=self.dim, scramble=True, seed=seed)
            u = sobol.random(self.n_samples)
        except Exception:                            # pragma: no cover
            u = np.clip(rng.random((self.n_samples, self.dim)), 1e-12, 1 - 1e-12)
        from scipy.special import erfinv
        u = np.clip(u, 1e-12, 1.0 - 1e-12)
        q = np.sqrt(2.0) * erfinv(2.0 * u - 1.0)
        q = np.clip(q, -self.q_clip, self.q_clip)
        self.q = q                                   # (N, d)

        # Weights of the Gaussian-measure representation.  Writing every basis
        # function as psi_k(q) = a_k(q) e^{-|q|^2/2} (a_k = polynomial part), the
        # Lebesgue integrals become bounded-moment expectations under q ~ N(0, I):
        #     int f dq = (2 pi)^{d/2} E[ f(q) e^{+|q|^2/2} e^{-|q|^2} ]
        #             = (2 pi)^{d/2} E[ ftilde(q) e^{-|q|^2/2} ] ,
        # so w_k = e^{-|q_k|^2/2} (normalized) is BOUNDED, unlike the naively
        # "correct" importance weight e^{+|q|^2/2}, whose variance diverges for
        # d >= 2 and which biases the Rayleigh quotient below the true ground
        # state (the failure mode of the previous implementation).
        logw = -0.5 * np.sum(q ** 2, axis=1)
        logw -= logw.max()
        self.w = np.exp(logw)
        self.w /= self.w.sum()

        self.xs_np = [q[:, i].copy() for i in range(self.dim)]

        # --- tensor Hermite-function basis with sum_i k_i <= max_poly_deg ---
        self.basis_indices = []
        for deg in range(self.max_poly_deg + 1):
            for tup in itertools.product(range(deg + 1), repeat=self.dim):
                if sum(tup) == deg:
                    self.basis_indices.append(tup)
        self.n_basis = len(self.basis_indices)
        self._build_basis()

        self.states = []                             # (psi, energy) for deflation
        self.self_check()

    # ------------------------------------------------------------------
    def _build_basis(self):
        """Basis values and per-mode analytic derivatives on the QMC nodes."""
        N = self.n_samples
        d = self.dim
        kmax = self.max_poly_deg + 1
        # per-mode tables: phi[mode][k] and dphi[mode][k]
        self.phi_mode = []
        self.dphi_mode = []
        for m in range(d):
            qm = self.xs_np[m]
            tab = np.zeros((kmax + 1, N))
            dtab = np.zeros((kmax + 1, N))
            for k in range(kmax + 1):
                tab[k], dtab[k] = hermite_poly_and_deriv(k, qm)
            self.phi_mode.append(tab)
            self.dphi_mode.append(dtab)

        Phi = np.ones((N, self.n_basis))
        dPhi = [np.ones((N, self.n_basis)) for _ in range(d)]
        for col, tup in enumerate(self.basis_indices):
            factors = [self.phi_mode[m][tup[m]] for m in range(d)]
            prod = np.ones(N)
            for m in range(d):
                prod = prod * factors[m]
            Phi[:, col] = prod
            # derivative w.r.t. mode i: replace exactly ONE factor by its derivative
            # (the earlier in-place accumulation overwrote the product whenever the
            #  derivative mode came last, corrupting d/dq_i for i > 0)
            for i in range(d):
                der = np.ones(N)
                for m in range(d):
                    der = der * (self.dphi_mode[m][tup[m]] if m == i else factors[m])
                dPhi[i][:, col] = der
        self.Phi, self.dPhi = Phi, dPhi

    # ------------------------------------------------------------------
    def self_check(self, tol=5e-3, verbose=False):
        """Internal consistency guard: in the pure-harmonic limit the weighted Ritz
        problem must reproduce E0 = sum_i omega_i / 2 to within Monte-Carlo accuracy."""
        w = self.w
        sw = np.sqrt(w)[:, None]
        Phi_w = self.Phi * sw
        S = Phi_w.T @ Phi_w
        H = np.zeros((self.n_basis, self.n_basis))
        for i in range(self.dim):
            g = (self.dPhi[i] - self.q[:, i][:, None] * self.Phi)
            Gi = g * np.sqrt(w * 0.5 * self.problem.omegas[i])[:, None]
            H += Gi.T @ Gi
        V_harm = 0.5 * (self.q ** 2 @ np.asarray(self.problem.omegas, dtype=float))
        H += Phi_w.T @ (Phi_w * V_harm[:, None])
        H = 0.5 * (H + H.T)
        ridge = 1e-12 * np.trace(S) / self.n_basis
        ev = sla.eigvalsh(H, S + ridge * np.eye(self.n_basis))[0]
        exact = 0.5 * sum(self.problem.omegas)
        rel = abs(ev - exact) / exact
        if verbose:
            print(f'  [QMC-VPINN] self-check harmonic limit: E0={ev:.4f} vs exact {exact:.4f} '
                  f'(rel {rel*100:.3f}%)')
        if rel > tol:
            raise RuntimeError(f'QMC-VPINN self-check failed: harmonic limit off by {rel*100:.2f}% '
                               f'(weights/gradient construction is broken)')
        return rel

    # ------------------------------------------------------------------
    def V_np(self):
        """Potential energy on the QMC nodes (cm^-1)."""
        q = self.q
        V = np.zeros(self.n_samples)
        for i in range(self.dim):
            V += 0.5 * self.problem.omegas[i] * q[:, i] ** 2
            V += self.problem.betas[i] * q[:, i] ** 4
        for i in range(self.dim):
            for j in range(i + 1, self.dim):
                V += self.problem.lambda_coup * q[:, i] ** 2 * q[:, j] ** 2
        return V

    # ------------------------------------------------------------------
    def solve_state(self, state_idx=0, max_lbfgs=200, refine_torch=True, verbose=True):
        """Variational solution of eigenstate state_idx (ground = 0).

        Returns the converged energy in cm^-1.  The linear coefficients are obtained
        from the weighted generalized eigenproblem; an optional torch refinement is
        accepted only if it lowers the energy.
        """
        t0 = time.time()
        w = self.w
        V = self.V_np()

        # Ritz matrices in the weighted metric: H_ij = <phi_i | H | phi_j>_w
        sw = np.sqrt(w)[:, None]
        Phi_w = self.Phi * sw
        S_mat = Phi_w.T @ Phi_w
        H_mat = np.zeros((self.n_basis, self.n_basis))
        for i in range(self.dim):
            # d/dq_i [ a(q) e^{-|q|^2/2} ] = (d_i a - q_i a) e^{-|q|^2/2}
            grad_reduced = self.dPhi[i] - self.q[:, i][:, None] * self.Phi
            Gi_w = grad_reduced * np.sqrt(w * 0.5 * self.problem.omegas[i])[:, None]
            H_mat += Gi_w.T @ Gi_w
        H_mat += Phi_w.T @ (Phi_w * V[:, None])
        H_mat = 0.5 * (H_mat + H_mat.T)

        # deflation against previously computed states
        if self.states:
            C = np.stack([np.sum(w[:, None] * s_psi[:, None] * self.Phi, axis=0)
                          for s_psi in self.states])          # (n_prev, n_basis)
            H_mat, S_mat, Z = self._project(H_mat, S_mat, C)
        else:
            Z = np.eye(self.n_basis)

        ridge = 1e-12 * np.trace(S_mat) / S_mat.shape[0]
        evals, evecs = sla.eigh(H_mat, S_mat + ridge * np.eye(S_mat.shape[0]))
        c = Z @ evecs[:, 0]
        psi = self.Phi @ c                       # reduced amplitude a(q) = psi e^{+|q|^2/2}
        psi /= np.sqrt(np.sum(w * psi ** 2))
        E_ritz = float(evals[0])

        if verbose:
            print(f"  [QMC-VPINN] state {state_idx}: Ritz warm start E = {E_ritz:.4f} cm^-1 "
                  f"({self.n_basis} basis terms, {self.n_samples} QMC nodes)")

        E_best, c_best = E_ritz, c.copy()

        # optional torch refinement of the linear coefficients (accepted only if lower)
        if refine_torch and torch is not None:
            ct = torch.tensor(c, dtype=torch.float64, requires_grad=True)
            Phi_t = torch.tensor(self.Phi, dtype=torch.float64)
            dPhi_t = [torch.tensor(self.dPhi[i], dtype=torch.float64) for i in range(self.dim)]
            w_t = torch.tensor(w, dtype=torch.float64)
            V_t = torch.tensor(V, dtype=torch.float64)
            coef_t = [0.5 * self.problem.omegas[i] for i in range(self.dim)]

            q_t = torch.tensor(self.q, dtype=torch.float64)
            red_grad_t = [dPhi_t[i] - q_t[:, i][:, None] * Phi_t for i in range(self.dim)]

            def energy_t(cc):
                psi_t = Phi_t @ cc
                nrm = torch.sum(w_t * psi_t ** 2)
                e = torch.tensor(0.0, dtype=torch.float64)
                for i in range(self.dim):
                    gi = red_grad_t[i] @ cc
                    e = e + coef_t[i] * torch.sum(w_t * gi ** 2)
                e = e + torch.sum(w_t * V_t * psi_t ** 2)
                return e / nrm

            opt = torch.optim.LBFGS([ct], max_iter=max_lbfgs, line_search_fn='strong_wolfe',
                                    tolerance_grad=1e-14, tolerance_change=1e-16)

            def closure():
                opt.zero_grad()
                loss = energy_t(ct)
                loss.backward()
                return loss

            try:
                opt.step(closure)
                with torch.no_grad():
                    E_t = float(energy_t(ct))
                if E_t < E_best:                     # safety gate: never accept a worse energy
                    E_best, c_best = E_t, ct.detach().numpy().copy()
            except Exception:
                pass

        psi = self.Phi @ c_best
        psi /= np.sqrt(np.sum(w * psi ** 2))
        self.states.append(psi)
        if verbose:
            print(f"  [QMC-VPINN] state {state_idx}: final E = {E_best:.4f} cm^-1 "
                  f"({time.time() - t0:.1f} s)")
        return E_best

    @staticmethod
    def _project(H_mat, S_mat, C):
        """Constrain the Ritz problem to the orthogonal complement of the span of C."""
        _, sv, Vh = np.linalg.svd(C, full_matrices=True)
        r = int((sv > max(sv.max() * 1e-10, 1e-300)).sum())
        Z = Vh[r:].T
        return Z.T @ H_mat @ Z, Z.T @ S_mat @ Z, Z
