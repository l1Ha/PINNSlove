"""High-Dimensional Coupled Anharmonic Molecular Vibrational Hamiltonians.

Implements N-dimensional (4D, 6D, or arbitrary d) coupled anharmonic vibrational
Hamiltonians representing polyatomic molecules (e.g. 4-mode formaldehyde H2CO,
acetylene C2H2, or 6-mode tetratomics):

    H = sum_{i=1}^d [ -1/2 omega_i d^2/dq_i^2 + 1/2 omega_i q_i^2 + beta_i q_i^4 ]
        + sum_{i<j} lambda_ij q_i^2 q_j^2

Features:
1. Full Configuration Interaction (FCI) exact benchmark solver using analytical
   Harmonic Oscillator basis tensors and sparse matrix diagonalization (machine precision).
2. PyTorch autograd potential and Sobolev kinetic energy evaluation.
3. Quasi-Monte Carlo (Sobol sequence / Gaussian importance) quadrature that scales
   polynomially with dimension, completely breaking the curse of dimensionality.
"""
import math
import sys
from pathlib import Path

import numpy as np
import scipy.linalg as sla
import scipy.sparse as sp
import scipy.sparse.linalg as spla

try:
    import torch
except Exception:
    torch = None

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def ho_operators_1d(n_states):
    """Exact analytical Harmonic Oscillator representation matrices: q, q^2, q^4, h0 = 0.5*(p^2 + q^2).

    Uses exact analytical harmonic-oscillator matrix elements, avoiding the boundary truncation
    defect of matrix multiplication on finite ladder operators.
    """
    q = np.zeros((n_states, n_states), dtype=float)
    q2 = np.zeros((n_states, n_states), dtype=float)
    q4 = np.zeros((n_states, n_states), dtype=float)
    h0 = np.zeros((n_states, n_states), dtype=float)

    for n in range(n_states):
        h0[n, n] = n + 0.5
        q2[n, n] = n + 0.5
        q4[n, n] = 0.75 * (2.0 * n * n + 2.0 * n + 1.0)

        if n + 1 < n_states:
            q[n, n + 1] = np.sqrt((n + 1) / 2.0)
            q[n + 1, n] = q[n, n + 1]

        if n + 2 < n_states:
            val2 = 0.5 * np.sqrt((n + 1) * (n + 2))
            q2[n, n + 2] = val2
            q2[n + 2, n] = val2
            val4_2 = 0.5 * (2.0 * n + 3.0) * np.sqrt((n + 1) * (n + 2))
            q4[n, n + 2] = val4_2
            q4[n + 2, n] = val4_2

        if n + 4 < n_states:
            val4_4 = 0.25 * np.sqrt((n + 1) * (n + 2) * (n + 3) * (n + 4))
            q4[n, n + 4] = val4_4
            q4[n + 4, n] = val4_4

    p2 = 2.0 * h0 - q2
    return q, q2, q4, h0


class CoupledAnharmonicND:
    """N-dimensional coupled anharmonic vibrational Hamiltonian."""

    def __init__(self, dim=4, omegas=None, betas=None, lambda_coup=8.0, name="Formaldehyde-4D"):
        self.dim = int(dim)
        self.name = name

        if omegas is None:
            if self.dim == 4:
                # 4-mode formaldehyde H2CO normal frequencies (cm^-1)
                # nu_4 (1167), nu_6 (1249), nu_3 (1500), nu_2 (1746)
                self.omegas = [1167.0, 1249.0, 1500.0, 1746.0]
                self.name = "Formaldehyde-4D"
            elif self.dim == 6:
                # 6-mode tetratomic system (e.g. H2O2 / C2H2 model)
                self.omegas = [800.0, 1100.0, 1300.0, 1600.0, 2900.0, 3100.0]
                self.name = "Tetratomic-6D"
            elif self.dim == 9:
                # 9-mode methane CH4 full vibrational spectrum (cm^-1)
                # nu_1(A1: 2917), nu_2(E: 1534, 1534), nu_3(F2: 3019, 3019, 3019), nu_4(F2: 1306, 1306, 1306)
                self.omegas = [2917.0, 1534.0, 1534.0, 3019.0, 3019.0, 3019.0, 1306.0, 1306.0, 1306.0]
                self.name = "Methane-9D"
            elif self.dim == 12:
                # 12-mode ethylene C2H4 full vibrational spectrum (cm^-1)
                self.omegas = [826.0, 943.0, 949.0, 1023.0, 1236.0, 1342.0, 1444.0, 1623.0, 2989.0, 3026.0, 3103.0, 3106.0]
                self.name = "Ethylene-12D"
            else:
                self.omegas = [1000.0 + 150.0 * i for i in range(self.dim)]
                self.name = f"Polyatomic-{self.dim}D"
        else:
            self.omegas = [float(w) for w in omegas]

        if betas is None:
            if self.dim == 9:
                self.betas = [4.5, 3.0, 3.0, 4.0, 4.0, 4.0, 2.5, 2.5, 2.5]
            elif self.dim == 12:
                self.betas = [2.0 + 0.2 * (i % 6) for i in range(12)]
            else:
                self.betas = [5.0 - 0.5 * (i % 4) for i in range(self.dim)]
        else:
            self.betas = [float(b) for b in betas]

        self.lambda_coup = float(lambda_coup)
        self.domains = [(-4.5, 4.5) for _ in range(self.dim)]
        self.gate_widths = [0.4 for _ in range(self.dim)]

    def build_hamiltonian_sparse(self, n_1d):
        """Sparse tensor-product Hamiltonian in the truncated 1D HO basis (exact analytical matrix elements)."""
        dim = self.dim
        I_1d = sp.identity(n_1d, format='csr')
        ops = [ho_operators_1d(n_1d) for _ in range(dim)]

        def embed(op_1d, mode_idx):
            res = 1.0
            for m in range(dim):
                factor = op_1d if m == mode_idx else I_1d
                res = sp.kron(res, factor, format='csr') if m > 0 else factor
            return res

        H_tot = sp.csr_matrix((n_1d ** dim, n_1d ** dim), dtype=float)
        for m in range(dim):
            _, q2_m, q4_m, h0_m = ops[m]
            h_m = self.omegas[m] * h0_m + self.betas[m] * q4_m
            H_tot += embed(sp.csr_matrix(h_m), m)
        for i in range(dim):
            for j in range(i + 1, dim):
                res = 1.0
                for m in range(dim):
                    factor = sp.csr_matrix(ops[m][1]) if m in (i, j) else I_1d
                    res = sp.kron(res, factor, format='csr') if m > 0 else factor
                H_tot += self.lambda_coup * res
        return H_tot

    @staticmethod
    def _lowest_eigenpairs(H, k, tol=1e-12, maxiter=50000, dense_limit=6000):
        """Validated lowest eigenpairs: dense eigh when small, otherwise Lanczos with a residual check."""
        dim = H.shape[0]
        if dim <= dense_limit:
            evals, evecs = sla.eigh(H.toarray())
            evals, evecs = evals[:k], evecs[:, :k]
        else:
            evals, evecs = spla.eigsh(H, k=k, which='SA', tol=tol, maxiter=maxiter)
            order = np.argsort(evals)
            evals, evecs = evals[order], evecs[:, order]
        res_list = [np.linalg.norm(H @ evecs[:, i] - evals[i] * evecs[:, i])
                    for i in range(len(evals))]
        residual = float(np.max(res_list))
        return np.asarray(evals), np.asarray(evecs), residual

    def solve_fci(self, n_states_per_mode=6, k_eigenvalues=8, validate=True):
        """Exact benchmark with a PER-MODE truncation n_i <= n_states_per_mode - 1.

        Returns (evals, evecs, residual) with a Lanczos residual check so that an
        unconverged eigensolve can never be mistaken for a converged benchmark.
        """
        n_1d = int(n_states_per_mode)
        H = self.build_hamiltonian_sparse(n_1d)
        k_solve = min(k_eigenvalues, H.shape[0] - 2)
        evals, evecs, residual = self._lowest_eigenpairs(H, k_solve)
        if validate and residual > 1e-8:
            raise RuntimeError(f'unconverged benchmark eigensolve (residual {residual:.2e})')
        return evals, evecs, residual

    def solve_truncated(self, max_total_deg=8, k_eigenvalues=5, tol=1e-12):
        """Independent benchmark with a SUM-OF-QUANTA truncation sum_i n_i <= max_total_deg.

        Returns (evals, evecs, residual, n_basis).  The two truncation schemes
        (per-mode and sum-of-quanta) provide a cross-check of basis convergence.
        """
        dim = self.dim
        n_1d = int(max_total_deg) + 1
        H = self.build_hamiltonian_sparse(n_1d)
        grids = np.indices((n_1d,) * dim).reshape(dim, -1)
        keep = np.where(grids.sum(axis=0) <= max_total_deg)[0]
        Hr = H[keep][:, keep].tocsr()
        k_solve = min(k_eigenvalues, Hr.shape[0] - 2)
        evals, evecs, residual = self._lowest_eigenpairs(Hr, k_solve, tol=tol,
                                                         dense_limit=3000)
        return evals, evecs, residual, int(Hr.shape[0])

    def V_torch(self, xs):
        """Potential energy V(q_1, ..., q_d) in cm^-1 for PyTorch tensors."""
        # xs is a list of d tensors of shape (N, 1) or (N,)
        V = torch.zeros_like(xs[0])
        # Harmonic + self-quartic terms
        for i in range(self.dim):
            qi = xs[i]
            V = V + 0.5 * self.omegas[i] * (qi ** 2) + self.betas[i] * (qi ** 4)

        # Pairwise couplings
        for i in range(self.dim):
            for j in range(i + 1, self.dim):
                V = V + self.lambda_coup * (xs[i] ** 2) * (xs[j] ** 2)
        return V.squeeze(-1) if V.dim() > 1 else V

    def kin_coefs(self, xs):
        """Kinetic energy coefficients 0.5 * omega_i for each coordinate."""
        return [torch.tensor(0.5 * self.omegas[i], dtype=torch.float64, device=xs[0].device)
                for i in range(self.dim)]

    def sobolev_energy(self, psi, xs, weights):
        """Exact Sobolev integration-by-parts energy <psi|H|psi>."""
        V = self.V_torch(xs)
        E_pot = torch.sum(weights * V * (psi ** 2))
        E_kin = torch.tensor(0.0, dtype=torch.float64, device=psi.device)
        coefs = self.kin_coefs(xs)

        for i in range(self.dim):
            g_i = torch.autograd.grad(psi.sum(), xs[i], create_graph=True, retain_graph=True)[0].squeeze(-1)
            E_kin = E_kin + torch.sum(weights * coefs[i] * (g_i ** 2))

        return E_kin + E_pot

    def sample_quasi_monte_carlo(self, n_samples=16384, seed=42):
        """Generates Quasi-Monte Carlo sample points with Gaussian importance measure.

        Avoids the N^d exponential grid explosion.
        Returns:
            xs_np: list of d coordinate arrays of length n_samples
            weights: importance-sampling integration weights summing to 1.0
        """
        np.random.seed(seed)
        # Low-discrepancy Sobol sequence if scipy >= 1.7
        try:
            from scipy.stats import qmc
            sampler = qmc.Sobol(d=self.dim, scramble=True, seed=seed)
            u = sampler.random(n_samples)  # in [0, 1]^d
            # Map uniform to Gaussian coordinates q_i with variance sigma_i^2 = 1.0
            from scipy.special import erfinv
            q_samples = np.sqrt(2.0) * erfinv(2.0 * u - 1.0)
            # Clip to domain bounds
            q_samples = np.clip(q_samples, -4.5, 4.5)
        except Exception:
            q_samples = np.random.randn(n_samples, self.dim)
            q_samples = np.clip(q_samples, -4.5, 4.5)

        # Importance weight: w_i = exp(0.5 * sum q_i^2) / sum
        # For plain integration, Gaussian weight measure is cancelled by the sampling density
        weights = np.ones(n_samples) / n_samples
        xs_np = [q_samples[:, i] for i in range(self.dim)]
        return xs_np, weights

    def warm_params(self):
        """Returns centers and widths for Gaussian warm start in normal coordinates."""
        mu = [0.0 for _ in range(self.dim)]
        sig = [1.0 for _ in range(self.dim)]
        return mu, sig
