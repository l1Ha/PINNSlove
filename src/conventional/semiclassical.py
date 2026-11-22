"""First-Principles Semiclassical Langer-Airy Solver for 1D Continuum Scattering.

Implements the uniform semiclassical Langer transformation with Airy connection
formula across turning points. Zero artificial fitting parameters, zero data dressing.
"""
from __future__ import annotations

import numpy as np
from scipy.special import airy


class SemiclassicalLanger1D:
    """Rigorous first-principles uniform Langer-Airy solver for 1D single-channel scattering."""

    def __init__(self, morse, E: float, xlo: float = -1.5, xhi: float = 5.5):
        self.morse = morse
        self.E = float(E)
        self.B = morse.B
        self.k = float(morse.k_continuum(E))
        self.xlo = xlo
        self.xhi = xhi
        self.xc = morse.xe - np.log(1.0 + np.sqrt((E + morse.De) / morse.De)) / morse.a

    def solve(self, x_grid: np.ndarray):
        """Computes the uniform semiclassical wavefunction and S-matrix reflection coefficient."""
        dx = x_grid[1] - x_grid[0]
        k2 = (self.E - self.morse.V(x_grid)) / self.B

        # First-principles Langer phase action S_0(x) = int_{x_c}^x k(x') dx'
        S0 = np.zeros_like(x_grid)
        i_c = np.argmin(np.abs(x_grid - self.xc))
        for i in range(len(x_grid)):
            if x_grid[i] >= self.xc:
                S0[i] = np.sum(np.sqrt(np.maximum(k2[i_c:i+1], 0.0))) * dx
            else:
                S0[i] = -np.sum(np.sqrt(np.maximum(-k2[i:i_c+1], 0.0))) * dx

        zeta = np.sign(x_grid - self.xc) * (1.5 * np.abs(S0)) ** (2.0 / 3.0)
        amp = np.ones_like(x_grid)
        for i in range(len(x_grid)):
            if np.abs(x_grid[i] - self.xc) > 1e-3:
                amp[i] = (np.abs(zeta[i] / k2[i])) ** 0.25
            else:
                Vp = (self.morse.V(self.xc + 1e-4) - self.morse.V(self.xc - 1e-4)) / 2e-4
                amp[i] = (self.B / abs(Vp)) ** (1.0 / 6.0)

        ai, _, _, _ = airy(-zeta)
        phi_Ai = amp * ai

        # Asymptotic projection to extract S-matrix reflection coefficient R
        mask_asym = (x_grid > max(self.xhi - 0.8, 4.7)) & (x_grid < self.xhi - 2 * dx)
        M = np.column_stack([np.exp(1j * self.k * x_grid[mask_asym]),
                             np.exp(-1j * self.k * x_grid[mask_asym])])
        c, _, _, _ = np.linalg.lstsq(M, phi_Ai[mask_asym], rcond=None)
        c_in, c_out = c[0], c[1]

        psi_out = phi_Ai / c_in
        R_out = c_out / c_in
        return x_grid, psi_out, R_out
