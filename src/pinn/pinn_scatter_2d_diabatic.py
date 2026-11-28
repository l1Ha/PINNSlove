"""Diabatic-channel carrier PINN for the 2D scattered wave.

Theory: docs/theory_scattering_2d.md section 12.  Subclasses CarrierPINN2D and
replaces the adiabatic carrier columns with diabatic ones: the parallel
transport (Eq. 27) annihilates the derivative couplings exactly, the carriers
are built on the effective diabatic curves Veff_jj = (U^T diag(eps) U)_jj -
B_R (U^T A^T A U)_jj (Eq. 30, the O(A^2) correction kept), and the residual of
each column is purely the multiplicative off-diagonal coupling (Eq. 33).
Everything else (same discrete operator, cross-validated ridge, wall window,
honest fallback) is inherited unchanged.
"""
from __future__ import annotations

import numpy as np

from src.pinn.pinn_scatter_2d_carrier import CarrierPINN2D
from src.conventional.diabatic_carrier import (DiabaticBasis,
                                               build_langer_carriers)


class DiabaticCarrierPINN2D(CarrierPINN2D):
    def __init__(self, pes, E_c, R_lo=0.6, R_hi=12.0, n_R=400, n_gamma=80,
                 cap_start=5.0, cap_eta=100.0, cap_width=2.0,
                 n_channels=8, n_corr=0, sigma_R=(8.0, 20.0), sigma_g=6.0,
                 seed=0, val_frac=0.2, wall_window=0.6, wall_weight=1e3):
        super().__init__(pes, E_c, R_lo=R_lo, R_hi=R_hi, n_R=n_R,
                         n_gamma=n_gamma, cap_start=cap_start,
                         cap_eta=cap_eta, cap_width=cap_width,
                         n_channels=n_channels, n_corr=0,
                         sigma_R=sigma_R, sigma_g=sigma_g, seed=seed,
                         val_frac=val_frac, wall_window=wall_window,
                         wall_weight=wall_weight, include_particular=False)
        h = self.R[1] - self.R[0]
        self.dia = DiabaticBasis(self.carrier.Phi, self.carrier.eps, h,
                                 self.pes.BR,
                                 H_gamma=self.carrier.H_gamma)
        curves = self.dia.curves()
        order = np.argsort(curves[-1])[:n_channels]
        cols, names = [], []
        self.carriers_dia = {}
        for jj in order:
            d = build_langer_carriers(curves[:, jj], self.R, self.E_tot,
                                      self.pes.BR)
            if d is None:
                continue
            self.carriers_dia[jj] = d
            Phi = self.dia.Phi_dia[:, :, jj]
            cols.append((d['U_plus'][:, None] * Phi).ravel())
            names.append(f'dU+_{jj}')
            cols.append((d['U_minus'][:, None] * Phi).ravel())
            names.append(f'dU-_{jj}')
        self.B = (np.column_stack(cols) if cols
                  else np.zeros((self.n_R * self.n_gamma, 0), complex))
        self.names = names

    def column_residuals(self):
        """Per-column residual r_j = ||A Psi_j|| / (||Psi_j|| |E_tot|)."""
        Ab = self.A @ self.B
        scale = abs(self.E_tot)
        return {name: float(np.linalg.norm(Ab[:, k])
                            / np.linalg.norm(self.B[:, k]) / scale)
                for k, name in enumerate(self.names)}
