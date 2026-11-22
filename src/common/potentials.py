"""Potentials: 1D Morse oscillator and a model water (H2O) PES.

The water PES is written in valence coordinates (r1, r2, theta) as
    V = sum_i De(1-exp(-a (ri-re)))^2 - 2 De                (two Morse OH stretches)
        + g(r1) g(r2) [ 1/2 k_theta (theta-theta_e)^2 + c4 (theta-theta_e)^4 ]
        + F12 (r1-re)(r2-re) exp(-b12[(r1-re)^2+(r2-re)^2])
        + f_rtheta (dr1+dr2)(theta-theta_e) g(r1) g(r2),
with g(r) = exp(-b_g (r-re)^2) damping the angular terms so that the bending
force vanishes when either OH bond is broken (H + OH free-rotation limit).
Parameters are calibrated so the vibrational spectrum is close to real water
(see scripts/calibrate_pes.py).  Energy zero: O + H + H at infinite separation.

PyTorch is optional: the NumPy path (used by all classical reference solvers)
imports and runs without it; the torch paths raise a clear error if it is absent.
"""
import numpy as np

try:                                    # optional dependency
    import torch
except Exception:                       # pragma: no cover
    torch = None

from src.common.units import M_O, M_H, kinetic_factor


def _backend(x):
    """Return the array module matching ``x`` (numpy for ndarrays, torch otherwise)."""
    if isinstance(x, np.ndarray) or np.isscalar(x):
        return np
    if torch is not None and isinstance(x, torch.Tensor):
        return torch
    raise TypeError(
        f"Unsupported array type {type(x)!r}; expected numpy.ndarray or torch.Tensor."
    )


class Morse:
    """Morse oscillator V(x) = De (1 - exp(-a (x - xe)))^2 - De.

    x may be a bond length in Angstrom; the potential minimum is -De at x = xe
    and the dissociation asymptote is 0.  The bound spectrum is exact:
        E_n = -De + we(n+1/2) - we_xe (n+1/2)^2,
        we = 2 a sqrt(De * B),  we_xe = B a^2,  B = hbar^2/(2 mu).
    """

    def __init__(self, De_cm1, a_A1, xe_A, mass_amu):
        self.De = float(De_cm1)
        self.a = float(a_A1)
        self.xe = float(xe_A)
        self.mass = float(mass_amu)
        self.B = float(kinetic_factor(self.mass))       # cm^-1 A^2
        self.we = 2.0 * self.a * np.sqrt(self.De * self.B)
        self.we_xe = self.B * self.a ** 2
        # number of bound levels: n < lambda - 1/2, lambda = sqrt(De/(B a^2))
        self.lambda_ = np.sqrt(self.De / (self.B * self.a ** 2))
        self.n_bound = int(np.floor(self.lambda_ - 0.5 + 1e-9)) + 1

    def V(self, x):
        x = np.asarray(x, dtype=float)
        return self.De * (1.0 - np.exp(-self.a * (x - self.xe))) ** 2 - self.De

    def level(self, n):
        n = np.asarray(n, dtype=float)
        return -self.De + self.we * (n + 0.5) - self.we_xe * (n + 0.5) ** 2

    def k_continuum(self, E_cm1):
        """Asymptotic wavenumber for E > 0 (right-hand open channel):
        k = sqrt(E / (hbar^2/2mu))."""
        return np.sqrt(np.asarray(E_cm1) / self.B)


class WaterPES:
    """Model water PES (valence coordinates) with Jacobi-coordinate interface."""

    def __init__(self, re=0.95785, theta_e_deg=104.50, De=40570.0, a=2.27,
                 k_theta=34200.0, c4=0.0, F12=-1600.0, b12=2.0,
                 b_g=2.5, f_rtheta=0.0):
        self.re = float(re)
        self.theta_e = np.deg2rad(theta_e_deg)
        self.De = float(De)
        self.a = float(a)
        self.k_theta = float(k_theta)
        self.c4 = float(c4)
        self.F12 = float(F12)
        self.b12 = float(b12)
        self.b_g = float(b_g)
        self.f_rtheta = float(f_rtheta)
        # pairwise Jacobi reduced masses (diatom = O-H1, atom = H2)
        self.mu_r = M_O * M_H / (M_O + M_H)
        self.mu_R = M_H * (M_O + M_H) / (M_O + 2.0 * M_H)
        self.Br = kinetic_factor(self.mu_r)     # cm^-1 A^2
        self.BR = kinetic_factor(self.mu_R)
        # Jacobi geometry at equilibrium
        alpha = M_H / (M_O + M_H)
        self.alpha = alpha
        self.R_e, self.gamma_e = self._jacobi_geom()

    def _jacobi_geom(self):
        re, th = self.re, self.theta_e
        com1 = alpha = self.alpha * re
        R2 = np.array([re * np.cos(th) - com1, re * np.sin(th)])
        R_e = np.linalg.norm(R2)
        gamma_e = np.arccos(np.clip(R2[0] / R_e, -1, 1))
        return R_e, gamma_e

    def V_min(self):
        return -2.0 * self.De

    # ---- valence coordinates -------------------------------------------
    def V_valence(self, r1, r2, theta):
        """Works with numpy arrays or torch tensors (no grad needed for V)."""
        m = _backend(r1)
        dr1, dr2 = r1 - self.re, r2 - self.re
        dth = theta - self.theta_e
        g1 = m.exp(-self.b_g * dr1 ** 2)
        g2 = m.exp(-self.b_g * dr2 ** 2)
        Vmorse = (self.De * (1.0 - m.exp(-self.a * dr1)) ** 2 - self.De
                  + self.De * (1.0 - m.exp(-self.a * dr2)) ** 2 - self.De)
        vbend = 0.5 * self.k_theta * dth ** 2 + self.c4 * dth ** 4
        v12 = self.F12 * dr1 * dr2 * m.exp(-self.b12 * (dr1 ** 2 + dr2 ** 2))
        vrt = self.f_rtheta * (dr1 + dr2) * dth * g1 * g2
        return Vmorse + g1 * g2 * vbend + v12 + vrt

    # ---- Jacobi coordinates --------------------------------------------
    def to_valence(self, r, R, gamma):
        """(r, R, gamma) -> (r1, r2, theta); numpy or torch (broadcast)."""
        m = _backend(r)
        arc = self.alpha * r
        hx = arc + R * m.cos(gamma)           # x of H2 relative to O
        hy = R * m.sin(gamma)
        r2 = m.hypot(hx, hy)
        if m is np:
            cos_t = np.clip(hx / r2, -1.0, 1.0)
        else:
            cos_t = m.clamp(hx / r2, -1.0, 1.0)
        return r, r2, m.arccos(cos_t)

    def V_jacobi(self, r, R, gamma):
        r1, r2, th = self.to_valence(r, R, gamma)
        return self.V_valence(r1, r2, th)
