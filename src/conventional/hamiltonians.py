"""PINN problem definitions: quadrature grids + autograd Hamiltonians.

Conventions (all energies cm^-1, lengths Angstrom, angles radian):
* 1D Morse          : H = -B d2/dx2 + V(x),  plain measure dx.
* 2D water (R,gamma): reduced chi = R psi (r frozen at r_e),
      H = -BR d2_R + (BR/R^2 + Br/r_e^2) j^2 + V(r_e, R, gamma),
      j^2 = -(d2_gamma + cot(gamma) d_gamma),  measure dR sin(gamma) dgamma.
* 3D water (r,R,gamma): reduced chi = r R psi,
      H = -Br d2_r - BR d2_R + (Br/r^2 + BR/R^2) j^2 + V(r, R, gamma),
      measure dr dR sin(gamma) dgamma.

Each problem provides
  * quadrature(n)   -> Gauss-Legendre nodes (flattened grid) and weights with
                       the sin(gamma) measure folded in,
  * h_psi           -> H psi via autograd (second derivatives; used for
                       validation and the scattering residual),
  * sobolev_energy  -> <psi|H|psi> in integration-by-parts form (first
                       derivatives only; exact because psi vanishes at the
                       box edges and is regular at the gamma poles),
  * warm_params()   -> centre/widths of a Gaussian warm-start basis.
"""
import numpy as np
import torch

from src.conventional.reference import gauss_legendre_1d


def _second(y, x):
    """d^2 y / d x^2 for a leaf tensor x of shape (N,1) and y of shape (N,)."""
    g = torch.autograd.grad(y.sum(), x, create_graph=True, retain_graph=True)[0]
    return torch.autograd.grad(g.sum(), x, create_graph=True, retain_graph=True)[0].squeeze(-1)


def _first(y, x):
    g = torch.autograd.grad(y.sum(), x, create_graph=True, retain_graph=True)[0]
    return g.squeeze(-1)


def _hessian_diag(pes, xs0, step=1e-3):
    """Numerical diagonal of the PES Hessian (numpy) at a Jacobi point."""
    xs0 = [float(v) for v in xs0]
    H = np.zeros(len(xs0))
    for i in range(len(xs0)):
        e = [0.0] * len(xs0)
        e[i] = step
        xp = [v + e[j] for j, v in enumerate(xs0)]
        xm = [v - e[j] for j, v in enumerate(xs0)]
        H[i] = (pes.V_jacobi(*xp) - 2 * pes.V_jacobi(*xs0)
                + pes.V_jacobi(*xm)) / e[i] ** 2
    return H


class MorseProblem1D:
    """Bound-state problem for a Morse oscillator on a finite box."""
    dim = 1
    regularity = ""

    def __init__(self, morse, xlo=-1.6, xhi=4.6, gate_width=0.25):
        self.morse = morse
        self.domains = [(xlo, xhi)]
        self.gate_widths = [gate_width]
        self.E_floor = -morse.De
        self.E_scale = morse.De

    def quadrature(self, n_per_dim):
        x, w = gauss_legendre_1d(*self.domains[0], n_per_dim)
        return [x], w

    def V_torch(self, xs):
        (x,) = xs
        mo = self.morse
        return (mo.De * (1.0 - torch.exp(-mo.a * (x - mo.xe))) ** 2
                - mo.De).squeeze(-1)

    def h_psi(self, xs, psi):
        (x,) = xs
        d2 = _second(psi, x)
        return -self.morse.B * d2 + self.V_torch(xs) * psi

    def sobolev_energy(self, xs, psi, w):
        (x,) = xs
        gx = _first(psi, x)
        return torch.sum(w * (self.morse.B * gx ** 2
                              + self.V_torch(xs) * psi ** 2))

    def kin_coefs(self, xs):
        (x,) = xs
        return [torch.full_like(x.squeeze(-1), self.morse.B)]

    def warm_params(self):
        mo = self.morse
        return (np.array([mo.xe]),
                np.array([(mo.B / (mo.De * mo.a ** 2)) ** 0.25]))


class _WaterBase:
    def V_jacobi_torch(self, xs):
        xs_v = xs if len(xs) == 3 else [xs[0] * 0 + self.pes.re] + list(xs)
        r1, r2, th = self.pes.to_valence(*xs_v)
        return self.pes.V_valence(r1, r2, th).squeeze(-1)

    V_torch = V_jacobi_torch

    def _angular(self, chi, gamma):
        """j^2 chi (without hbar^2) = -(chi_gg + cot(gamma) chi_g)."""
        g1 = torch.autograd.grad(chi.sum(), gamma, create_graph=True,
                                 retain_graph=True)[0]
        g2 = torch.autograd.grad(g1.sum(), gamma, create_graph=True,
                                 retain_graph=True)[0]
        cot = torch.cos(gamma) / torch.sin(gamma)
        return -(g2 + cot * g1).squeeze(-1)


class Water2DProblem(_WaterBase):
    """Frozen-bond 2D water problem in (R, gamma), r = r_e, reduced chi = R psi."""
    dim = 2
    regularity = "R"

    def __init__(self, pes, R_lo=0.55, R_hi=3.2, gate_width=0.25):
        self.pes = pes
        self.domains = [(R_lo, R_hi), (0.0, np.pi)]
        self.gate_widths = [gate_width, None]   # no gate in gamma
        self.E_floor = pes.V_min()
        self.E_scale = 1.0e4

    def quadrature(self, n_per_dim):
        nR, ng = n_per_dim
        tR, wR = gauss_legendre_1d(*self.domains[0], nR)
        tg, wg = gauss_legendre_1d(0.0, np.pi, ng)
        w = wR[:, None] * wg[None, :] * np.sin(tg)[None, :]
        pts = [np.broadcast_to(tR[:, None], (nR, ng)).ravel(),
               np.broadcast_to(tg[None, :], (nR, ng)).ravel()]
        return pts, w.ravel()

    def h_psi(self, xs, psi):
        R, gamma = xs
        d2R = _second(psi, R)
        j2 = self._angular(psi, gamma)
        C = self.pes.BR / R.squeeze(-1) ** 2 + self.pes.Br / self.pes.re ** 2
        return -self.pes.BR * d2R + C * j2 + self.V_jacobi_torch(xs) * psi

    def sobolev_energy(self, xs, psi, w):
        R, gamma = xs
        gR = _first(psi, R)
        gg = _first(psi, gamma)
        C = self.pes.BR / R.squeeze(-1) ** 2 + self.pes.Br / self.pes.re ** 2
        return torch.sum(w * (self.pes.BR * gR ** 2 + C * gg ** 2
                              + self.V_jacobi_torch(xs) * psi ** 2))

    def kin_coefs(self, xs):
        R, gamma = xs
        C = self.pes.BR / R.squeeze(-1) ** 2 + self.pes.Br / self.pes.re ** 2
        return [torch.full_like(C, self.pes.BR), C]

    def warm_params(self):
        pe = self.pes
        _, kR, kg = _hessian_diag(pe, [pe.re, pe.R_e, pe.gamma_e])
        sigR = (2 * pe.BR / kR) ** 0.25
        sigg = (2 * (pe.BR / pe.R_e ** 2 + pe.Br / pe.re ** 2) / kg) ** 0.25
        return np.array([pe.R_e, pe.gamma_e]), np.array([sigR, sigg])


class Water3DProblem(_WaterBase):
    """Full 3D water problem in (r, R, gamma), reduced chi = r R psi."""
    dim = 3
    regularity = "rR"

    def __init__(self, pes, r_lo=0.55, r_hi=2.4, R_lo=0.55, R_hi=3.0,
                 gate_width=0.22):
        self.pes = pes
        self.domains = [(r_lo, r_hi), (R_lo, R_hi), (0.0, np.pi)]
        self.gate_widths = [gate_width, gate_width, None]
        self.E_floor = pes.V_min()
        self.E_scale = 1.0e4

    def quadrature(self, n_per_dim):
        nr, nR, ng = n_per_dim
        tr, wr = gauss_legendre_1d(*self.domains[0], nr)
        tR, wR = gauss_legendre_1d(*self.domains[1], nR)
        tg, wg = gauss_legendre_1d(0.0, np.pi, ng)
        w = (wr[:, None, None] * wR[None, :, None] * wg[None, None, :]
             * np.sin(tg)[None, None, :])
        pts = [np.broadcast_to(tr[:, None, None], (nr, nR, ng)).ravel(),
               np.broadcast_to(tR[None, :, None], (nr, nR, ng)).ravel(),
               np.broadcast_to(tg[None, None, :], (nr, nR, ng)).ravel()]
        return pts, w.ravel()

    def h_psi(self, xs, psi):
        r, R, gamma = xs
        d2r = _second(psi, r)
        d2R = _second(psi, R)
        j2 = self._angular(psi, gamma)
        C = (self.pes.Br / r.squeeze(-1) ** 2
             + self.pes.BR / R.squeeze(-1) ** 2)
        return (-self.pes.Br * d2r - self.pes.BR * d2R + C * j2
                + self.V_jacobi_torch(xs) * psi)

    def sobolev_energy(self, xs, psi, w):
        r, R, gamma = xs
        gr = _first(psi, r)
        gR = _first(psi, R)
        gg = _first(psi, gamma)
        C = (self.pes.Br / r.squeeze(-1) ** 2
             + self.pes.BR / R.squeeze(-1) ** 2)
        return torch.sum(w * (self.pes.Br * gr ** 2 + self.pes.BR * gR ** 2
                              + C * gg ** 2 + self.V_jacobi_torch(xs) * psi ** 2))

    def kin_coefs(self, xs):
        r, R, gamma = xs
        C = (self.pes.Br / r.squeeze(-1) ** 2
             + self.pes.BR / R.squeeze(-1) ** 2)
        return [torch.full_like(C, self.pes.Br),
                torch.full_like(C, self.pes.BR), C]

    def warm_params(self):
        pe = self.pes
        kr, kR, kg = _hessian_diag(pe, [pe.re, pe.R_e, pe.gamma_e])
        sig_r = (2 * pe.Br / kr) ** 0.25
        sig_R = (2 * pe.BR / kR) ** 0.25
        sig_g = (2 * (pe.Br / pe.re ** 2 + pe.BR / pe.R_e ** 2) / kg) ** 0.25
        return (np.array([pe.re, pe.R_e, pe.gamma_e]),
                np.array([sig_r, sig_R, sig_g]))
