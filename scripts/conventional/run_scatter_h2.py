"""A+BC = H + H2(rigid, v=0) fixed-energy scattering on the BKMP2 PES.
2D (R, gamma) model, J=0-like:

  H = -BR d2_R + (BR/R^2 + Brot) j^2(gamma) + V(re, R, gamma)
  E_tot = E_c + V_H2(re)
  BVP: chi(R_min)=0 ; (d_R - iK) chi(R_max) = -2 i K e^{-i K R_max}

.. note::
   **Legacy / failure-analysis script.**  This is the ORIGINAL plain Close-coupling (channel-projected)
   ansatz (incoming wave gated at a fixed R_gate, soft Robin penalty) that
   produced the 100-156% relative L2 errors documented in the paper
   (Helmholtz phase-slip + deep-well wavelength compression).
   For the corrected, < 1% solver use ``scripts/run_scatter_server.py``
   (channel-coupled scheme) or ``scripts/run_scatter_3d.py`` (full 3D).

FD reference: complex sparse Kronecker system (Robin with source).
Close-coupling (channel-projected): incoming wave built into the ansatz, wall gate, Robin/Dirichlet
penalties; the BKMP2 PES is precomputed once on the quadrature.
"""
import sys
import time
from pathlib import Path
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.compat import get_torch
from src.reference import LegendreAngularDVR, _fd_second, gauss_legendre_1d
from src import bkmp2

torch = get_torch()                     # Close-coupling (channel-projected) solvers need PyTorch
torch.set_num_threads(3)

BOHR = bkmp2.BOHR
re_H2_A = 0.74144                       # H2 equilibrium (Angstrom)
BR = 16.857633 / 0.671883               # cm^-1 A^2 (mu_R = 2mH/3)
Brot = 16.857633 / (0.503967 * re_H2_A**2)   # H2 rotational constant, cm^-1
M_H = 1.007825

bkmp2.get_lib(verbose=True)             # locate or auto-build the library


def V_2d(R_A, gamma):
    """BKMP2 with the H2 bond frozen at re; (R, gamma) -> Eh (broadcast)."""
    R_a = np.asarray(R_A, float)
    G_a = np.asarray(gamma, float)
    Rb, Gg = np.broadcast_arrays(R_a, G_a)
    Rb = Rb / BOHR
    re_b = re_H2_A / BOHR
    d2 = np.sqrt(Rb ** 2 + re_b ** 2 / 4 + Rb * re_b * np.cos(Gg))
    d3 = np.sqrt(Rb ** 2 + re_b ** 2 / 4 - Rb * re_b * np.cos(Gg))
    out = np.empty(d2.shape, dtype=float)
    lib = bkmp2.get_lib(verbose=False)
    it = np.nditer([d2, d3, out], [], [['readonly']] * 2 + [['readwrite']])
    for b_, c_, z in it:
        z[...] = bkmp2.V_h3_bohr(re_b, float(b_), float(c_), lib=lib)
    return out

V_asym = float(V_2d(30.0, 0.0))          # isolated H2 at re (R=30 A), Eh
print(f'BKMP2: V_H2(re) = {V_asym*219474.63:.1f} cm^-1, Brot = {Brot:.2f} cm^-1',
      flush=True)

def fd_scatter(E_c, R_lo=0.7, R_hi=9.0, n_R=450, n_gamma=60):
    R = np.linspace(R_lo, R_hi, n_R); h = R[1]-R[0]
    ang = LegendreAngularDVR(n_gamma); xg = ang.x
    gamma = np.arccos(xg)
    RR, GG = np.meshgrid(R, gamma, indexing='ij')
    Vg = V_2d(RR, GG) * 219474.63137   # Eh -> cm^-1
    E_tot = E_c + V_asym*219474.63
    K = float(np.sqrt(E_c / BR))
    H = sp.kron(-BR*_fd_second(n_R, h, 4), sp.identity(n_gamma), format='csr')
    C = BR/R**2 + Brot
    H += sp.kron(sp.diags(C), sp.csr_matrix(ang.j2), format='csr')
    H += sp.diags(Vg.ravel())
    A = (H - E_tot*sp.identity(H.shape[0])).astype(complex).tolil()
    b = np.zeros(n_R*n_gamma, dtype=complex)
    A[0, :] = 0.0; A[0, 0] = 1.0; b[0] = 0.0
    j = n_R-1; coef = np.array([3.0, -4.0, 1.0])/(2.0*h)
    for g in range(n_gamma):
        row = j*n_gamma + g; A[row, :] = 0.0
        A[row, j*n_gamma+g] = coef[0] - 1j*K
        A[row, (j-1)*n_gamma+g] = coef[1]
        A[row, (j-2)*n_gamma+g] = coef[2]
        b[row] = -2j*K*np.exp(-1j*K*R[j])
    chi = spla.spsolve(A.tocsc(), b).reshape(n_R, n_gamma)
    return chi, R, xg, K, E_tot

class Pinns:
    """Close-coupling (channel-projected) for the H+H2 2D scattering BVP (complex, quadrature-cached V)."""
    def __init__(self, E_c, R_lo=0.7, R_hi=9.0, R_gate=1.6, gate_w=0.6,
                 n_quad=(250, 50), hidden=(3, 64), n_fourier=32, seed=0):
        self.E_c = float(E_c); self.E_tot = E_c + V_asym*219474.63
        self.K = float(np.sqrt(E_c / BR))
        self.R_lo, self.R_hi, self.R_gate, self.gate_w = R_lo, R_hi, R_gate, gate_w
        torch.manual_seed(seed)
        self.nf = n_fourier
        mods, d = [], 2 + 4*n_fourier
        for hh in hidden:
            lin = torch.nn.Linear(d, hh, dtype=torch.float64)
            torch.nn.init.xavier_normal_(lin.weight, gain=1.6)
            mods += [lin, torch.nn.Tanh()]; d = hh
        out = torch.nn.Linear(d, 2, dtype=torch.float64)
        torch.nn.init.xavier_normal_(out.weight, gain=0.3)
        self.net = torch.nn.Sequential(*mods)
        ptsR, wR = gauss_legendre_1d(R_lo, R_hi, n_quad[0])
        ptsG, wG = gauss_legendre_1d(0.0, np.pi, n_quad[1])
        RR, GG = np.meshgrid(ptsR, ptsG, indexing='ij')
        self.pts = [RR.ravel(), GG.ravel()]
        self.w = torch.tensor((wR[:, None]*wG[None, :]).ravel(),
                              dtype=torch.float64)
        Vq = V_2d(self.pts[0], self.pts[1]) * 219474.63137
        self.Vq = torch.tensor(Vq.ravel(), dtype=torch.float64)
        self.x = [torch.tensor(p, dtype=torch.float64,
                               requires_grad=True).unsqueeze(-1)
                  for p in self.pts]
        self.Rg = self.x[0].squeeze(-1); self.Gg = self.x[1].squeeze(-1)
        self.scale = abs(self.E_tot) + BR*self.K**2
        self.bc_w = 30.0
        # Robin boundary collocation
        xb, _ = gauss_legendre_1d(0.0, np.pi, 80)
        self.xb = [torch.tensor(np.full(len(xb), R_hi), dtype=torch.float64,
                                requires_grad=True).unsqueeze(-1),
                   torch.tensor(xb, dtype=torch.float64,
                                device='cpu').unsqueeze(-1)]

    def _feat(self, xs):
        f = [2.0*(t-lo)/(hi-lo)-1.0 for t, lo, hi in
             zip(xs, (self.R_lo, 0.0), (self.R_hi, np.pi))]
        ms = torch.arange(1, self.nf+1, dtype=torch.float64)*np.pi
        for q in f:
            f += [torch.sin(ms*q), torch.cos(ms*q)]
        return torch.cat(f, dim=-1)

    def chi(self, xs):
        R, gam = xs[0].squeeze(-1), xs[1].squeeze(-1)
        u = self.net(self._feat(xs)); u = u[:, 0] + 1j*u[:, 1]
        W = torch.sigmoid((R - self.R_gate)/self.gate_w)
        return torch.exp(-1j*self.K*R)*(1.0-W) + (R-self.R_lo)*u*W

    def residual(self, xs):
        p = self.chi(xs); re, im = p.real, p.imag
        parts = []
        for comp in (re, im):
            d2R = torch.autograd.grad(comp.sum(), xs[0], create_graph=True)[0]
            d2R = torch.autograd.grad(d2R.sum(), xs[0], create_graph=True)[0]
            g1 = torch.autograd.grad(comp.sum(), xs[1], create_graph=True)[0]
            g2 = torch.autograd.grad(g1.sum(), xs[1], create_graph=True)[0]
            cot = torch.cos(xs[1])/torch.sin(xs[1])
            j2 = -(g2 + cot*g1).squeeze(-1)
            d2R = d2R.squeeze(-1)
            C = BR/self.Rg**2 + Brot
            parts.append(-BR*d2R + C*j2 + self.Vq*comp - self.E_tot*comp)
        return parts[0] + 1j*parts[1]

    def loss(self):
        res = self.residual(self.x)
        bulk = torch.mean((res/self.scale).abs()**2)
        p0 = self.chi([torch.tensor([[self.R_lo]], dtype=torch.float64,
                                    requires_grad=True), self.xb[1][:1]])
        Rb = self.xb[0]
        pb = self.chi(self.xb)
        dR = torch.autograd.grad(pb.sum(), Rb, create_graph=True)[0].squeeze(-1)
        rob = dR - 1j*self.K*pb + 2j*self.K*torch.exp(-1j*self.K*self.R_hi)
        bc = (p0.abs()**2).squeeze() + (rob.abs()**2).mean()
        loss = bulk + self.bc_w*bc/self.scale**2
        return loss, float(bulk)

    def params(self):
        return list(self.net.parameters())

    def train(self, steps=3000, lr0=1.5e-3, lbfgs=True):
        opt = torch.optim.Adam(self.params(), lr=lr0)
        sch = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=steps,
                                                         eta_min=lr0/50)
        for it in range(steps):
            opt.zero_grad(); l, _ = self.loss(); l.backward()
            opt.step(); sch.step()
        if lbfgs:
            best, bst, stall, done = np.inf, None, 0, 0
            while done < 500 and stall < 4:
                opt = torch.optim.LBFGS(self.params(), lr=0.15, max_iter=50,
                                        history_size=60,
                                        line_search_fn='strong_wolfe')
                def cl():
                    opt.zero_grad(); l, _ = self.loss(); l.backward()
                    return l
                opt.step(cl); done += 50
                l, _ = self.loss()
                if not np.isfinite(l): break
                if l < best-1e-8:
                    best, stall = l, 0
                    bst = [t.detach().clone() for t in self.params()]
                else:
                    stall += 1
            if bst is not None:
                with torch.no_grad():
                    for t, b_ in zip(self.params(), bst): t.copy_(b_)

    def evaluate(self, R_grid, xg_grid):
        RR, GG = np.meshgrid(R_grid, xg_grid, indexing='ij')
        xs = [torch.tensor(RR.ravel(), dtype=torch.float64,
                           requires_grad=True).unsqueeze(-1),
              torch.tensor(GG.ravel(), dtype=torch.float64).unsqueeze(-1)]
        with torch.no_grad():
            p = self.chi(xs).cpu().numpy()
        return p.reshape(RR.shape)

if __name__ == '__main__':
    OUT = Path(__file__).resolve().parents[2]/'results'; OUT.mkdir(exist_ok=True)
    rows = []
    for E_c in [400.0, 800.0]:
        t0 = time.time()
        chi_ref, R, xg, K, E_tot = fd_scatter(E_c)
        print(f'E_c={E_c:.0f} FD: K={K:.3f} asym|chi|={np.abs(chi_ref[R>8]).mean():.3f} '
              f'max={np.abs(chi_ref).max():.2f} ({time.time()-t0:.0f}s)', flush=True)
        pinn = Pinns(E_c)
        pinn.train(steps=3000, lr0=1.5e-3, lbfgs=True)
        ps = pinn.evaluate(R, xg)
        l2 = float(np.sqrt(np.sum(np.abs(ps-chi_ref)**2))
                   / np.sqrt(np.sum(np.abs(chi_ref)**2)))
        rows.append(dict(E_c=E_c, K=K, l2=l2))
        print(f'E_c={E_c:.0f}: Close-coupling (channel-projected) L2 = {l2:.4f}', flush=True)
        np.savez(OUT/f'_h2scatter_{int(E_c)}.npz', R=R, xg=xg,
                 chi_ref=chi_ref, chi_cc=ps)
    import json
    json.dump(rows, open(OUT/'scatter_h2_summary.json', 'w'), indent=2)
    print('done', flush=True)
