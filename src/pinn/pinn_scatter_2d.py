"""Physics-informed neural network for the 2D A+BC scattering BVP (scattered wave).

Genuine torch neural network trained on the PDE residual.  Formulation
identical to the FD reference (src/conventional/scatter_2d.py,
assemble_2d_scatter): with time convention e^{+iEt} the full wave is
chi = e^{+iKR} + chi_sc (unit isotropic incoming VALUE amplitude) and
chi_sc solves

    (H + CAP - E_tot) chi_sc = -(H - E_tot) chi_in
    chi_sc(R_lo) = chi_sc(R_hi) = 0,

where (H - E) chi_in = (E_c + V - E_tot) e^{+iKR} is localized to the
interaction region and the CAP absorber -i*eta*s^2 (R > cap_start) removes
the deep-well box resonance.  The ansatz chi_sc = (R - R_lo)(R_hi - R) u(x)
satisfies both Dirichlet conditions exactly, so the loss is the bulk PDE
residual only.  The FD reference is converged ONLY with the CAP active
(run_cap_scan.py); accuracy claims are meaningful only for that case.
"""
import numpy as np
import torch

from src.conventional.hamiltonians import Water2DProblem
from src.conventional.reference import gauss_legendre_1d


class ScatteringPINN2D:
    """PINN for the scattered 2D wave chi_sc defined above."""
    def __init__(self, pes, E_c, R_lo=0.6, R_hi=9.0,
                 n_quad=(300, 60), hidden=(3, 64), n_fourier=32, seed=0,
                 device='cpu', cap_start=None, cap_eta=None, cap_width=2.0):
        self.pes = pes
        self.E_c = float(E_c)
        self.E_tot = E_c + float(pes.V_valence(pes.re, pes.re, pes.theta_e)) + pes.De
        self.R_lo, self.R_hi = R_lo, R_hi
        # CAP absorber, identical to the FD reference (assemble_2d_scatter)
        self.cap_start, self.cap_eta, self.cap_width = cap_start, cap_eta, cap_width
        self.K = float(np.sqrt(E_c / pes.BR))
        self.problem = Water2DProblem(pes, R_lo=R_lo, R_hi=R_hi, gate_width=0.3)
        self.device = device
        torch.manual_seed(seed)
        self.n_fourier = n_fourier
        d_in = 2 + 2 * 2 * n_fourier
        mods, d = [], d_in
        for h in hidden:
            lin = torch.nn.Linear(d, h, dtype=torch.float64)
            torch.nn.init.xavier_normal_(lin.weight,
                                         gain=torch.nn.init.calculate_gain('tanh'))
            mods += [lin, torch.nn.Tanh()]
            d = h
        out = torch.nn.Linear(d, 2, dtype=torch.float64)
        torch.nn.init.xavier_normal_(out.weight, gain=0.3)
        mods.append(out)
        self.net = torch.nn.Sequential(*mods).to(device)

        pts, w = self.problem.quadrature(n_quad)
        self.pts_np = [np.ascontiguousarray(p, float) for p in pts]
        self.w = torch.tensor(w, dtype=torch.float64, device=device)
        self.lo = np.array([d[0] for d in self.problem.domains])
        self.hi = np.array([d[1] for d in self.problem.domains])
        self.x_train = [torch.tensor(p, dtype=torch.float64, device=device,
                                     requires_grad=True).unsqueeze(-1)
                        for p in self.pts_np]
        self.scale = abs(self.E_tot) + pes.BR * self.K ** 2

    def _features(self, xs):
        f = [2.0 * (t - lo) / (hi - lo) - 1.0
             for t, lo, hi in zip(xs, self.lo, self.hi)]
        ms = torch.arange(1, self.n_fourier + 1, dtype=torch.float64,
                          device=self.device) * np.pi
        base = list(f)          # iterate over the coordinate features only;
        for q in base:          # appending to f while looping over it made the
            f += [torch.sin(ms * q), torch.cos(ms * q)]   # list grow forever
        return torch.cat(f, dim=-1)

    def chi_sc(self, xs):
        R = xs[0].squeeze(-1)
        u = self.net(self._features(xs))
        u = u[:, 0] + 1j * u[:, 1]
        # hard Dirichlet at both walls: chi_sc(R_lo) = chi_sc(R_hi) = 0
        return (R - self.R_lo) * (self.R_hi - R) * u

    def _residual(self, xs):
        p = self.chi_sc(xs)
        Rv = xs[0].squeeze(-1)
        accs = []
        for comp in (p.real, p.imag):
            d2R = torch.autograd.grad(comp.sum(), xs[0], create_graph=True)[0]
            d2R = torch.autograd.grad(d2R.sum(), xs[0], create_graph=True)[0]
            g1 = torch.autograd.grad(comp.sum(), xs[1], create_graph=True)[0]
            g2 = torch.autograd.grad(g1.sum(), xs[1], create_graph=True)[0]
            cot = torch.cos(xs[1]) / torch.sin(xs[1])
            j2 = -(g2 + cot * g1).squeeze(-1)
            d2R = d2R.squeeze(-1)
            C = self.pes.BR / Rv ** 2 + self.pes.Br / self.pes.re ** 2
            V = self.problem.V_jacobi_torch(xs)
            if self.cap_eta:
                s = torch.clamp((Rv - self.cap_start) / self.cap_width,
                                min=0.0)
                V = V - 1j * self.cap_eta * s ** 2
            accs.append(-self.pes.BR * d2R + C * j2 + V * comp
                        - self.E_tot * comp)
        # source: (H - E) chi_in = (E_c + V - E_tot) e^{+iKR}, with the
        # PHYSICAL V (no CAP: the incoming wave is given analytically)
        V_phys = self.problem.V_jacobi_torch(xs)
        src = (self.E_c + V_phys - self.E_tot) * torch.exp(
            1j * self.K * Rv)
        return accs[0] + 1j * accs[1] + src

    def loss(self):
        res = self._residual(self.x_train)
        bulk = torch.mean((res / self.scale).abs() ** 2)
        return bulk, float(bulk)

    def train(self, steps=3000, lr0=1.5e-3, verbose_every=1000, lbfgs=True):
        opt = torch.optim.Adam(self._all_params(), lr=lr0)
        sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=steps,
                                                           eta_min=lr0 / 50)
        for it in range(steps):
            opt.zero_grad()
            loss, bulk = self.loss()
            loss.backward()
            opt.step()
            sched.step()
            if verbose_every and (it % verbose_every == 0 or it == steps - 1):
                print(f'  it {it}: bulk {bulk:.3e}  loss {float(loss):.3e}',
                      flush=True)
        if lbfgs:
            best, best_state, stall, done = np.inf, None, 0, 0
            while done < 500 and stall < 4:
                opt = torch.optim.LBFGS(self._all_params(), lr=0.15,
                                        max_iter=50, history_size=60,
                                        line_search_fn='strong_wolfe')
                def closure():
                    opt.zero_grad()
                    l, _ = self.loss()
                    l.backward()
                    return l
                opt.step(closure)
                done += 50
                l, _ = self.loss()
                if not np.isfinite(float(l)):
                    break
                if float(l) < best - 1e-8:
                    best, stall = float(l), 0
                    best_state = [t.detach().clone() for t in self._all_params()]
                else:
                    stall += 1
            if best_state is not None:
                with torch.no_grad():
                    for t, b_ in zip(self._all_params(), best_state):
                        t.copy_(b_)

    def _all_params(self):
        return list(self.net.parameters())

    def evaluate(self, R_grid, xg_grid):
        RR, GG = np.meshgrid(R_grid, xg_grid, indexing='ij')
        xs = [torch.tensor(RR.ravel(), dtype=torch.float64,
                           device=self.device).unsqueeze(-1),
              torch.tensor(GG.ravel(), dtype=torch.float64,
                           device=self.device).unsqueeze(-1)]
        with torch.no_grad():
            p = self.chi_sc(xs).cpu().numpy()
        return p.reshape(RR.shape)
