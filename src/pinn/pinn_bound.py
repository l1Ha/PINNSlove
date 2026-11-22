"""Variational PINN for bound states (Rayleigh-Ritz / Deep-Ritz style).

For state n the network ansatz is
    chi(x) = regularity(x) * gates(x) * NN_theta(x),
    psi    = chi / sqrt(<chi|chi>)          (hard normalisation),
and the training loss is
    L = ( <psi|H|psi>_IBP - E_floor ) / E_scale
        + ortho_weight * sum_{m<n} <psi|psi_m>^2 ,
where <.|.>_IBP is the integration-by-parts (Sobolev) form of the Rayleigh
quotient  sum_i B_i <(d_i psi)^2> + <V psi^2>  -- mathematically identical to
<psi|H|psi> for functions vanishing at the box edges, but needing only first
derivatives (no double backprop).  The normalisation constant is detached:
its gradient is purely radial, so removing it implements projected gradient
descent on the unit sphere.

Because a randomly initialised network makes the Rayleigh quotient descend
extremely slowly (spectral bias), each state is warm-started from a mini
Rayleigh-Ritz solve in a Gaussian x polynomial basis built from the local
PES curvature (no reference data involved): the Ritz vector n is used as a
pretraining target for the network, after which the variational descent only
refines it.  Orthogonality penalties to the converged lower states force the
sequential solve onto successive excited states.
"""
import numpy as np
import torch

from src.conventional.hamiltonians import _first


def _mlp(dim_in, hidden, seed, out_dim=1, arch='modified'):
    """Plain tanh MLP or the two-stream gated MLP of Wang, Teng & Perdikaris
    (2023), U = act(W_s z + b_s) * act(W_x z + b_x), which trains more stably
    than the plain tanh stack at equal parameter count.
    """
    torch.manual_seed(seed)
    if arch == 'plain':
        layers, d = [], dim_in
        for h in hidden:
            lin = torch.nn.Linear(d, h, dtype=torch.float64)
            torch.nn.init.xavier_normal_(
                lin.weight, gain=torch.nn.init.calculate_gain('tanh'))
            layers += [lin, torch.nn.Tanh()]
            d = h
        out = torch.nn.Linear(d, out_dim, dtype=torch.float64)
        torch.nn.init.xavier_normal_(out.weight, gain=0.5)
        layers.append(out)
        return torch.nn.Sequential(*layers)
    # modified (two-stream) gated MLP
    net = _ModifiedMLP(dim_in, hidden, out_dim)
    return net


class _GatePair(torch.nn.Module):
    """U = act(W_s z + b_s) * act(W_x z + b_x) (Wang, Teng & Perdikaris,
    'An Expert's Guide to Training PINNs', arXiv:2308.08468)."""
    def __init__(self, d_in, d_out):
        super().__init__()
        self.ws = torch.nn.Linear(d_in, d_out, dtype=torch.float64)
        self.wx = torch.nn.Linear(d_in, d_out, dtype=torch.float64)
        g = torch.nn.init.calculate_gain('tanh')
        torch.nn.init.xavier_normal_(self.ws.weight, gain=g)
        torch.nn.init.xavier_normal_(self.wx.weight, gain=g)

    def forward(self, z):
        return torch.tanh(self.ws(z)) * torch.tanh(self.wx(z))


class _ModifiedMLP(torch.nn.Sequential):
    def __init__(self, dim_in, hidden, out_dim):
        super().__init__()
        d = dim_in
        for h in hidden:
            self.add_module(f'gate{len(self)}', _GatePair(d, h))
            d = h
        out = torch.nn.Linear(d, out_dim, dtype=torch.float64)
        torch.nn.init.xavier_normal_(out.weight, gain=0.5)
        self.add_module('out', out)


class VariationalPINN:
    def __init__(self, problem, hidden=(3, 64), n_quad=400, seed=0,
                 device='cpu', ortho_weight=300.0, n_fourier=0,
                 gaussian_env=True, env_scale=1.0, arch='plain',
                 learn_fft=True, basis_scales=(1.0,), k_extra=6,
                 spectral_layer=False, spec_kmax=12):
        self.problem = problem
        self.device = device
        self.ortho_weight = ortho_weight
        self.n_fourier = n_fourier
        self.gaussian_env = gaussian_env
        self.env_scale = float(env_scale)
        self.basis_scales = tuple(basis_scales)
        self.k_extra = int(k_extra)
        self.spectral_layer = bool(spectral_layer)
        self.spec_kmax = int(spec_kmax)
        # count the spectral basis: 1 + d * kmax * n_scales
        n_spec = 1 + problem.dim * self.spec_kmax * len(self.basis_scales)
        self.spec_coeffs = (torch.nn.Parameter(
            torch.zeros(n_spec, dtype=torch.float64))
            if spectral_layer else None)
        dim_in = problem.dim + 2 * problem.dim * n_fourier
        self.net = _mlp(dim_in, hidden, seed, arch=arch).to(device)
        # learnable global frequency scale for the Fourier features
        # (Tancik et al. 2021; the scale is adapted during training)
        self.fft_log = (torch.nn.Parameter(torch.zeros(1, dtype=torch.float64))
                        if (learn_fft and n_fourier > 0) else None)
        if gaussian_env:
            mu0, sig0 = problem.warm_params()
            self.mu = torch.nn.Parameter(torch.tensor(mu0, dtype=torch.float64))
            self.log_sig = torch.nn.Parameter(
                torch.log(torch.tensor(sig0, dtype=torch.float64)))
        self.E_param = torch.nn.Parameter(torch.tensor(0.0, dtype=torch.float64))
        self._ritz_energy = None
        pts, w = problem.quadrature(n_quad)
        self.pts_np = [np.ascontiguousarray(p, dtype=float) for p in pts]
        self.w_np = w
        self.w = torch.tensor(w, dtype=torch.float64, device=device)
        self.lo = np.array([d[0] for d in problem.domains])
        self.hi = np.array([d[1] for d in problem.domains])
        self.states = []          # list of dict(psi, E, E_pde)

    # ------------------------------------------------------------------
    def _leaf_inputs(self, pts_list):
        return [torch.tensor(np.asarray(p, float), dtype=torch.float64,
                             device=self.device, requires_grad=True).unsqueeze(-1)
                for p in pts_list]

    def _scaled(self, x, lo, hi):
        return 2.0 * (x - lo) / (hi - lo) - 1.0

    def _features(self, xs):
        feats = [self._scaled(t, lo, hi) for t, lo, hi in zip(xs, self.lo, self.hi)]
        if self.n_fourier:
            ms = torch.arange(1, self.n_fourier + 1, dtype=torch.float64,
                              device=self.device) * np.pi
            if self.fft_log is not None:
                ms = ms * torch.exp(self.fft_log.clamp(-2.0, 2.0))
            fourier = []
            for f in feats:
                fourier += [torch.sin(ms * f), torch.cos(ms * f)]
            feats = feats + fourier
        return torch.cat(feats, dim=-1)

    def envelope(self, xs):
        """Multiplicative envelope: trainable Gaussians x regularity x gates."""
        pr = self.problem
        env = torch.ones_like(xs[0].squeeze(-1))
        if self.gaussian_env:
            sig = torch.exp(self.log_sig.clamp(-3.0, 1.5))
            for i, (x, lg) in enumerate(zip(xs, sig)):
                lo, hi = self.problem.domains[i]
                m = self.mu[i].clamp(lo + 0.05, hi - 0.05)
                env = env * torch.exp(-0.5 * ((x.squeeze(-1) - m) / lg) ** 2)
        if getattr(pr, 'regularity', '') == 'R':      # 2D: chi = R u
            env = env * xs[0].squeeze(-1)
        elif getattr(pr, 'regularity', '') == 'rR':   # 3D: chi = r R u
            env = env * xs[0].squeeze(-1) * xs[1].squeeze(-1)
        for t, (lo, hi), wd in zip(xs, pr.domains, pr.gate_widths):
            if wd is not None:
                x = t.squeeze(-1)
                env = env * torch.tanh((x - lo) / wd) * torch.tanh((hi - x) / wd)
        return env

    def _all_params(self):
        ps = list(self.net.parameters())
        if self.gaussian_env:
            ps += [self.mu, self.log_sig]
        if self.fft_log is not None:
            ps += [self.fft_log]
        if self.spec_coeffs is not None:
            ps += [self.spec_coeffs]
        ps += [self.E_param]
        return ps

    def _spectral_part(self, xs):
        """Linear spectral layer: sum_k c_k phi_k(x). The basis carries the
        same regularity factors and boundary gates as the ansatz, so the
        spectral part alone is admissible.  With u == 0 the ansatz reduces
        exactly to the (trainable-width) Gauss-Hermite Ritz subspace."""
        cols = self._basis_torch(xs, self._spec_mu, self._spec_sig,
                                 self.spec_kmax - self.k_extra,
                                 scales=self.basis_scales)
        return (self.spec_coeffs * torch.stack(cols, dim=-1)).sum(-1)

    def chi_with_env(self, xs):
        pr = self.problem
        if getattr(self, 'spec_coeffs', None) is None:
            return self.net(self._features(xs)).squeeze(-1) * self.envelope(xs)
        if not hasattr(self, '_spec_mu'):
            mu0, sig0 = pr.warm_params()
            sig0 = sig0 * self.env_scale
            self._spec_mu = torch.tensor(mu0, dtype=torch.float64,
                                         device=self.device)
            self._spec_sig = torch.tensor(sig0, dtype=torch.float64,
                                          device=self.device)
        chi = self.net(self._features(xs)).squeeze(-1) * self.envelope(xs)
        return chi + self._spectral_part(xs)

    def _norm_full_pool(self):
        with torch.no_grad():
            xs = self._leaf_inputs(self.pts_np)
            chi = self.chi_with_env(xs)
            return float(torch.sqrt(torch.sum(self.w * chi ** 2)))

    # ------------------------------------------------------------------
    def warm_start(self, n_state, n_sub=5000, pre_steps=600, lr=5e-3,
                   verbose=True):
        """Mini Rayleigh-Ritz in a Gaussian x polynomial basis; pretrain the
        network on Ritz vector n.  Uses only PES information."""
        pr = self.problem
        rng = np.random.default_rng(1234 + n_state)
        N = len(self.w_np)
        idx = np.sort(rng.choice(N, min(n_sub, N), replace=False))
        pts_sub = [p[idx] for p in self.pts_np]
        w_sub = self.w_np[idx]
        # fresh per-state envelope at the equilibrium geometry; the network
        # (Fourier features) handles any extra width of excited states
        if self.gaussian_env:
            mu0, sig0 = pr.warm_params()
            sig0 = sig0 * self.env_scale
            with torch.no_grad():
                self.mu.copy_(torch.tensor(mu0, dtype=torch.float64))
                self.log_sig.copy_(torch.log(torch.tensor(
                    sig0, dtype=torch.float64)))
            mu, sig = mu0, sig0
        else:
            mu, sig = pr.warm_params()
            sig = sig * 1.2
        xtil = [2.0 * (p - lo) / (hi - lo) - 1.0
                for p, lo, hi in zip(pts_sub, self.lo, self.hi)]
        G = np.exp(-0.5 * sum(((p - m) / s) ** 2
                              for p, m, s in zip(pts_sub, mu, sig)))
        # Hermite-polynomial x Gaussian basis: H_k((x_i-mu_i)/sigma_i) are
        # orthogonal under the Gaussian measure, keeping S well conditioned
        # (plain powers x^k become near-collinear for narrow Gaussians).
        from numpy.polynomial.hermite import hermval
        cols = [G]
        k_max = n_state + self.k_extra   # keep the constrained span rich
        for sc in self.basis_scales:
            Gs = np.exp(-0.5 * sum(((p - m) / (s * sc)) ** 2
                                   for p, m, s in zip(pts_sub, mu, sig)))
            if pr.dim == 1:
                t_i = (pts_sub[0] - mu[0]) / (sig[0] * sc)
                for k in range(1, k_max + 1):
                    cols.append(Gs * hermval(t_i, [0.0] * k + [1.0]))
            elif pr.dim == 2:
                t0 = (pts_sub[0] - mu[0]) / (sig[0] * sc)
                t1 = (pts_sub[1] - mu[1]) / (sig[1] * sc)
                for k1 in range(k_max + 1):
                    for k2 in range(k_max + 1 - k1):
                        if k1 == 0 and k2 == 0:
                            continue
                        cols.append(Gs * hermval(t0, [0.0] * k1 + [1.0])
                                       * hermval(t1, [0.0] * k2 + [1.0]))
            elif pr.dim == 3:
                k_3d = min(k_max, 5)
                t0 = (pts_sub[0] - mu[0]) / (sig[0] * sc)
                t1 = (pts_sub[1] - mu[1]) / (sig[1] * sc)
                t2 = (pts_sub[2] - mu[2]) / (sig[2] * sc)
                for k1 in range(k_3d + 1):
                    for k2 in range(k_3d + 1 - k1):
                        for k3 in range(k_3d + 1 - k1 - k2):
                            if k1 == 0 and k2 == 0 and k3 == 0:
                                continue
                            cols.append(Gs * hermval(t0, [0.0] * k1 + [1.0])
                                           * hermval(t1, [0.0] * k2 + [1.0])
                                           * hermval(t2, [0.0] * k3 + [1.0]))
            else:
                for i in range(len(mu)):
                    t_i = (pts_sub[i] - mu[i]) / (sig[i] * sc)
                    for k in range(1, k_max + 1):
                        cols.append(Gs * hermval(t_i, [0.0] * k + [1.0]))
        B = np.stack(cols, axis=1)

        # mini-Ritz (generalized): H v = E S v in the raw basis, with H in the
        # integration-by-parts form  H_ij = sum_p w[sum_c coef_c df_i df_j + V f_i f_j]
        import scipy.linalg as sla
        xs = self._leaf_inputs(pts_sub)
        idx_t = torch.tensor(idx, dtype=torch.long, device=self.device)
        cols_t = self._basis_torch(xs, mu, sig, n_state)      # torch, with graph
        n_b = len(cols_t)
        D = []
        for c in range(pr.dim):
            Dc = torch.stack([_first(f, xs[c]) for f in cols_t], dim=1)
            D.append(Dc.detach().cpu().numpy())
        Vp = pr.V_torch(xs).detach().cpu().numpy()
        coefs = [c.detach().cpu().numpy() for c in pr.kin_coefs(xs)]
        H = np.zeros((n_b, n_b))
        for c in range(pr.dim):
            Dc = D[c] * np.sqrt(w_sub * coefs[c])[:, None]
            H += Dc.T @ Dc
        Bw = B * np.sqrt(w_sub)[:, None]
        H += Bw.T @ (Bw * Vp[:, None])          # sum_p w_p V_p f_i f_j
        H = 0.5 * (H + H.T)
        S = (B * np.sqrt(w_sub)[:, None]).T @ (B * np.sqrt(w_sub)[:, None])
        # constrain the Ritz problem to functions orthogonal to the converged
        # lower states: project the coefficient space on null(C)
        if self.states:
            C = np.stack([w_sub * pm[idx] for pm in
                          (s['psi'] for s in self.states)]) @ B        # (n_prev, n_b)
            _, sv, Vh = np.linalg.svd(C)
            r = int((sv > max(sv.max() * 1e-10, 1e-30)).sum())
            Z = Vh[r:].T                    # orthonormal basis of null(C)
            Hz, Sz = Z.T @ H @ Z, Z.T @ S @ Z
        else:
            Hz, Sz, Z = H, S, np.eye(n_b)
        ridge = 1e-10 * max(float(np.trace(Sz)) / Sz.shape[0], 1e-30)
        try:
            evals, evecs = sla.eigh(Hz, Sz + ridge * np.eye(Sz.shape[0]))
        except np.linalg.LinAlgError:
            evals, evecs = sla.eigh(Hz + ridge*np.eye(Sz.shape[0])*max(1.0, abs(Hz).max()),
                                    Sz + 1e-8 * np.eye(Sz.shape[0]) * np.trace(Sz) / Sz.shape[0])
        if verbose:
            print(f'  mini-Ritz estimates (cm-1): '
                  f'{np.round(evals[:min(n_state + 3, n_b)], 1)}', flush=True)
        # in the constrained (deflated) problem the target state is the
        # LOWEST Ritz value: all lower true states are excluded by the
        # orthogonality constraints.
        i_use = 0 if self.states else min(n_state, len(evals) - 1)
        self._ritz_energy = float(evals[i_use])
        with torch.no_grad():
            self.E_param.fill_(evals[i_use])
        v_full = Z @ evecs[:, i_use]        # back to raw coefficients
        target = B @ v_full
        target /= np.sqrt(w_sub @ target ** 2)

        # Pretrain by fitting chi directly to the Ritz target.  In 3D the
        # r*R regularity factor times two Gaussians makes target/envelope
        # blow up near the box corners, which a smooth tanh net cannot fit;
        # therefore fit chi only where the envelope is appreciable and pull
        # the network's u to the interior median elsewhere (keeps the fit
        # well scaled and u smooth).
        tt = torch.tensor(target, dtype=torch.float64, device=self.device)
        with torch.no_grad():
            env_np = self.envelope(xs).cpu().numpy()
        mask_np = env_np > 0.02 * env_np.max()
        if mask_np.sum() < 10:
            mask_np = env_np > 0.005 * env_np.max()
        u_int = target[mask_np] / env_np[mask_np]
        u_out_med = float(np.median(u_int))
        # cap the dynamic range: a nodal excited-state target divided by a
        # narrow envelope spans orders of magnitude and breaks the fit
        umax = 30.0 * max(abs(u_out_med), 1e-3)
        u_int_c = np.clip(u_int, -umax, umax)
        loss_w = np.where(mask_np, 1.0, 0.1) * w_sub
        lw = torch.tensor(loss_w, dtype=torch.float64, device=self.device)
        mt = torch.tensor(mask_np, dtype=torch.bool, device=self.device)
        u_med_t = torch.tensor(u_out_med, dtype=torch.float64,
                               device=self.device)

        u_int_t = torch.tensor(u_int_c, dtype=torch.float64,
                               device=self.device)

        def pre_closure():
            opt.zero_grad()
            chi = self.chi_with_env(xs)
            loss = torch.sum(lw * (chi - tt) ** 2)
            u = self.net(self._features(xs)).squeeze(-1)
            loss = loss + 0.1 * lw.sum() / max(float(mt.sum()), 1) * \
                torch.sum(lw[~mt] * ((u - u_med_t) ** 2)[~mt])
            loss.backward()
            return loss

        # Analytic starting point: constant u = median of the interior
        # target values.  This is a well-conditioned Gaussian trial function;
        # the network pretrain is only kept if it lowers the Rayleigh
        # quotient (the masked fit can otherwise degrade the start).
        last = [m for m in self.net.modules()
                if isinstance(m, torch.nn.Linear)][-1]
        with torch.no_grad():
            last.weight.zero_()
            last.bias.fill_(u_out_med)

        def rayleigh_now():
            with torch.enable_grad():
                xx = self._leaf_inputs(self.pts_np)
                ch = self.chi_with_env(xx)
                nn_ = max(float(torch.sqrt(torch.sum(
                    self.w * ch ** 2))), 1e-12)
                return float(pr.sobolev_energy(xx, ch / nn_, self.w))

        E_const = rayleigh_now()
        state_const = [t.detach().clone() for t in self._all_params()]
        opt = torch.optim.LBFGS(self._all_params(), lr=0.15,
                                max_iter=pre_steps, history_size=60,
                                line_search_fn='strong_wolfe')
        opt.step(pre_closure)
        E_pre = rayleigh_now()
        # the constant-u ansatz is a ground-state-like (nodeless) function:
        # only fall back to it for the ground state
        if n_state == 0 and (not np.isfinite(E_pre) or E_pre > E_const):
            with torch.no_grad():
                for t, b in zip(self._all_params(), state_const):
                    t.copy_(b)
        elif verbose:
            print(f'  pretrain: E {E_pre:.1f} (const {E_const:.1f})', flush=True)
        return None

    def _basis_torch(self, xs, mu, sig, n_state, scales=None):
        """Gaussian x Hermite basis as torch functions of xs (with graph)."""
        scales = self.basis_scales if scales is None else scales
        def hermite(t, k):
            # physicists' Hermite H_k via H_{n+1} = 2t H_n - 2n H_{n-1}
            h_prev, h_cur = torch.ones_like(t), 2.0 * t
            if k == 0:
                return h_prev
            for n in range(1, k):
                h_prev, h_cur = h_cur, 2.0 * t * h_cur - 2.0 * n * h_prev
            return h_cur
        cols = []
        G0 = None
        for x, m, s in zip(xs, mu, sig):
            g = torch.exp(-0.5 * ((x - m) / s) ** 2)
            G0 = g if G0 is None else G0 * g
        cols.append(G0.squeeze(-1) if G0.dim() > 1 else G0)
        k_max = n_state + self.k_extra
        pr = self.problem
        for sc in scales:
            Gt = None
            for x, m, s in zip(xs, mu, sig):
                g = torch.exp(-0.5 * ((x - m) / (s * sc)) ** 2)
                Gt = g if Gt is None else Gt * g
            Gs = Gt.squeeze(-1)
            if pr.dim == 1:
                t_i = (xs[0].squeeze(-1) - mu[0]) / (sig[0] * sc)
                for k in range(1, k_max + 1):
                    cols.append(Gs * hermite(t_i, k))
            elif pr.dim == 2:
                t0 = (xs[0].squeeze(-1) - mu[0]) / (sig[0] * sc)
                t1 = (xs[1].squeeze(-1) - mu[1]) / (sig[1] * sc)
                for k1 in range(k_max + 1):
                    for k2 in range(k_max + 1 - k1):
                        if k1 == 0 and k2 == 0:
                            continue
                        cols.append(Gs * hermite(t0, k1) * hermite(t1, k2))
            elif pr.dim == 3:
                k_3d = min(k_max, 5)
                t0 = (xs[0].squeeze(-1) - mu[0]) / (sig[0] * sc)
                t1 = (xs[1].squeeze(-1) - mu[1]) / (sig[1] * sc)
                t2 = (xs[2].squeeze(-1) - mu[2]) / (sig[2] * sc)
                for k1 in range(k_3d + 1):
                    for k2 in range(k_3d + 1 - k1):
                        for k3 in range(k_3d + 1 - k1 - k2):
                            if k1 == 0 and k2 == 0 and k3 == 0:
                                continue
                            cols.append(Gs * hermite(t0, k1) * hermite(t1, k2) * hermite(t2, k3))
            else:
                for i, x in enumerate(xs):
                    t_i = (x.squeeze(-1) - mu[i]) / (sig[i] * sc)
                    for k in range(1, k_max + 1):
                        cols.append(Gs * hermite(t_i, k))
        return cols

    # ------------------------------------------------------------------
    def train_state(self, n_state, steps=3000, lr0=1e-4, batch_size=None,
                    verbose_every=1000, seed=0, warm=True, lbfgs_steps=600,
                    lbfgs_chunk=40, lbfgs_lr=0.15, patience=3,
                    select_grid=None):
        """Adam refinement (small lr; the Sobolev loss is stiff) followed by
        full-batch LBFGS with a Wolfe line search (monotone descent)."""
        pr = self.problem
        rng = np.random.default_rng(seed + 917 * (n_state + 1))
        n_pool = len(self.w_np)
        prev_np = [s['psi'] for s in self.states]

        if warm:
            self.warm_start(n_state)
        use_batch = batch_size is not None and batch_size < n_pool

        def make_loss(xs, w, prev_idx=None):
            chi = self.chi_with_env(xs)
            psi = chi / self._norm_scalar()
            res = pr.h_psi(xs, psi) - self.E_param * psi
            # variance form: sum w |(H-E)psi|^2 (E trainable); zero iff psi is
            # an eigenfunction.  Hard normalisation removes the trivial sink.
            loss = torch.sum(w * res ** 2) / pr.E_scale ** 2
            if prev_np:
                if prev_idx is None:
                    pm_np = np.stack(prev_np)
                else:
                    pm_np = np.stack([p[prev_idx] for p in prev_np])
                pm = torch.tensor(pm_np, dtype=torch.float64,
                                  device=self.device)
                ov = torch.sum(w[None, :] * psi[None, :] * pm, dim=1)
                loss = loss + self.ortho_weight * torch.sum(ov * ov)
            return loss, self.E_param.detach()

        if steps:
            opt = torch.optim.Adam(self._all_params(), lr=lr0)
            sched = torch.optim.lr_scheduler.CosineAnnealingLR(
                opt, T_max=steps, eta_min=lr0 / 20)
            for it in range(steps):
                opt.zero_grad()
                if use_batch:
                    idx = rng.choice(n_pool, size=batch_size, replace=False)
                    xs = self._leaf_inputs([p[idx] for p in self.pts_np])
                    w = self.w[idx]
                else:
                    idx = None
                    xs = self._leaf_inputs(self.pts_np)
                    w = self.w
                loss, E = make_loss(xs, w, idx)
                loss.backward()
                opt.step()
                sched.step()
                if verbose_every and (it % verbose_every == 0 or it == steps - 1):
                    print(f'  state {n_state} [adam] step {it:5d}  E-Efloor = '
                          f'{(E.item()-pr.E_floor):12.2f} cm-1', flush=True)

        if lbfgs_steps and not use_batch:
            xs = self._leaf_inputs(self.pts_np)

            def closure():
                opt.zero_grad()
                loss, _ = make_loss(xs, self.w)
                loss.backward()
                return loss

            # chunked LBFGS with best-parameter tracking: near the optimum a
            # long LBFGS run can blow up through an ill-conditioned history,
            # so we keep the best Rayleigh quotient seen and stop on stalling.
            # seed the best-state tracker with the post-pretrain Ritz state:
            # chunked LBFGS may diverge, but can then never do worse than the
            # variational Ritz warm start
            def loss_now():
                with torch.enable_grad():
                    l, E = make_loss(xs, self.w)
                    return float(l), float(E)
            L0, E0_ray = loss_now()
            done, best_L, stall = 0, L0, 0
            if select_grid is not None:
                best_sel = self.evaluate_energies(
                    select_grid, state=len(self.states))['E_pde']
            else:
                best_sel = L0
            best_state = [t.detach().clone() for t in self._all_params()]
            if verbose_every:
                print(f'  state {n_state} [post-pretrain]  '
                      f'loss {L0:.5f}  E-Efloor = '
                      f'{E0_ray - pr.E_floor:12.3f} cm-1', flush=True)
            while done < lbfgs_steps and stall < patience:
                opt = torch.optim.LBFGS(self._all_params(), lr=lbfgs_lr,
                                        max_iter=min(lbfgs_chunk,
                                                     lbfgs_steps - done),
                                        history_size=60,
                                        line_search_fn='strong_wolfe')
                opt.step(closure)
                done += lbfgs_chunk
                L, E_ray = loss_now()
                if not np.isfinite(L) or not np.isfinite(E_ray):
                    with torch.no_grad():
                        for t, b in zip(self._all_params(), best_state):
                            t.copy_(b)
                    break
                with torch.no_grad():
                    self.E_param.fill_(E_ray)   # closed-form variance min
                if select_grid is not None:
                    E_sel = self.evaluate_energies(select_grid, state=len(
                        self.states))['E_pde']
                else:
                    E_sel = L
                if (E_sel if select_grid is not None else L) < best_sel - 1e-8:
                    best_sel = E_sel
                    best_L, stall = L, 0
                    best_state = [t.detach().clone() for t in self._all_params()]
                else:
                    stall += 1
                # deflate in function space: project the current function
                # orthogonal to the converged lower states and refit the net
                # to the projection, so the next chunk optimises inside the
                # admissible (orthogonal) subspace
                if prev_np and done < lbfgs_steps:
                    with torch.enable_grad():
                        ch = self.chi_with_env(xs)
                        nn_ = max(float(torch.sqrt(torch.sum(
                            self.w * ch ** 2))), 1e-12)
                        chn = (ch / nn_).detach().numpy()
                    for pm in prev_np:
                        ovp = float(np.sum(w_sub_full := self.w_np * chn * pm))
                        chn = chn - ovp * pm
                    nrm_p = float(np.sqrt(np.sum(self.w_np * chn ** 2)))
                    if np.isfinite(nrm_p) and nrm_p > 1e-8:
                        tgt = torch.tensor(chn / nrm_p, dtype=torch.float64,
                                           device=self.device)
                        txs = self._leaf_inputs(self.pts_np)
                        lw_t = self.w

                        def defl_closure():
                            opt2.zero_grad()
                            c = self.chi_with_env(txs)
                            l = torch.sum(lw_t * (c - tgt) ** 2)
                            l.backward()
                            return l
                        opt2 = torch.optim.LBFGS(self._all_params(), lr=0.15,
                                                max_iter=200, history_size=60,
                                                line_search_fn='strong_wolfe')
                        opt2.step(defl_closure)
                if verbose_every:
                    print(f'  state {n_state} [lbfgs] ~{done:5d} it  '
                          f'loss {L:.6f}  E-Efloor = '
                          f'{E_ray - pr.E_floor:12.3f} cm-1', flush=True)
            if best_state is not None:
                with torch.no_grad():
                    for t, b in zip(self._all_params(), best_state):
                        t.copy_(b)
                # closed-form E refresh for the restored best state
                with torch.enable_grad():
                    chi_b = self.chi_with_env(xs)
                    nrm_b = max(float(torch.sqrt(torch.sum(
                        self.w * chi_b ** 2))), 1e-12)
                    E_b = float(pr.sobolev_energy(xs, chi_b / nrm_b, self.w))
                    with torch.no_grad():
                        self.E_param.fill_(E_b)

        return self.finalize_state(n_state)

    def _norm_scalar(self):
        with torch.no_grad():
            xs = self._leaf_inputs(self.pts_np)
            return max(float(torch.sqrt(torch.sum(
                self.w * self.chi_with_env(xs) ** 2))), 1e-12)

    def finalize_state(self, n_state):
        """Full-pool energies (Sobolev and PDE forms) and stored psi.

        The normalisation constant MUST be detached to a float: keeping it as
        a graph tensor makes d(psi)/dx pick up a spurious d(norm)/dx term and
        corrupts the kinetic energy (values are unaffected)."""
        pr = self.problem
        xs = self._leaf_inputs(self.pts_np)
        chi = self.chi_with_env(xs)
        with torch.no_grad():
            nrm = float(torch.sqrt(torch.sum(self.w * chi ** 2)))
        psi = chi / nrm
        E_sob = float(pr.sobolev_energy(xs, psi, self.w))
        hpsi = pr.h_psi(xs, psi)
        E_pde = float(torch.sum(self.w * psi * hpsi))
        psi_d = psi.detach().cpu().numpy()
        j = int(np.argmax(np.abs(psi_d)))
        if psi_d[j] < 0:
            psi_d = -psi_d
        self.states.append(dict(psi=psi_d, E=E_sob, E_pde=E_pde))
        print(f'  state {n_state}: E(Sobolev) = {E_sob:.2f}, '
              f'E(PDE form) = {E_pde:.2f} cm-1', flush=True)
        return E_sob, psi_d

    def evaluate_energies(self, n_per_dim, state=-1):
        """Independent-grid evaluation: rebuild a Gauss-Legendre quadrature
        of the requested size (finer than the training pool) and report both
        the Sobolev and PDE-form energies.  Guards against the training
        quadrature being too coarse for the converged wavefunction."""
        pr = self.problem
        pts_np, w_np = pr.quadrature(n_per_dim)
        xs = [torch.tensor(np.asarray(p, float), dtype=torch.float64,
                           device=self.device, requires_grad=True).unsqueeze(-1)
              for p in pts_np]
        chi = self.chi_with_env(xs)
        nrm = max(float(torch.sqrt(torch.sum(
            torch.tensor(w_np, dtype=torch.float64, device=self.device)
            * chi ** 2))), 1e-12)
        psi = chi / nrm
        w = torch.tensor(w_np, dtype=torch.float64, device=self.device)
        E_sob = float(pr.sobolev_energy(xs, psi, w))
        E_pde = float(torch.sum(w * psi * pr.h_psi(xs, psi)))
        return dict(E_sob=E_sob, E_pde=E_pde)

    # ------------------------------------------------------------------
    def evaluate(self, pts_list, state=-1):
        nrm = np.sqrt(np.sum(self.w_np * self.states[state]['psi'] ** 2))
        xs = []
        for p, lo, hi in zip(pts_list, self.lo, self.hi):
            s = torch.tensor(np.asarray(p, float), dtype=torch.float64,
                             device=self.device).unsqueeze(-1)
            xs.append(s)
        with torch.no_grad():
            chi = self.chi_with_env(xs).cpu().numpy()
        return chi / nrm
