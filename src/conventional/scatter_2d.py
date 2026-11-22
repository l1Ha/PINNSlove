"""A+BC (H + OH rigid-bond) fixed-energy scattering on the water PES.

KNOWN LIMITATION (honest negative result, 2026-10-11 audit): the model water
PES has a deep well (V_min = -2 D_e = -81 140 cm^-1) and the scattering
energy E_tot = E_c - D_e sits amid a dense manifold of well levels.  The
finite-box BVP solution therefore varies strongly with the box size and grid
(the interior solution changes by O(1) when R_hi moves 9 -> 15) because the
crude local Robin boundary cannot absorb the resonant multichannel flux.
A converged reference requires CAP absorbing potentials or log-derivative
matching.  Until such a reference is implemented, PINN-vs-reference accuracy
claims for this 2D model are NOT meaningful and are not made.

2D model in (R, gamma) with the reduced wavefunction chi = R psi:
  H = -BR d2_R + (BR/R^2 + Br/re^2) j^2 + V(re, R, gamma),
  E_tot = E_c + V(re, R->inf)  (OH at its Morse minimum).
Scattered-wave BVP (see assemble_2d_scatter): chi = e^{+iKR} + chi_sc with
  (H - E_tot) chi_sc = -(H - E_tot) e^{+iKR},  chi_sc(R_lo) = 0,
  outgoing Robin or Dirichlet-under-CAP at R_hi.  scatter_fd_2d returns
  chi_sc; the full wave is e^{+iKR} + chi_sc.
"""
import numpy as np
import torch
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import scipy.linalg as sla

from src.conventional.reference import LegendreAngularDVR, _fd_second


def assemble_2d_scatter(pes, E_c, R_lo=0.6, R_hi=9.0, n_R=400, n_gamma=60,
                        fd_order=4, cap_start=None, cap_eta=None,
                        cap_width=2.0):
    """Complex sparse BVP system for the SCATTERED wave; returns
    (A, b, R grid, x=cos(gamma) grid).

    Formulation (time convention e^{+iEt}): the full wave is
        chi = chi_in + chi_sc,   chi_in = e^{+iKR}   (unit isotropic
        incoming VALUE amplitude, moving toward -R).
    chi_sc solves the INHOMOGENEOUS problem
        (H - E_tot) chi_sc = -(H - E_tot) chi_in      (source = interaction
        region only; (H - E) chi_in = (E_c + V - E_tot) e^{+iKR} vanishes in
        the flat asymptote V = -D_e),
        chi_sc(R_lo) = 0,
        R = R_hi:  channel-projected OUTGOING Robin
            (d_R - i k_ch) [chi_sc]_ch = 0   (no CAP),   or
            chi_sc(R_hi) = 0                  (CAP active).
    With a quadratic CAP absorber -i*eta*s^2 (s = max(0, (R-cap_start)/
    cap_width), fixed cap_width independent of R_hi) added to H, the outgoing
    flux is absorbed inside the box and the deep-well box resonance documented
    below is removed.  cap_eta is chosen from the plateau of the interior
    solution under box/grid change (run_cap_scan.py), a standard
    regularization study, not a fitted parameter.

    KNOWN LIMITATION (honest negative result, 2026-10-11 audit): WITHOUT the
    CAP the finite-box solution is NOT grid/box convergent (deep well
    V_min = -2 D_e = -81 140 cm^-1, dense manifold of well levels; verified
    quantitatively by scripts/conventional/run_scatter2d_convergence.py:
    interior L2 differences 128-156% between grids, 83-136% between box
    sizes).  The CAP path exists precisely to cure this.
    """
    R = np.linspace(R_lo, R_hi, n_R)
    h = R[1] - R[0]
    ang = LegendreAngularDVR(n_gamma)
    xg = ang.x
    gamma = np.arccos(xg)
    RR, GG = np.meshgrid(R, gamma, indexing='ij')
    Vg = pes.V_jacobi(pes.re, RR, GG)
    # H + OH(v=0) asymptote: one OH Morse at its minimum, the other bond broken
    E_tot = E_c + float(pes.V_valence(pes.re, pes.re, pes.theta_e)) + pes.De
    K = float(np.sqrt(E_c / pes.BR))

    H = sp.kron(-pes.BR * _fd_second(n_R, h, fd_order),
                sp.identity(n_gamma), format='csr')
    C = pes.BR / R ** 2 + pes.Br / pes.re ** 2
    H += sp.kron(sp.diags(C), sp.csr_matrix(ang.j2), format='csr')
    H += sp.diags(Vg.ravel())
    # RHS: source that drives the scattered wave.  chi_in = e^{+iKR} is
    # gamma-independent, so the angular j^2 term vanishes on it exactly and
    # the discrete (H - E) chi_in needs only the matvec below.
    chi_in = np.exp(1j * K * RR)
    b = -(H @ chi_in.ravel() - E_tot * chi_in.ravel())
    # CAP absorber on the physical H (after the source is formed from the
    # physical H, so the incoming wave is not attenuated)
    if cap_start is not None and cap_eta:
        s = np.clip((R - cap_start) / cap_width, 0.0, None)
        cap_R = -1j * cap_eta * s ** 2
        H = H + sp.kron(sp.diags(cap_R), sp.identity(n_gamma), format='csr')
    A = (H - E_tot * sp.identity(H.shape[0])).astype(complex).tolil()

    # Dirichlet at R_min: chi_sc = 0
    j0 = 0
    A[j0, :] = 0.0
    A[j0, j0] = 1.0
    b[j0] = 0.0

    # Outer boundary at R_max.
    B_rot = pes.Br / pes.re ** 2
    evals, evecs = np.linalg.eigh(np.asarray(ang.j2))
    coef = np.array([3.0, -4.0, 1.0]) / (2.0 * h)
    j = n_R - 1
    if cap_start is not None and cap_eta:
        # CAP has absorbed the outgoing wave; plain Dirichlet.  Any spurious
        # reflection from this wall travels back through the absorber.
        for ch in range(n_gamma):
            row = j * n_gamma + ch
            A[row, :] = 0.0
            A[row, row] = 1.0
            b[row] = 0.0
    else:
        # Channel-projected OUTGOING Robin: (d_R - i k_ch) [chi_sc]_ch = 0.
        # Time e^{+iEt}: outgoing = e^{-iKR}; k_ch = sqrt((E_c - B_rot j(j+1))
        # / BR) is real for open channels, imaginary (decaying) for closed.
        for ch in range(n_gamma):
            eval_ch = max(evals[ch], 0.0)
            E_ch = E_c - B_rot * eval_ch
            k_ch = np.sqrt(E_ch / pes.BR) if E_ch > 0 else 1j * np.sqrt(-E_ch / pes.BR)
            row = j * n_gamma + ch
            A[row, :] = 0.0
            for t, cv in enumerate(coef):
                A[row, (j - t) * n_gamma:(j - t + 1) * n_gamma] += cv * evecs[:, ch]
            A[row, j * n_gamma:(j + 1) * n_gamma] += -1j * k_ch * evecs[:, ch]
            b[row] = 0.0
    return A.tocsc(), b, R, xg, K, E_tot


def scatter_fd_2d(pes, E_c, **kw):
    A, b, R, xg, K, E_tot = assemble_2d_scatter(pes, E_c, **kw)
    chi = spla.spsolve(A, b)
    n_gamma = len(xg)
    return chi.reshape(len(R), n_gamma), R, xg, K, E_tot


