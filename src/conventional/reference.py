"""Classical reference solvers used to validate the PINN results.

* fd_bound_1d      : finite-difference bound states in 1D
* scatter_fd_1d    : fixed-energy 1D scattering state.  Complex tridiagonal
                     solve with a Robin outgoing-wave condition at the right
                     boundary (incoming amplitude normalised to 1); the
                     reflection amplitude R is extracted by a WKB-phase
                     plane-wave projection on the asymptotic plateau.
* LegendreAngularDVR : Gauss-Legendre DVR in x = cos(gamma) with the exact j^2
                     matrix (eigenvalues j(j+1) on Legendre polynomials).
* dvr_2d, dvr_3d   : product FD/DVR bound-state solvers for the water
                     Hamiltonian in Jacobi coordinates (J = 0, reduced
                     wavefunction chi = r R psi, plain L^2 measure).
"""
import numpy as np
import scipy.linalg as sla
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.special import eval_legendre, roots_legendre

try:                                     # NumPy >= 2.0 / < 2.0
    from src.common.compat import trapezoid
except ImportError:                      # run as a standalone script
    trapezoid = getattr(np, 'trapezoid', np.trapz)


# --------------------------------------------------------------------------
# 1D reference solvers
# --------------------------------------------------------------------------
def fd_bound_1d(V, xlo, xhi, n_states, npts=8000, fd_order=4):
    """Lowest bound states of -B psi'' + V psi = E psi on [xlo, xhi]."""
    x = np.linspace(xlo, xhi, npts)
    h = x[1] - x[0]
    B = V.B
    n = npts
    if fd_order == 2:
        main = 2.0 * B / h ** 2 + V.V(x)
        off1 = -B / h ** 2 * np.ones(n - 1)
        H = sp.diags([off1, main, off1], [-1, 0, 1], format='csc')
    else:
        c2, c1, c0 = 1.0, -16.0, 30.0
        f = B / (12.0 * h ** 2)
        main = c0 * f + V.V(x)
        off1 = c1 * f * np.ones(n - 1)
        off2 = c2 * f * np.ones(n - 2)
        H = sp.diags([off2, off1, main, off1, off2], [-2, -1, 0, 1, 2],
                     format='csc')
    vals, vecs = spla.eigsh(H, k=n_states, sigma=float(V.V(x).min()) - 1.0,
                            which='LM')
    order = np.argsort(vals)
    vals, vecs = vals[order], vecs[:, order]
    for i in range(vecs.shape[1]):                      # fix sign at max |psi|
        j = np.argmax(np.abs(vecs[:, i]))
        if vecs[j, i] < 0:
            vecs[:, i] *= -1
    nrm = np.sqrt(trapezoid(vecs ** 2, x, axis=0))
    return vals, vecs / nrm, x


def scatter_fd_1d(V, E, xlo, xhi, npts=16000, fit_length=None):
    """Fixed-energy (E > 0) scattering state of the Morse potential.

    Incoming wave from the right (open side):  psi ~ e^{+ikx} + R e^{-ikx}
    for x -> +inf; psi -> 0 into the repulsive wall (x -> -inf).

    Method: complex tridiagonal FD solve.  The last grid equation is replaced
    by the 2nd-order backward-difference Robin condition
        psi'(x_R) + i k psi(x_R) = 2 i k e^{i k x_R}
    which enforces unit incoming amplitude.  R is then obtained by projecting
    the solution on WKB waves k(x)^{-1/2} e^{+-i theta(x)} over an asymptotic
    window where |V| < 0.01 E (theta = integral of the local wavenumber).

    Returns x, psi (complex, incoming amp = 1), R (complex), k.
    """
    x = np.linspace(xlo, xhi, npts)
    h = x[1] - x[0]
    n = npts
    B = V.B
    k = float(V.k_continuum(E))
    pot = V.V(x)
    diag = 2.0 * B / h ** 2 + pot - E
    off = -B / h ** 2 * np.ones(n - 1)
    A = sp.diags([off, diag, off], [-1, 0, 1], shape=(n, n), format='lil', dtype=complex)
    b = np.zeros(n, dtype=complex)

    # left Dirichlet psi_0 = 0
    A[0, :] = 0.0
    A[0, 0] = 1.0
    b[0] = 0.0

    # right Robin: (3 psi_{n-1} - 4 psi_{n-2} + psi_{n-3})/(2h) + ik psi = 2ik e^{ikx}
    j = n - 1
    A[j, :] = 0.0
    A[j, j] = 3.0 / (2.0 * h) + 1j * k
    A[j, j - 1] = -4.0 / (2.0 * h)
    A[j, j - 2] = 1.0 / (2.0 * h)
    b[j] = 2j * k * np.exp(1j * k * x[j])

    psi = spla.spsolve(A.tocsc(), b)

    # Plane-wave projection for R on the asymptotic plateau (where |V| << E)
    if fit_length is None:
        fit_length = 0.8
    mask = (x > max(xhi - fit_length, xlo + 2.0)) & (x < xhi - 2 * h)
    M = np.column_stack([np.exp(1j * k * x[mask]), np.exp(-1j * k * x[mask])])
    (c_in, c_out), *_ = np.linalg.lstsq(M, psi[mask], rcond=None)

    # Rescale to unit incoming plane-wave amplitude
    psi = psi / c_in
    R = c_out / c_in
    return x, psi, R, k


# --------------------------------------------------------------------------
# quadrature / angular DVR
# --------------------------------------------------------------------------
def gauss_legendre_1d(a, b, n):
    """Gauss-Legendre nodes and weights on [a, b]."""
    t, w = roots_legendre(n)
    return 0.5 * (b - a) * t + 0.5 * (a + b), 0.5 * (b - a) * w


class LegendreAngularDVR:
    """Gauss-Legendre DVR in x = cos(gamma) with the exact j^2 matrix.

    j^2 = -[d2/dgamma^2 + cot(gamma) d/dgamma] (units of hbar^2 dropped) has
    eigenvalues j(j+1) on P_j(x).  DVR matrix: M = T diag(j(j+1)) T^T with
    T_{b j} = sqrt(w_b) sqrt((2j+1)/2) P_j(x_b)  (orthogonal for j <= n-1).
    """

    def __init__(self, n_gamma):
        self.n = n_gamma
        self.x, self.w = roots_legendre(n_gamma)
        j = np.arange(n_gamma)
        P = np.array([eval_legendre(jj, self.x) for jj in j])
        T = np.sqrt(self.w)[None, :] * (P * np.sqrt((2 * j + 1) / 2.0)[:, None])
        self.T = T
        M = (T.T * (j * (j + 1))[None, :]) @ T          # T^T diag(j(j+1)) T
        self.j2 = 0.5 * (M + M.T)


def _fd_second(n, h, order=4):
    """Matrix approximating +d^2/dx^2 (Dirichlet ghost = 0), 2- or 5-point."""
    if order == 2:
        return sp.diags([np.ones(n - 1), -2.0 * np.ones(n), np.ones(n - 1)],
                        [-1, 0, 1], format='csr') / h ** 2
    f = 1.0 / (12.0 * h ** 2)
    return sp.diags([-np.ones(n - 2), 16.0 * np.ones(n - 1),
                     -30.0 * np.ones(n), 16.0 * np.ones(n - 1),
                     -np.ones(n - 2)], [-2, -1, 0, 1, 2], format='csr') * f


# --------------------------------------------------------------------------
# water reference solvers (Jacobi coordinates, J = 0, chi = r R psi)
# --------------------------------------------------------------------------
def dvr_2d(pes, r_fix, R_lo, R_hi, n_R, n_gamma, k=10, fd_order=4):
    """2D frozen-r water problem  H chi = E chi  on (R, x = cos gamma):

    H = -BR d2_R + (BR/R^2 + Br/r_fix^2) j^2 + V(r_fix, R, gamma).
    Norm measure: dR dx.  Returns (E, vecs, R, x, w, Vgrid).
    """
    R = np.linspace(R_lo, R_hi, n_R)
    hR = R[1] - R[0]
    ang = LegendreAngularDVR(n_gamma)
    xg, wg = ang.x, ang.w
    gamma = np.arccos(xg)
    RR, GG = np.meshgrid(R, gamma, indexing='ij')
    Vg = pes.V_jacobi(r_fix, RR, GG)
    C = pes.BR / R ** 2 + pes.Br / r_fix ** 2

    H = sp.kron(-pes.BR * _fd_second(n_R, hR, fd_order),
                sp.identity(n_gamma), format='csr')
    H += sp.kron(sp.diags(C), sp.csr_matrix(ang.j2), format='csr')
    H += sp.diags(Vg.ravel())

    vals, vecs = spla.eigsh(H, k=k, sigma=float(Vg.min()) - 1.0, which='LM')
    order = np.argsort(vals)
    return vals[order], vecs[:, order], R, xg, wg, Vg


def dvr_3d(pes, r_lo, r_hi, n_r, R_lo, R_hi, n_R, n_gamma, k=8,
           fd_order=4, tol=0.0, verbose=False):
    """3D water bound states, matrix assembled sparsely.

    H = -Br d2_r - BR d2_R + (Br/r^2 + BR/R^2) j^2 + V(r, R, gamma)
    acting on chi(r, R, x = cos gamma); plain L^2(dr dx) measure.
    Returns (E, vecs, (r, R, x, w), Vgrid).
    """
    r = np.linspace(r_lo, r_hi, n_r)
    R = np.linspace(R_lo, R_hi, n_R)
    hr, hR = r[1] - r[0], R[1] - R[0]
    ang = LegendreAngularDVR(n_gamma)
    xg, wg = ang.x, ang.w
    gamma = np.arccos(xg)
    rr, RR, GG = np.meshgrid(r, R, gamma, indexing='ij')
    Vg = pes.V_jacobi(rr, RR, GG)
    C = (pes.Br / r ** 2)[:, None] + (pes.BR / R ** 2)[None, :]   # (n_r, n_R)
    n_r, n_R, n_g = len(r), len(R), n_gamma

    Ir = sp.identity(n_r, format='csr')
    IG = sp.identity(n_g, format='csr')
    H = sp.kron(-pes.Br * _fd_second(n_r, hr, fd_order),
                sp.identity(n_R * n_g), format='csr')
    H += sp.kron(Ir, sp.kron(-pes.BR * _fd_second(n_R, hR, fd_order), IG),
                 format='csr')
    H += sp.kron(sp.diags(C.ravel()), sp.csr_matrix(ang.j2), format='csr')
    H += sp.diags(Vg.ravel())

    k_eff = min(k + 4, n_r * n_R * n_g - 2)
    vals, vecs = spla.eigsh(H, k=k_eff, which='SA', tol=tol)
    order = np.argsort(vals)
    vals, vecs = vals[order][:k], vecs[:, order][:, :k]
    if verbose:
        print('  dvr_3d:', np.round(vals, 3))
    return vals, vecs, (r, R, xg, wg), Vg
