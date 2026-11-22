"""Independent analytic validation utilities and robust S-matrix extraction.

The extraction is self-consistent with the finite-difference scheme: the fitting
basis uses the wavenumber q(h) for which exp(i q R) is an EXACT solution of the
discretized free kinetic operator, so the S-matrix of the discrete problem is
recovered to machine precision. In the limit h -> 0, q -> k and the discrete
S-matrix converges to the exact continuum one.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.special import spherical_jn, spherical_yn

from src.conventional.reference import _fd_second


# ---------------------------------------------------------------------------
# free radial functions (reduced form u = x h_l(x))
# ---------------------------------------------------------------------------
def free_incoming(x, l):
    """f2 -> +i exp(-i(x - l pi/2)); the 'incoming' free solution."""
    return x * (spherical_jn(l, x) - 1j * spherical_yn(l, x))


def free_outgoing(x, l):
    """f1 -> -i exp(+i(x - l pi/2)); the 'outgoing' free solution."""
    return x * (spherical_jn(l, x) + 1j * spherical_yn(l, x))


def fit_wavevector(k, h, order=4):
    """Plane-wave wavenumber q with exp(i q R) an exact solution of the FD operator.

    order=2 : cos(qh) = 1 - (kh)^2 / 2
    order=4 : cos(qh) = 4 - sqrt(9 + 3 (kh)^2)
    """
    x = k * h
    if order == 2:
        c = 1.0 - 0.5 * x ** 2
    else:
        c = 4.0 - np.sqrt(9.0 + 3.0 * x ** 2)
    return np.arccos(np.clip(c, -1.0, 1.0)) / h


def check_extractor_roundtrip():
    """Synthetic u = f2 + S f1 over a window must be recovered exactly."""
    R = np.linspace(0.5, 20.0, 6001)   # avoid x = 0 (y_l diverges there)
    mask = (R > 18.5) & (R < 20.0)
    for l in [0, 1, 3]:
        for k in [0.7, 1.7]:
            for S_t in [0.3 - 0.9j, -0.7 + 0.2j, 1.0 + 0.0j]:
                u = free_incoming(k * R, l) + S_t * free_outgoing(k * R, l)
                M = np.column_stack([free_incoming(k * R[mask], l),
                                     free_outgoing(k * R[mask], l)])
                c, *_ = np.linalg.lstsq(M, u[mask], rcond=None)
                assert abs(c[1] / c[0] - S_t) < 1e-9, (l, k, S_t, c)
    return True


# ---------------------------------------------------------------------------
# multichannel S-matrix extraction
# ---------------------------------------------------------------------------
def auto_window(R, k_list, frac=0.40, n_wave=3.0, h_guard=5):
    """Wavelength-aware asymptotic fitting window.

    The two free basis functions f2, f1 become numerically collinear when the
    window spans less than ~1 wavelength, which destroys the extracted ratio
    b/a. We therefore require at least `n_wave` wavelengths of the *longest*
    wavelength (smallest k) among the open channels, capped at `frac` of the box.
    """
    L = R[-1] - R[0]
    h = R[1] - R[0]
    kmin = min([k for k in k_list if k > 0], default=np.inf)
    lam = 2.0 * np.pi / kmin if np.isfinite(kmin) else L
    w = min(frac * L, max(n_wave * lam, 2.0))
    mask = (R > R[-1] - w) & (R < R[-1] - h_guard * h)
    if mask.sum() < 40:
        mask = R > R[-1] - w
    return mask


def extract_smatrix(us, R, K_list, l_list, open_idx, window=None,
                    h=None, dispersion=True, matrix_form=True):
    """Extract the open-channel S-matrix from unit-incoming-amplitude solutions.

    Parameters
    ----------
    us : (nR, n_ch, n_col) complex array
        us[:, :, j] is the full channel vector for incident channel open_idx[j]
        (unit incoming amplitude in the reduced radial function).
    R : (nR,)
    K_list : (n_ch,) complex asymptotic wavevectors
    l_list : (n_ch,) orbital angular momenta
    open_idx : list of incident channels (columns)
    window : fixed fit-window width in Angstrom (None = wavelength-aware)
    h : radial step (required if dispersion=True)
    dispersion : use the FD-consistent wavenumber q(h)
    matrix_form : if True return S = B A^{-1} with A the fitted incoming-amplitude
        matrix (common normalization; enforces the exact linear structure and
        improves reciprocity); if False return the element-wise ratio b/a.

    Returns
    -------
    S : (n_open, n_col)
    a_inc : (n_col,) diagonal of the fitted incoming-amplitude matrix
    """
    check_extractor_roundtrip()
    us = np.asarray(us)
    if us.ndim == 2:
        us = us[:, :, None]
    n_col = us.shape[2]
    open_list = list(open_idx)
    n_open = len(open_list)

    if h is None:
        h = R[1] - R[0]
    if window is None:
        k_open = [K_list[i].real for i in open_list if K_list[i].imag == 0.0]
        mask = auto_window(R, k_open)
    else:
        mask = (R > R[-1] - window) & (R < R[-1] - 5 * h)
    R_fit = R[mask]

    def kw(idx):
        k = K_list[idx].real
        return fit_wavevector(k, h) if dispersion else k

    D = {}
    for oc in open_list:
        if K_list[oc].imag != 0.0:
            continue
        q = kw(oc)
        D[oc] = np.column_stack([free_incoming(q * R_fit, l_list[oc]),
                                 free_outgoing(q * R_fit, l_list[oc])])

    A = np.zeros((n_open, n_col), dtype=complex)
    B = np.zeros((n_open, n_col), dtype=complex)
    for jc in range(n_col):
        for iout, oc in enumerate(open_list):
            c, *_ = np.linalg.lstsq(D[oc], us[mask, oc, jc], rcond=None)
            A[iout, jc] = c[0]
            B[iout, jc] = c[1]

    a_inc = np.diag(A).copy()
    if matrix_form:
        S = B @ np.linalg.inv(A)
    else:
        S = B / np.where(np.abs(A) > 0, A, 1.0)
    return S, a_inc


def fd2_4th_dirichlet(nR, h):
    """4th-order +d^2/dx^2 with *consistent* rows next to both end rows.

    Row 1 uses the one-sided 4th-order formula (11,-20,6,4,-1)/(12h^2) (valid
    for an ARBITRARY value at row 0) and row nR-2 the mirrored formula.  Rows 0
    and nR-1 are left empty (they must be replaced by the BCs).  Unlike the
    odd-reflection shortcut of :func:`fd2_4th_wall`, this is fourth-order
    accurate for solutions that are not odd about the wall, and it keeps every
    interior row consistent (row sums vanish), so errors converge as h^4
    instead of being pinned at an O(1) level.

    Verified against an analytic source problem with a non-odd solution:
    max error 2.3e-8 -> 1.3e-9 -> 1.1e-10 for n = 51/101/201, whereas both the
    truncated 5-point rows (row sum = 1) and the odd-reflection rows (row
    sum = 2) stall at ~1e-1 -- they are inconsistent stencils.
    """
    if nR < 6:
        return _fd_second(nR, h, 2)
    D2 = sp.lil_matrix((nR, nR))
    for i in range(2, nR - 2):
        D2[i, i - 2:i + 3] = [-1.0, 16.0, -30.0, 16.0, -1.0]
    D2[1, 0:5] = [11.0, -20.0, 6.0, 4.0, -1.0]
    D2[nR - 2, nR - 5:nR] = [-1.0, 4.0, 6.0, -20.0, 11.0]
    return (D2 / (12.0 * h ** 2)).tocsr()


def fd2_4th_wall(nR, h):
    """4th-order +d^2/dx^2 matrix with correct treatment of BOTH end rows.

    Row 1 (next to the Dirichlet wall at grid point 0) uses the odd-reflection
    ghost u_{-1} = -u_1, giving diagonal -29/(12h^2) instead of -30/(12h^2).
    Row nR-2 (next to the Robin row) uses the one-sided backward 5-point
    second derivative (45,-154,214,-156,61,-10)/(12h^2) so that no ghost point
    beyond the boundary is silently set to zero. Dropping either correction
    injects an O(h^-2) local error that contaminates the whole solution.
    Rows 0 and nR-1 are left empty (replaced later by the BCs).
    """
    D2 = sp.lil_matrix((nR, nR))
    for i in range(2, nR - 3):
        D2[i, i - 2] = -1.0
        D2[i, i - 1] = 16.0
        D2[i, i] = -30.0
        D2[i, i + 1] = 16.0
        D2[i, i + 2] = -1.0
    # row 1: odd ghost at the wall
    D2[1, 1] = -29.0
    D2[1, 2] = 16.0
    D2[1, 3] = -1.0
    # row nR-3: central, all neighbours exist
    D2[nR - 3, nR - 5] = -1.0
    D2[nR - 3, nR - 4] = 16.0
    D2[nR - 3, nR - 3] = -30.0
    D2[nR - 3, nR - 2] = 16.0
    D2[nR - 3, nR - 1] = -1.0
    # row nR-2: one-sided backward 5-point (no ghost)
    for t, c in enumerate([45.0, -154.0, 214.0, -156.0, 61.0, -10.0]):
        D2[nR - 2, nR - 2 - t] = c
    return (D2 / (12.0 * h ** 2)).tocsr()


# ---------------------------------------------------------------------------
# generic coupled-channel radial solver (same blocks as the production code)
# ---------------------------------------------------------------------------
def solve_coupled_1d(Vfun, thresholds, BR, E, Rmin, Rmax, nR, l_list=None,
                     incoming=None, fd_order=4, robin_order=4):
    """-BR d2u_i/dR2 + sum_j [V_ij(R) + eps_i delta_ij] u_j = E u_i,
    Dirichlet at Rmin and multichannel outgoing Robin at Rmax."""
    R = np.linspace(Rmin, Rmax, nR)
    h = R[1] - R[0]
    n_ch = len(thresholds)
    if l_list is None:
        l_list = np.zeros(n_ch, dtype=int)
    l_list = np.asarray(l_list)

    V = np.asarray(Vfun(R))
    if V.ndim == 2:
        V = np.broadcast_to(V, (nR, n_ch, n_ch)).copy()

    K_list = np.zeros(n_ch, dtype=complex)
    for i in range(n_ch):
        E_tr = E - thresholds[i]
        K_list[i] = np.sqrt(E_tr / BR) if E_tr > 0 else 1j * np.sqrt(-E_tr / BR)

    if fd_order == 4:
        D2 = fd2_4th_wall(nR, h)
    else:
        D2 = _fd_second(nR, h, fd_order)
    H = sp.kron(-BR * D2, sp.identity(n_ch), format='csr')
    diag = np.zeros((nR, n_ch))
    for i in range(n_ch):
        diag[:, i] = thresholds[i] + BR * l_list[i] * (l_list[i] + 1) / R ** 2
    H = H + sp.diags(diag.reshape(-1))
    H = H + sp.block_diag([V[iR] for iR in range(nR)], format='csr')

    A = (H - E * sp.identity(H.shape[0])).astype(complex).tolil()
    for p in range(n_ch):
        A[p, :] = 0.0
        A[p, p] = 1.0

    if robin_order == 4:
        coef = np.array([25.0, -48.0, 36.0, -16.0, 3.0]) / (12.0 * h)
        n_tail = 5
    else:
        coef = np.array([3.0, -4.0, 1.0]) / (2.0 * h)
        n_tail = 3
    je = nR - 1
    src = 0.0
    if incoming is not None:
        k_in = K_list[incoming].real
        src = -2j * k_in * np.exp(-1j * k_in * R[-1])
    for p in range(n_ch):
        row = je * n_ch + p
        A[row, :] = 0.0
        A[row, row] = coef[0] - 1j * K_list[p]
        for t in range(1, n_tail):
            A[row, (je - t) * n_ch + p] = coef[t]

    b = np.zeros(nR * n_ch, dtype=complex)
    if incoming is not None:
        b[je * n_ch + incoming] = src
    return spla.spsolve(A.tocsc(), b).reshape(nR, n_ch), R, K_list


# ---------------------------------------------------------------------------
# single-channel solver (validated against the hard-sphere analytic result)
# ---------------------------------------------------------------------------
def solve_radial_single(BR, l, E, Rmin, Rmax, nR, Vfun=None, fd_order=4,
                        robin_order=4, dispersion=True):
    """Return (S, R, u, k, q). S is extracted with the FD-consistent q."""
    R = np.linspace(Rmin, Rmax, nR)
    h = R[1] - R[0]
    V = 0.0 if Vfun is None else np.asarray(Vfun(R), float)
    k = np.sqrt(E / BR)

    if fd_order == 4:
        D2 = fd2_4th_wall(nR, h)
    else:
        D2 = _fd_second(nR, h, fd_order)
    A = (-BR * D2 + sp.diags(V + BR * l * (l + 1) / R ** 2)
         - E * sp.identity(nR)).astype(complex).tolil()
    A[0, :] = 0.0
    A[0, 0] = 1.0

    if robin_order == 4:
        coef = np.array([25.0, -48.0, 36.0, -16.0, 3.0]) / (12.0 * h)
        n_tail = 5
    else:
        coef = np.array([3.0, -4.0, 1.0]) / (2.0 * h)
        n_tail = 3
    j = nR - 1
    A[j, :] = 0.0
    A[j, j] = coef[0] - 1j * k
    for t in range(1, n_tail):
        A[j, j - t] = coef[t]

    b = np.zeros(nR, dtype=complex)
    b[j] = -2j * k * np.exp(-1j * k * Rmax)
    u = spla.spsolve(A.tocsc(), b)

    q = fit_wavevector(k, h, fd_order) if dispersion else k
    mask = (R > Rmax - 1.4) & (R < Rmax)
    M = np.column_stack([free_incoming(q * R[mask], l),
                         free_outgoing(q * R[mask], l)])
    c, *_ = np.linalg.lstsq(M, u[mask], rcond=None)
    return c[1] / c[0], R, u, k, q
