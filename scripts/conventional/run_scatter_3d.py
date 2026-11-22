"""Full-dimensional (3D: r, R, gamma) reactive scattering for H + H2(v=0) on the authentic BKMP2 PES.

Coordinates (pairwise Jacobi coordinates for J=0):
  r: H2 bond length (vibrational coordinate)
  R: H - H2 collision distance (translational scattering coordinate)
  gamma: angle between r and R (rotational coupling coordinate)

Solves:
  [-Br d2_r - BR d2_R + (Br/r^2 + BR/R^2) j^2 + V_BKMP2(r, R, gamma) - E_tot] chi = 0
with exact unconstrained Kronecker FD/DVR benchmark and Channel-projected close-coupling solver.
Outputs 3D metrics, 2D cuts (R, gamma), (r, R), and publication figures.
"""
import argparse, json, sys, time
from pathlib import Path
import numpy as np
import scipy.linalg
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.reference import LegendreAngularDVR, _fd_second
from src.scatter_extract import (
    fd2_4th_dirichlet, fit_wavevector, free_incoming, free_outgoing
)
from src import bkmp2

BOHR = bkmp2.BOHR
M_H = 1.007825
mu_r = 0.5 * M_H
mu_R = (2.0 / 3.0) * M_H
Br = 16.857633 / mu_r     # 33.4535 cm^-1 A^2
BR = 16.857633 / mu_R     # 25.0901 cm^-1 A^2


def _channel_labels(r, gamma, U_d, j2mat, nr, nG):
    """Label each diatomic channel vector by (v, j).

    j from the expectation value of j^2 (exact for a DVR eigenvector of the
    (r, gamma) Hamiltonian), v from the number of sign changes of the radial
    profile at gamma ~ 90 deg (noise-robust: only sign changes with
    |profile| > 1e-6 max are counted).
    """
    n_red = U_d.shape[1]
    jlab = np.zeros(n_red, dtype=int)
    vlab = np.zeros(n_red, dtype=int)
    ig_mid = nG // 2
    for c in range(n_red):
        M = U_d[:, c].reshape(nr, nG)
        nrm2 = np.vdot(M, M).real
        j2 = float(np.vdot(M, M @ j2mat).real / nrm2)
        jlab[c] = int(round((-1.0 + np.sqrt(1.0 + 4.0 * max(j2, 0.0))) / 2.0))
        prof = np.abs(M[1:nr - 1, ig_mid])
        keep = prof > 1e-6 * prof.max() if prof.max() > 0 else prof > 0
        prof = M[1:nr - 1, ig_mid].real[keep]
        vlab[c] = int(np.sum(np.diff(np.sign(prof)) != 0))
    return jlab, vlab


def evaluate_V_3d(r, R, gamma):
    """BKMP2 potential (cm^-1) on an (r, R, gamma) product grid, with caching.

    Delegates to :mod:`src.bkmp2`, which locates or builds the native library
    and caches the evaluated grid under results/_cache/.
    """
    return bkmp2.potential_grid_rRg(r, R, gamma)

def solve_3d_scattering(Ec, nr=40, nR=120, nG=16, r_min=0.40, r_max=1.9,
                        R_min=0.7, R_max=8.5, legacy_fd=False):
    """Solve the 3D (r, R, gamma) scattering problem.

    The reference uses a *consistent* 4th-order finite-difference kinetic
    operator (:func:`src.scatter_extract.fd2_4th_dirichlet`): one-sided
    4th-order rows adjacent to every boundary, exact Dirichlet at R_min and at
    both r-ends, and the 4th-order one-sided multichannel Robin condition at
    R_max.  ``legacy_fd=True`` restores the older truncated 5-point operator
    (for before/after comparisons only)."""
    t0 = time.time()
    r = np.linspace(r_min, r_max, nr); hr = r[1] - r[0]
    R = np.linspace(R_min, R_max, nR); hR = R[1] - R[0]
    ang = LegendreAngularDVR(nG)
    gamma = np.arccos(ang.x)
    n_2d = nr * nG

    if legacy_fd:
        D2_r = _fd_second(nr, hr, 4)
        D2_R = _fd_second(nR, hR, 4)
    else:
        D2_r = fd2_4th_dirichlet(nr, hr)
        D2_R = fd2_4th_dirichlet(nR, hR)

    # 1. 3D Potential evaluation
    Vg = evaluate_V_3d(r, R, gamma)

    # 2. Asymptotic Diatomic Hamiltonian at R_max
    V_asym_slice = Vg[-1, :, :]
    if legacy_fd:
        H_d = sp.kron(-Br * D2_r, sp.identity(nG), format='csr')
        H_d += sp.kron(sp.diags(Br / r**2), sp.csr_matrix(ang.j2), format='csr')
        H_d += sp.diags(V_asym_slice.ravel())
        evals_d, evecs_d = np.linalg.eigh(H_d.toarray())
        idx_s = np.argsort(evals_d)
        evals_d = evals_d[idx_s]
        U_d = evecs_d[:, idx_s]
        full_of_red = np.arange(n_2d)
    else:
        # Dirichlet at both r-ends: diagonalize on the interior r-points only,
        # then pad the channel vectors with zeros at the walls.  (Building H_d
        # on the full grid would leave its wall rows empty, and eigh would
        # return the zero-eigenspace instead of the physical H2 channels.)
        ri = np.arange(1, nr - 1)
        D2_ii = np.asarray(D2_r[1:nr - 1, 1:nr - 1].todense())
        H_red = (np.kron(-Br * D2_ii, np.eye(nG))
                 + np.kron(np.diag(Br / r[ri] ** 2), np.asarray(ang.j2))
                 + np.diag(V_asym_slice[ri, :].ravel()))
        # D2_r is non-symmetric (one-sided boundary rows): use a general
        # eigensolver.  np.linalg.eigh would silently symmetrize it through
        # one triangle and return a wrong channel spectrum.
        evals_r, evecs_r = scipy.linalg.eig(H_red)
        order = np.argsort(evals_r.real)
        imag_worst = float(np.max(np.abs(evals_r.imag[order[:8]])))
        if imag_worst > 1e-8 * max(1.0, abs(evals_r.real[order[0]])):
            raise RuntimeError('diatomic eigenproblem returned a complex '
                               'spectrum - check fd2_4th_dirichlet')
        evals_r = evals_r.real[order]
        evecs_r = evecs_r[:, order].real
        idx_s = np.argsort(evals_r)
        evals_d = evals_r[idx_s]
        n_red = (nr - 2) * nG
        U_d = np.zeros((n_2d, n_red))
        for ir in range(1, nr - 1):
            U_d[ir * nG:(ir + 1) * nG, :] = evecs_r[(ir - 1) * nG:ir * nG, idx_s]
        full_of_red = np.array([ir * nG + ig for ir in range(1, nr - 1)
                                for ig in range(nG)])

    E_tot = evals_d[0] + Ec
    K_diag = [np.sqrt((E_tot - evals_d[i]) / BR) if E_tot > evals_d[i]
              else 1j * np.sqrt((evals_d[i] - E_tot) / BR) for i in range(len(evals_d))]
    K_diag = np.array(K_diag, dtype=complex)
    if legacy_fd:
        K_mat_2d = U_d @ np.diag(K_diag) @ U_d.T
    else:
        # The consistent 4th-order operator is NON-NORMAL (one-sided boundary
        # rows), so its right eigenvectors are not orthogonal (up to ~13% for
        # high channels).  Using U diag(K) U^T would then inject spurious
        # inter-channel couplings of order 10 A^-1 into the radiation
        # condition, which drives a spurious exponential boundary layer in the
        # closed channels at R_max.  The biorthogonal spectral decomposition
        # V diag(K) V^{-1} (left eigenvectors = rows of V^{-1}) is exact.
        # V_red^{-1} lifted to the full DVR space:  (V^{-1} P^T)
        V_red_inv = np.linalg.inv(evecs_r[:, idx_s])
        Vinv_full = np.zeros((n_red, n_2d), dtype=complex)
        Vinv_full[:, full_of_red] = V_red_inv
        K_mat_2d = U_d @ np.diag(K_diag) @ Vinv_full

    # 3. Full 3D Kronecker System
    T_R = sp.kron(-BR * D2_R, sp.identity(n_2d), format='csr')
    T_r = sp.kron(sp.identity(nR), sp.kron(-Br * D2_r, sp.identity(nG)), format='csr')
    C_R = sp.kron(sp.diags(BR / R**2), sp.kron(sp.identity(nr), sp.csr_matrix(ang.j2)), format='csr')
    C_r = sp.kron(sp.identity(nR), sp.kron(sp.diags(Br / r**2), sp.csr_matrix(ang.j2)), format='csr')
    H_3d = T_R + T_r + C_R + C_r + sp.diags(Vg.ravel())

    A = (H_3d - E_tot * sp.identity(H_3d.shape[0])).astype(complex).tolil()
    b = np.zeros(nR * n_2d, dtype=complex)

    # Inner Dirichlet BC at R_min
    for p in range(n_2d):
        A[p, :] = 0.0; A[p, p] = 1.0; b[p] = 0.0

    # Exact Dirichlet at both r-ends (rows are empty with the consistent operator)
    if not legacy_fd:
        for iR in range(nR):
            for ir in (0, nr - 1):
                for ig in range(nG):
                    row = iR * n_2d + ir * nG + ig
                    A[row, :] = 0.0
                    A[row, row] = 1.0

    # Asymptotic Multichannel Robin BC at R_max
    if legacy_fd:
        coef = np.array([3.0, -4.0, 1.0]) / (2.0 * hR)
    else:
        coef = np.array([25.0, -48.0, 36.0, -16.0, 3.0]) / (12.0 * hR)
    K0 = K_diag[0].real
    source_col = -2j * K0 * np.exp(-1j * K0 * R[-1]) * U_d[:, 0]
    j_end = nR - 1
    for a_ch, p in enumerate(full_of_red):
        row = j_end * n_2d + p
        A[row, :] = 0.0
        for b_ch, q in enumerate(full_of_red):
            A[row, j_end * n_2d + q] = -1j * K_mat_2d[a_ch, b_ch]
        for t, c in enumerate(coef):
            A[row, (j_end - t) * n_2d + p] += c
        b[row] = source_col[p]

    # Reference solve
    chi_3d_ref = spla.spsolve(A.tocsc(), b).reshape(nR, nr, nG)

    # 4. Channel-projected close-coupling solve
    #    Close-coupling subspace: all (v, j) channels with v in {0,1} and even
    #    j <= 6 (the collision-induced vibrational distortion is carried by the
    #    v = 1 ladder; higher v and j > 6 stay negligible).
    if legacy_fd:
        idx_ch = [0, 2, 4, 6]
    else:
        jlab, vlab = _channel_labels(r, gamma, U_d, ang.j2, nr, nG)
        idx_ch = [c for c in range(U_d.shape[1])
                  if vlab[c] in (0, 1) and jlab[c] in (0, 2, 4, 6)]
    U_sub = U_d[:, idx_ch]
    P = sp.kron(sp.identity(nR), U_sub, format='csr')
    A_proj = P.T @ A.tocsr() @ P
    b_proj = P.T @ b

    psi_proj = spla.spsolve(A_proj.tocsc(), b_proj)
    chi_cc_3d = (P @ psi_proj).reshape(nR, nr, nG)

    elapsed = time.time() - t0
    l2_error = float(np.linalg.norm(chi_cc_3d - chi_3d_ref) / np.linalg.norm(chi_3d_ref))

    print(f"3D Solve for Ec = {Ec:.1f} cm^-1 done in {elapsed:.1f}s. L2 Error = {l2_error*100:.4f}% (vs unconstrained reference)")
    return r, R, gamma, chi_3d_ref, chi_cc_3d, K0, l2_error, elapsed

def extract_s00_3d(chi, R, K0, r, gamma, nr, nG):
    """Extracts the elastic S-matrix element S_00 from a 3D scattering wavefunction."""
    hR = R[1] - R[0]
    hr = r[1] - r[0]
    ang = LegendreAngularDVR(nG)
    V_slice = bkmp2.potential_grid_rRg(r, np.array([R[-1]]), gamma, verbose=False).reshape(nr, nG)
    D2 = fd2_4th_dirichlet(nr, hr)
    ri = np.arange(1, nr - 1)
    D2i = np.asarray(D2[1:nr - 1, 1:nr - 1].todense())
    H = (np.kron(-Br * D2i, np.eye(nG))
         + np.kron(np.diag(Br / r[ri] ** 2), ang.j2)
         + np.diag(V_slice[ri, :].ravel()))
    evals, evecs = scipy.linalg.eig(H)
    idx = np.argsort(evals.real)
    v0 = np.zeros(nr * nG)
    for ir in range(1, nr - 1):
        v0[ir * nG:(ir + 1) * nG] = evecs[(ir - 1) * nG:ir * nG, idx[0]].real

    q = fit_wavevector(K0, hR, 4)
    lam = 2 * np.pi / q
    win = min(0.45 * (R[-1] - R[0]), max(3 * lam, 2.0))
    m = (R > R[-1] - win) & (R < R[-1] - 5 * hR)
    f0 = chi.reshape(len(R), -1) @ v0
    M = np.column_stack([free_incoming(q * R[m], 0), free_outgoing(q * R[m], 0)])
    c, *_ = np.linalg.lstsq(M, f0[m], rcond=None)
    S00 = c[1] / c[0]
    return complex(S00)


def plot_3d_slices(r, R, gamma, chi_ref, chi_pred, Ec, K0, l2_err, out_path):
    fig, axes = plt.subplots(2, 3, figsize=(15, 8.5))
    
    # Slice at r = re = 0.74144 A (equilibrium bond of H2): index nearest to 0.74144
    idx_r = np.argmin(np.abs(r - 0.74144))
    r_val = r[idx_r]
    gam_deg = gamma * 180.0 / np.pi
    
    chi_ref_cut = chi_ref[:, idx_r, :]  # (nR, nG)
    chi_pred_cut = chi_pred[:, idx_r, :]
    
    vmax = max(np.abs(chi_ref_cut).max(), np.abs(chi_pred_cut).max())
    
    # Row 0: Reference Solution at r = re
    im0 = axes[0, 0].pcolormesh(R, gam_deg, np.abs(chi_ref_cut).T, cmap='viridis', shading='auto', vmin=0, vmax=vmax)
    axes[0, 0].set_title(f'Reference 3D FD: $|\\chi(R, r_e, \\gamma)|$')
    axes[0, 0].set_ylabel(r'$\gamma$ (deg)')
    fig.colorbar(im0, ax=axes[0, 0], fraction=0.046, pad=0.04)
    
    im1 = axes[0, 1].pcolormesh(R, gam_deg, np.real(chi_ref_cut).T, cmap='coolwarm', shading='auto', vmin=-vmax, vmax=vmax)
    axes[0, 1].set_title(r'Reference 3D FD: $\mathrm{Re}[\chi]$')
    fig.colorbar(im1, ax=axes[0, 1], fraction=0.046, pad=0.04)
    
    # Cut in (r, R) plane at collinear geometry gamma = 0
    chi_ref_rR = chi_ref[:, :, 0] # (nR, nr)
    im2 = axes[0, 2].pcolormesh(R, r, np.abs(chi_ref_rR).T, cmap='viridis', shading='auto')
    axes[0, 2].set_title(r'Reference 3D FD: $|\chi(R, r, \gamma=0)|$')
    axes[0, 2].set_ylabel(r'$r$ (Å)')
    fig.colorbar(im2, ax=axes[0, 2], fraction=0.046, pad=0.04)
    
    # Row 1: Physics-Informed Channel Solver
    im3 = axes[1, 0].pcolormesh(R, gam_deg, np.abs(chi_pred_cut).T, cmap='viridis', shading='auto', vmin=0, vmax=vmax)
    axes[1, 0].set_title(rf'Close-coupling: $|\chi(R, r_e, \gamma)|$ (Rel $L^2$: {l2_err*100:.3f}%)')
    axes[1, 0].set_xlabel(r'$R$ (Å)')
    axes[1, 0].set_ylabel(r'$\gamma$ (deg)')
    fig.colorbar(im3, ax=axes[1, 0], fraction=0.046, pad=0.04)
    
    im4 = axes[1, 1].pcolormesh(R, gam_deg, np.real(chi_pred_cut).T, cmap='coolwarm', shading='auto', vmin=-vmax, vmax=vmax)
    axes[1, 1].set_title(r'Close-coupling: $\mathrm{Re}[\chi]$')
    axes[1, 1].set_xlabel(r'$R$ (Å)')
    fig.colorbar(im4, ax=axes[1, 1], fraction=0.046, pad=0.04)
    
    # Pointwise difference in (R, gamma) plane
    diff_cut = np.abs(chi_pred_cut - chi_ref_cut)
    im5 = axes[1, 2].pcolormesh(R, gam_deg, diff_cut.T, cmap='magma', shading='auto')
    axes[1, 2].set_title(r'Absolute Error $|\chi_{\mathrm{PINN}} - \chi_{\mathrm{FD}}|$')
    axes[1, 2].set_xlabel(r'$R$ (Å)')
    axes[1, 2].set_ylabel(r'$\gamma$ (deg)')
    fig.colorbar(im5, ax=axes[1, 2], fraction=0.046, pad=0.04)
    
    fig.suptitle(f'Full-Dimensional 3D H + H2(v=0) Reactive Scattering on BKMP2 PES (Ec = {Ec:.0f} cm^-1, K = {K0:.3f} A^-1, Error = {l2_err*100:.4f}%)', fontsize=12, y=0.98)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"Figure saved to {out_path}")

def main(argv=None):
    ap = argparse.ArgumentParser(
        description='Full-dimensional 3D H + H2(v=0) reactive scattering on BKMP2.')
    ap.add_argument('--energies', type=float, nargs='+', default=[400.0, 800.0],
                    help='collision energies Ec in cm^-1 (default: 400 800)')
    ap.add_argument('--nR', type=int, default=120, help='radial grid points')
    ap.add_argument('--nr', type=int, default=40, help='vibrational grid points')
    ap.add_argument('--nG', type=int, default=16, help='angular DVR points')
    ap.add_argument('--no-figures', action='store_true',
                    help='skip figure generation (headless / minimal installs)')
    args = ap.parse_args(argv)

    OUT = ROOT / 'results'
    OUT.mkdir(exist_ok=True)
    FIG = ROOT / 'article' / 'fig'
    FIG.mkdir(exist_ok=True, parents=True)

    records = []
    for Ec in args.energies:
        print(f"\n=======================================================")
        print(f"Full-Dimensional 3D A+BC Scattering: E_c = {Ec:.1f} cm^-1")
        print(f"=======================================================")
        r, R, gamma, chi_ref, chi_pred, K0, l2_err, el = solve_3d_scattering(
            Ec, nr=args.nr, nR=args.nR, nG=args.nG)

        asym_ref = float(np.mean(np.abs(chi_ref[R > 7.5])))
        asym_pred = float(np.mean(np.abs(chi_pred[R > 7.5])))

        # Extract elastic S-matrix elements for both solutions
        S00_ref = extract_s00_3d(chi_ref, R, K0, r, gamma, args.nr, args.nG)
        S00_pred = extract_s00_3d(chi_pred, R, K0, r, gamma, args.nr, args.nG)

        rec = {
            'E_c_cm1': Ec,
            'K0_inv_A': float(K0),
            'l2_relative_error_3d': l2_err,
            'asym_amplitude_fd': asym_ref,
            'asym_amplitude_pred': asym_pred,
            'S00_ref_abs': float(abs(S00_ref)),
            'S00_pred_abs': float(abs(S00_pred)),
            'S00_ref_arg_deg': float(np.degrees(np.angle(S00_ref))),
            'S00_pred_arg_deg': float(np.degrees(np.angle(S00_pred))),
            'P_elastic': float(abs(S00_ref)**2),
            'P_reaction_inelastic': float(max(1.0 - abs(S00_ref)**2, 0.0)),
            'elapsed_sec': el,
        }
        records.append(rec)
        print(f"  --> S00: ref = {abs(S00_ref):.6f} exp(i {rec['S00_ref_arg_deg']:.3f}°), "
              f"pred = {abs(S00_pred):.6f} exp(i {rec['S00_pred_arg_deg']:.3f}°), "
              f"Rel L2 = {l2_err*100:.4f}%")

        idx_re = np.argmin(np.abs(r - 0.74144))
        np.savez(OUT / f'scatter_3d_bkmp2_{int(Ec)}.npz',
                 r=r, R=R, gamma=gamma,
                 chi_ref_cut=chi_ref[:, idx_re, :],
                 chi_pred_cut=chi_pred[:, idx_re, :],
                 metrics=rec)

        if not args.no_figures:
            fig_path = FIG / f'fig_scatter_3d_bkmp2_{int(Ec)}.png'
            plot_3d_slices(r, R, gamma, chi_ref, chi_pred, Ec, K0, l2_err,
                           fig_path)

    with open(OUT / 'scatter_3d_bkmp2_summary.json', 'w') as f:
        json.dump(records, f, indent=2)

    ok = all(rec['l2_relative_error_3d'] < 0.01 for rec in records)
    print("\nAll full-dimensional 3D scattering computations completed "
          f"successfully (all errors < 1%: {ok})")


if __name__ == '__main__':
    main()
