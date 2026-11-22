"""Angular momentum coupling and frame transformations for A + BC scattering with J > 0.

Provides rigorous representation theory for triatomic scattering with non-zero
total angular momentum J:
1. Body-Fixed (BF) Representation:
     Quantum numbers: (E, J, M, v, j, K, p)
     where K is the projection of J (and j) onto the body-fixed Jacobi vector R.
     Hamiltonian includes diagonal centrifugal terms and off-diagonal Coriolis
     coupling terms (-2 J . j) between K and K +/- 1.

2. Space-Fixed (SF) Representation:
     Quantum numbers: (E, J, M, v, j, l, p)
     where l is the orbital angular momentum of A relative to BC.
     Centrifugal barrier is strictly diagonal: B_R * l(l+1) / R^2.

3. Unitary Frame Transformation (Pack-Parker / Arthurs-Dalgarno):
     U(J, j, p) diagonalizes the Coriolis L^2 operator in the parity-adapted BF basis,
     yielding eigenvalues l(l+1) and the exact frame transformation matrix:
         |J M j l p> = sum_K U_{K, l}(J, j, p) |J M j K p>
         chi_SF = U^dagger chi_BF
         S_SF = U^dagger S_BF U.
"""
import numpy as np


class JAngularCoupling:
    """Coupled angular momentum representation for A + BC scattering with total J >= 0."""

    def __init__(self, J=0, j_max=4, mass_A=1.007825, mass_B=1.007825, mass_C=1.007825):
        self.J = int(J)
        self.j_max = int(j_max)
        self.mA = float(mass_A)
        self.mB = float(mass_B)
        self.mC = float(mass_C)

        # Reduced masses in amu
        self.mu_r = (self.mB * self.mC) / (self.mB + self.mC)
        self.mu_R = (self.mA * (self.mB + self.mC)) / (self.mA + self.mB + self.mC)
        H2_OVER_2AMU = 16.857633
        self.Br = H2_OVER_2AMU / self.mu_r   # cm^-1 A^2
        self.BR = H2_OVER_2AMU / self.mu_R   # cm^-1 A^2

    def get_channels(self, parity=+1, j_step=1):
        """Returns the list of SF channels (j, l) and BF channels (j, K) for given parity p (+1 or -1).

        Parity is the total space-inversion parity p = (-1)^(j+l).
        ``j_step=1`` includes all j; ``j_step=2`` includes only even j (homonuclear para-H2).
        """
        p = int(parity)
        sf_channels = []
        bf_channels = []

        for j in range(0, self.j_max + 1, j_step):
            # SF channels: |J, j, l> with |J - j| <= l <= J + j and (-1)^(j+l) == p
            for l in range(abs(self.J - j), self.J + j + 1):
                if (-1)**(j + l) == p:
                    sf_channels.append((j, l))

            # BF parity-adapted channels: |J, j, K, p>
            # K ranges from k_min to min(j, J)
            k_min = 0 if (p == (-1)**self.J) else 1
            k_max = min(j, self.J)
            for K in range(k_min, k_max + 1):
                bf_channels.append((j, K))

        return sf_channels, bf_channels

    def coriolis_matrix_j(self, j, parity=+1):
        """Constructs the L^2 = (J - j)^2 Coriolis matrix in the parity-adapted BF basis for a single j block.

        Returns (L2_matrix, k_list, U_matrix, l_list) where:
          - L2_matrix has diagonal J(J+1) + j(j+1) - 2K^2
          - Off-diagonal Coriolis coupling to K+1
          - U_matrix is the exact orthogonal matrix diagonalizing L^2:
              U.T @ L2 @ U = diag(l(l+1))
        """
        p = int(parity)
        k_min = 0 if (p == (-1)**self.J) else 1
        k_max = min(j, self.J)

        if k_min > k_max:
            return np.empty((0, 0)), [], np.empty((0, 0)), []

        k_list = list(range(k_min, k_max + 1))
        n_k = len(k_list)

        ctot = self.J * (self.J + 1)
        cent = j * (j + 1)

        L2 = np.zeros((n_k, n_k), dtype=float)
        for idx, k in enumerate(k_list):
            # Diagonal: J(J+1) + j(j+1) - 2K^2
            L2[idx, idx] = ctot + cent - 2.0 * k * k

            # Off-diagonal Coriolis coupling to K+1
            if idx + 1 < n_k and k_list[idx + 1] == k + 1:
                if k == 0:
                    off = -np.sqrt(2.0 * ctot * cent)
                else:
                    ckp = k * (k + 1)
                    off = -np.sqrt(max((ctot - ckp) * (cent - ckp), 0.0))
                L2[idx, idx + 1] = off
                L2[idx + 1, idx] = off

        # Diagonalize L2
        evals, U = np.linalg.eigh(L2)
        l_list = []
        for ev in evals:
            # ev = l(l+1) => l = sqrt(ev + 0.25) - 0.5
            l_val = int(np.round(np.sqrt(max(ev + 0.25, 0.0)) - 0.5))
            l_list.append(l_val)

        # Consistent phase convention (same as Pack-Parker and Fortran ABC code):
        # phase = (-1)**(j + k_max)
        phase = (-1)**(j + k_max)
        for col in range(n_k):
            if phase * U[-1, col] < 0.0:
                U[:, col] = -U[:, col]

        return L2, k_list, U, l_list

    def full_frame_transformation(self, parity=+1, j_step=1):
        """Builds the block-diagonal unitary frame transformation matrix between BF and SF

        for all channels up to j_max in the specified parity subspace.
        Returns:
            U_total: block-diagonal orthogonal matrix of shape (n_channels, n_channels)
            sf_channels: list of (j, l)
            bf_channels: list of (j, K)
            l_eigenvalues: array of l(l+1) values
        """
        sf_ch, bf_ch = self.get_channels(parity, j_step=j_step)
        n_tot = len(bf_ch)
        U_total = np.zeros((n_tot, n_tot), dtype=float)
        l_evals = np.zeros(n_tot, dtype=float)

        offset = 0
        for j in range(0, self.j_max + 1, j_step):
            L2, k_list, U_block, l_list = self.coriolis_matrix_j(j, parity)
            n_block = len(k_list)
            if n_block > 0:
                U_total[offset:offset + n_block, offset:offset + n_block] = U_block
                for b_idx, l_val in enumerate(l_list):
                    l_evals[offset + b_idx] = l_val * (l_val + 1)
                offset += n_block

        return U_total, sf_ch, bf_ch, l_evals

    def transform_bf_to_sf(self, chi_bf, U):
        """Transforms a radial multi-channel wavefunction from BF representation to SF:

        chi_SF = U.T @ chi_BF
        """
        # chi_bf can have shape (nR, n_channels) or (n_channels,)
        if chi_bf.ndim == 1:
            return U.T @ chi_bf
        elif chi_bf.ndim == 2:
            return chi_bf @ U
        else:
            raise ValueError(f"Unsupported dimension: {chi_bf.ndim}")

    def transform_s_matrix_bf_to_sf(self, S_bf, U):
        """Transforms the scattering S-matrix from BF representation to SF:

        S_SF = U.T @ S_BF @ U
        """
        return U.T @ S_bf @ U
