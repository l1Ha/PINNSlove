"""Test suite and numerical verification for J > 0 angular momentum coupling
and Body-Fixed (BF) vs Space-Fixed (SF) frame transformations.
"""
import numpy as np
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.j_frames import JAngularCoupling


def test_j_angular_momentum():
    print("================================================================")
    print("Testing J > 0 Coupled-Channel Representation & Frame Transforms")
    print("================================================================")

    for J in [0, 1, 2, 3]:
        print(f"\n--- Total Angular Momentum J = {J} ---")
        coupling = JAngularCoupling(J=J, j_max=4)

        for p in [+1, -1]:
            sf_ch, bf_ch = coupling.get_channels(parity=p)
            print(f"  Parity p = {p:+d}:")
            print(f"    SF channels (j, l) count={len(sf_ch)}: {sf_ch}")
            print(f"    BF channels (j, K) count={len(bf_ch)}: {bf_ch}")

            # Verify identical channel count
            assert len(sf_ch) == len(bf_ch), f"Count mismatch for J={J}, p={p}!"

            U, sf_c, bf_c, l_evals = coupling.full_frame_transformation(parity=p)
            if len(bf_ch) > 0:
                # 1. Verify strict unitarity / orthogonality of U
                eye_err = np.max(np.abs(U.T @ U - np.eye(len(bf_ch))))
                assert eye_err < 1e-12, f"U is not strictly orthogonal! Error: {eye_err}"

                # 2. Verify eigenvalues match l(l+1)
                for idx, (j, l) in enumerate(sf_c):
                    target_ev = l * (l + 1)
                    actual_ev = l_evals[idx]
                    assert abs(actual_ev - target_ev) < 1e-10, (
                        f"L^2 eigenvalue mismatch: got {actual_ev}, expected {target_ev}"
                    )

                print(f"    -> Orthogonality error: {eye_err:.2e}, l(l+1) spectrum verified!")

    print("\n[PASSED] All channel counts, L^2 spectra, and unitary transforms verified!")


def test_hamiltonian_equivalence():
    print("\n================================================================")
    print("Verifying Hamiltonian Operator Equivalence (H_SF == U^T H_BF U)")
    print("================================================================")

    J = 1
    p = -1
    coupling = JAngularCoupling(J=J, j_max=2)
    U, sf_ch, bf_ch, l_evals = coupling.full_frame_transformation(parity=p)
    n_ch = len(bf_ch)

    print(f"For J={J}, parity={p:+d}: {n_ch} channels")
    print(f"  BF channels: {bf_ch}")
    print(f"  SF channels: {sf_ch}")

    # Build Coriolis matrix in BF
    L2_bf = np.zeros((n_ch, n_ch))
    offset = 0
    for j in range(3):
        L2_block, k_list, _, _ = coupling.coriolis_matrix_j(j, p)
        nb = len(k_list)
        if nb > 0:
            L2_bf[offset:offset + nb, offset:offset + nb] = L2_block
            offset += nb

    # Centrifugal in SF
    L2_sf = np.diag(l_evals)

    # Check that U.T @ L2_bf @ U == L2_sf exactly!
    transformed = U.T @ L2_bf @ U
    diff = np.max(np.abs(transformed - L2_sf))
    print(f"Max difference between U.T @ L2_bf @ U and diagonal L2_sf: {diff:.3e}")
    assert diff < 1e-12, "Coriolis L^2 transformation failed!"

    # Test mock potential matrix V (diagonal in K in BF frame)
    # V_bf(j, K; j', K)
    V_bf = np.zeros((n_ch, n_ch))
    for idx1, (j1, K1) in enumerate(bf_ch):
        for idx2, (j2, K2) in enumerate(bf_ch):
            if K1 == K2:
                V_bf[idx1, idx2] = 100.0 / (1.0 + abs(j1 - j2))

    # Transform V to SF representation:
    V_sf = U.T @ V_bf @ U

    # Total coupled operator in BF and SF at R = 2.0 A:
    R = 2.0
    BR = coupling.BR
    H_bf = (BR / R**2) * L2_bf + V_bf
    H_sf = (BR / R**2) * L2_sf + V_sf

    # Compare eigenvalues:
    ev_bf = np.linalg.eigvalsh(H_bf)
    ev_sf = np.linalg.eigvalsh(H_sf)
    spec_diff = np.max(np.abs(ev_bf - ev_sf))
    print(f"Max difference between H_bf and H_sf energy eigenvalues: {spec_diff:.3e}")
    assert spec_diff < 1e-12, "Hamiltonian eigenvalues differ between representations!"
    print("[PASSED] Exact mathematical equivalence between BF and SF proven!")


if __name__ == '__main__':
    test_j_angular_momentum()
    test_hamiltonian_equivalence()
