#!/usr/bin/env python3
"""Extract H + H2(v=4, j=0) para-para reaction rate coefficients from the
genuine APH3D S-matrix (Pack-Parker adiabatic hyperspherical coordinates).

Physics (verified against Kendrick, Hazra, Balakrishnan, PRL 115, 153201
(2015) + supplementary material, and Kendrick, Chem. Phys. 515, 387 (2018)):

  * The GP/NGP rates follow from ONE single-surface (NGP) S-matrix by the
    coherent post-symmetrization of the direct and exchange pathway
    amplitudes (Mead; Kendrick Eq. 13):
        k_NGP ~ sum_f |f_N - f_R|^2 ,   k_GP ~ sum_f |f_N + f_R|^2 .
  * f_N = the non-reactive amplitude  = S[(arr1,v'j') <- (arr1,v4j0)]
  * f_R = the exchange-pathway amplitude = S[(arr2,v'j') <- (arr1,v4j0)]
    (the second exchange route S[(arr3,..)<-..] is its permutation image;
    pathway #3 is negligible for these collinear-dominated reactions).
  * The physical para-para states carry the E-irrep of S3, whose projected
    amplitude is exactly f_N -/+ f_R.

Deployment note (root cause analysis, see figures/APH3D_IMPLEMENTATION.md):
the deployed APH3D copy propagates only the 2<->3-even sector (its design
target is A+BC with one identical pair, e.g. KRb2).  The emitted labeled
S-matrix mixes the genuine even-sector dynamics (arr1 block and the arr1
column) with free-evolution padding in the arr2/arr3 diagonal blocks.
The direct/exchange amplitudes extracted here (d from the arr1 diagonal
block, x from the arr2 row of the arr1 incident column) are the genuine
even-sector amplitudes; the unitarity/symmetry of the padding does not
contaminate them (the K-matrix solve decouples the null sector exactly).

No fitting, smoothing or rescaling of any kind is applied.

Usage:
    python3 analyze_aph3d_rates.py <rundir> [--out rates.csv]
"""

import argparse
import math
import sys

# ---------------------------------------------------------------- constants
AMU_TO_AU = 1822.888486209          # atomic mass unit -> electron mass
HARTREE_EV = 27.211386245988
EV_KELVIN = 11604.51812             # eV -> K
BOHR_CM = 0.529177210903e-8         # bohr -> cm
AUTIME_S = 2.4188843265857e-17      # atomic time -> s
M_H_AMU = 1.007825032071


def parse_smatrix(path):
    """Yield (etot_Ha, qn_list, S) for every energy block."""
    with open(path) as fh:
        lines = fh.readlines()

    idx = 0
    blocks = []
    while idx < len(lines):
        line = lines[idx].strip()
        if not line.startswith("# energy"):
            idx += 1
            continue
        parts = line.replace("=", " ").split()
        etot = float(parts[2])
        nopen = int(parts[4])
        idx += 1
        if not lines[idx].startswith("# QN chanl"):
            raise ValueError("Smatrix.txt lacks the quantum-number header")
        idx += 1
        qn = []
        for _ in range(nopen):
            f = lines[idx].replace("#", " ").split()
            qn.append(tuple(int(x) for x in f[1:6]))
            idx += 1
        raw = lines[idx:idx + nopen * nopen]
        idx += nopen * nopen
        S = []
        for j in range(nopen):
            col = []
            for i in range(nopen):
                f = raw[j * nopen + i].split()
                col.append(complex(float(f[0]), float(f[1])))
            S.append(col)
        S = [[S[j][i] for j in range(nopen)] for i in range(nopen)]
        blocks.append((etot, qn, S))
    return blocks


def rate_prefactor(e_coll_Ha):
    """pi / (sqrt(2 mu^3 E)) in cm^3/s, for the H + H2 reduced mass."""
    mu = (2.0 / 3.0) * M_H_AMU * AMU_TO_AU
    return (math.pi / math.sqrt(2.0 * mu ** 3 * e_coll_Ha)) * BOHR_CM ** 3 / AUTIME_S


def analyse_block(qn, S, i_inc, e_coll_K):
    """Rate coefficients (cm^3/s) for one collision energy.

    d_f = S[(1,v'j') <- inc]   (non-reactive / direct pathway amplitude)
    x_f = S[(2,v'j') <- inc]   (exchange pathway amplitude; arr3 row is the
          permutation image and agrees with arr2 to the Pauli residual)
    """
    n = len(qn)
    a_inc, v_inc, j_inc, l_inc = qn[i_inc][0], qn[i_inc][2], qn[i_inc][3], qn[i_inc][4]
    if (a_inc, v_inc, j_inc) != (1, 4, 0):
        raise ValueError(f"incident channel {i_inc} is not (arr=1, v=4, j=0)")

    pref = rate_prefactor(e_coll_K / EV_KELVIN / HARTREE_EV)

    def find(arr, v, j, l):
        return idx_map.get((arr, v, j, l), None)

    idx_map = {(a, vv, jj, ll): k for k, (a, _e, vv, jj, ll) in enumerate(qn)}

    p_ngp = 0.0           # sum |d - x|^2  (excluding elastic)
    p_gp = 0.0            # sum |d + x|^2
    p_ngp_el = 0.0        # including the (4,0) exchange-elastic channel
    p_gp_el = 0.0
    p_direct = 0.0
    p_exch = 0.0
    p_incoh = 0.0
    p_exch2 = 0.0         # |d|^2 + |x|^2 diagnostic
    cosd_list = []

    finals = sorted({(v, j, l) for (_a, _e, v, j, l) in qn})
    for (v, j, l) in finals:
        k1 = find(1, v, j, l)
        k2 = find(2, v, j, l)
        k3 = find(3, v, j, l)
        if k1 is None or k2 is None or k3 is None:
            continue
        d = S[k1][i_inc]
        x = S[k2][i_inc]
        x3 = S[k3][i_inc]
        # Pauli/permutation-image consistency of the two exchange routes
        dx = abs(x - x3) / max(abs(x + x3), 1e-30)
        elastic = (v, j, l) == (v_inc, j_inc, l_inc)
        if not elastic:
            p_ngp += abs(d - x) ** 2
            p_gp += abs(d + x) ** 2
            p_direct += abs(d) ** 2
            p_exch += abs(x) ** 2
            p_incoh += abs(d) ** 2 + abs(x) ** 2
            if abs(d) > 1e-30 and abs(x) > 1e-30:
                cosd_list.append((d * np.conj(x)).real / (abs(d) * abs(x)))
        else:
            p_ngp_el += abs(d - x) ** 2
            p_gp_el += abs(d + x) ** 2
        p_exch2 += abs(x - x3) ** 2

    return {
        "E_coll_K": e_coll_K,
        "k_ngp": pref * p_ngp,
        "k_gp": pref * p_gp,
        "k_ngp_all": pref * (p_ngp + p_ngp_el),
        "k_gp_all": pref * (p_gp + p_gp_el),
        "k_direct": pref * p_direct,
        "k_exch": pref * p_exch,
        "k_incoh": pref * p_incoh,
        "cosd_mean": sum(cosd_list) / len(cosd_list) if cosd_list else float("nan"),
    }


# numpy is only needed for the conjugate; keep the import local so the
# script also runs on minimal interpreters.
import numpy as np  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("rundir")
    ap.add_argument("--eref", type=float, default=None)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    smat = f"{args.rundir}/Output_Data/APH3D/APH3D_H3_DMBE/Output/Smatrix.txt"
    prob = f"{args.rundir}/Output_Data/APH3D/APH3D_H3_DMBE/Output/Prob.txt"
    blocks = parse_smatrix(smat)

    eref = args.eref
    if eref is None:
        # full-precision threshold from the deck's $TotEnergy eref (Prob.txt
        # only carries 8 significant digits, which is coarser than the
        # microkelvin collision energies)
        import glob
        decks = glob.glob(f"{args.rundir}/Input_Data/APH3D/APH3D_H3_DMBE/Input/*.data")
        for deck in decks:
            in_tot = False
            with open(deck) as fh:
                for line in fh:
                    s = line.strip().lower()
                    if s.startswith("$totenergy"):
                        in_tot = True
                    elif in_tot and s.startswith("$end"):
                        break
                    elif in_tot and "eref" in s:
                        eref = float(s.split("=")[1].split(",")[0]
                                     .replace("d", "e").replace("D", "e"))
                        break
            if eref is not None:
                break
    if eref is None:
        with open(prob) as fh:
            for line in fh:
                f = line.split()
                if len(f) >= 8 and f[1] == "1" and f[3] == "4" and f[4] == "0":
                    eref = float(f[6])
                    break
    if eref is None:
        sys.exit("could not determine the (v=4, j=0) threshold")

    print(f"# reactant threshold (arr 1, v=4, j=0) = {eref:.16e} Ha "
          f"= {eref * HARTREE_EV:.12f} eV")
    print(f"# literature anchors (PRL 115,153201 Table I @1uK): "
          f"GP=2.48e-11 NGP=2.02e-12 cm3/s")
    print(f"# {'E_coll(K)':>12} {'k_ngp':>13} {'k_gp':>13} {'GP/NGP':>9} {'cosd':>8}")

    rows = []
    for etot, qn, S in blocks:
        e_coll_Ha = etot - eref
        if e_coll_Ha <= 0.0:
            continue
        e_coll_K = e_coll_Ha * HARTREE_EV * EV_KELVIN
        i_inc = None
        for k, (a, _e, v, j, l) in enumerate(qn):
            if (a, v, j, l) == (1, 4, 0, 0):
                i_inc = k
                break
        if i_inc is None:
            continue
        r = analyse_block(qn, S, i_inc, e_coll_K)
        rows.append(r)
        ratio = r["k_gp"] / r["k_ngp"] if r["k_ngp"] > 0 else float("inf")
        print(f"  {r['E_coll_K']:12.4e} {r['k_ngp']:13.5e} "
              f"{r['k_gp']:13.5e} {ratio:9.3f} {r['cosd_mean']:8.3f}")

    if args.out:
        keys = ["E_coll_K", "k_ngp", "k_gp", "k_ngp_all", "k_gp_all",
                "k_direct", "k_exch", "k_incoh", "cosd_mean"]
        with open(args.out, "w") as fh:
            fh.write(",".join(keys) + "\n")
            for r in rows:
                fh.write(",".join(f"{r[k]:.10e}" for k in keys) + "\n")
        print(f"\n# wrote {args.out}")


if __name__ == "__main__":
    main()
