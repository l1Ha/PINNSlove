#!/usr/bin/env python3
"""
Exact Channel-Resolved Quantum Scattering Extraction & Rate Computation.

Eliminates the critical discrepancy of confusing total loss (1 - |S_ii|^2) with pure reaction.
Distinguishes:
  1. P_elastic: |S_{ii}|^2 (Elastic scattering)
  2. P_quench / P_inelastic: sum_{f in arr 1, f != i} |S_{fi}|^2 (Pure non-reactive inelastic relaxation)
  3. P_react: sum_{f in arr 2, 3} |S_{fi}|^2 (True reactive exchange: A + BC -> B + AC or C + AB)
  4. P_loss: 1 - |S_{ii}|^2 = P_inelastic + P_react

Converts probabilities to physical rate coefficients k = (pi / (mu * k_coll)) * P * 6.1263e-9 cm^3/s.
"""

from __future__ import annotations

import argparse
import glob
import math
import os
import re
import sys
from pathlib import Path
import numpy as np

AU2CC = 6.12627e-9       # atomic units to cm^3/s
AMU_AU = 1822.8885       # amu to a.u.
EV_TO_K = 11604.51812    # eV to Kelvin
H2_AU = 27.211386        # Hartree in eV

MASS_TABLE = {
    1: 1.00782503223,
    2: 2.01410177812,
    3: 3.016049281,
    4: 4.00260325413,
    39: 38.9637064864,
    87: 86.9091805310,
}


def get_mu_au(masses: list[int]) -> float:
    """Computes the reduced mass in atomic units for reactant A + BC."""
    mA = MASS_TABLE.get(masses[0], masses[0] * 1.007825) * AMU_AU
    mB = MASS_TABLE.get(masses[1], masses[1] * 1.007825) * AMU_AU
    mC = MASS_TABLE.get(masses[2], masses[2] * 1.007825) * AMU_AU
    mBC = mB + mC
    return (mA * mBC) / (mA + mBC)


def rate_from_prob(E_coll_eV: np.ndarray, prob: np.ndarray, mu: float, jtot: int = 0) -> np.ndarray:
    """Computes reaction rate coefficient k(E) = (2J+1) * (pi / (mu * k_in)) * P in cm^3/s."""
    E_au = np.maximum(np.asarray(E_coll_eV, float), 1e-30) / H2_AU
    k_au = np.sqrt(2.0 * mu * E_au)
    pref = float(2 * jtot + 1) * np.pi / (mu * k_au)
    return AU2CC * pref * np.asarray(prob, float)


def parse_energy_diatomic(filepath: str | Path) -> list[dict]:
    """Parses energy_diatomic.dat to map channel index to (arr, v, j, k, E_eV)."""
    channels = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 6:
                try:
                    idx = int(parts[0])
                    arr = int(parts[1])
                    v = int(parts[2])
                    j = int(parts[3])
                    k = int(parts[4])
                    E_eV = float(parts[5])
                    channels.append({
                        'idx': idx, 'arr': arr, 'v': v, 'j': j, 'k': k, 'E_eV': E_eV
                    })
                except ValueError:
                    continue
    return channels


def find_incident_channel(channels: list[dict], target_v: int = 4, target_j: int = 0, target_k: int = 0) -> dict | None:
    """Finds the incident channel matching target (arr=1, v, j, k)."""
    for ch in channels:
        if ch['arr'] == 1 and ch['v'] == target_v and ch['j'] == target_j and ch['k'] == target_k:
            return ch
    # Fallback to arr=1 with specified v, lowest j
    cands = [ch for ch in channels if ch['arr'] == 1 and ch['v'] == target_v]
    if cands:
        return min(cands, key=lambda c: (abs(c['j'] - target_j), c['E_eV']))
    return None


def extract_from_sopen(sopen_file: str | Path, incident_v: int = 4, incident_j: int = 0, incident_k: int = 0) -> dict:
    """
    Parses a single S-open*.dat file and calculates channel-resolved transition probabilities:
      - P_elastic
      - P_inelastic (arr 1, f != in)
      - P_react_arr2 (arr 2)
      - P_react_arr3 (arr 3)
      - P_react_total (arr 2 + 3)
      - P_total_loss (1 - P_elastic)
    """
    p_elastic = 0.0
    p_inelastic = 0.0
    p_arr2 = 0.0
    p_arr3 = 0.0
    found_incident = False

    with open(sopen_file, 'r') as f:
        for line in f:
            parts = line.split()
            if len(parts) < 11:
                continue
            ii, ni, ji, ki = int(parts[0]), int(parts[1]), int(parts[2]), int(parts[3])
            iif, nf, jf, kf = int(parts[4]), int(parts[5]), int(parts[6]), int(parts[7])
            tp = float(parts[10])

            # Filter for rows originating from the specified incident state
            if ii == 1 and ni == incident_v and ji == incident_j and ki == incident_k:
                found_incident = True
                if iif == 1:
                    if nf == incident_v and jf == incident_j and kf == incident_k:
                        p_elastic += tp
                    else:
                        p_inelastic += tp
                elif iif == 2:
                    p_arr2 += tp
                elif iif == 3:
                    p_arr3 += tp

    p_react = p_arr2 + p_arr3
    p_total_loss = 1.0 - p_elastic if found_incident else 0.0

    return {
        'found': found_incident,
        'P_elastic': p_elastic,
        'P_inelastic': p_inelastic,
        'P_react_arr2': p_arr2,
        'P_react_arr3': p_arr3,
        'P_react': p_react,
        'P_total_loss': p_total_loss,
        'unitarity_sum': p_elastic + p_inelastic + p_react
    }


def extract_from_smatrix_complex(smat_file: str | Path, channels: list[dict], npes: int = 2,
                                 incident_v: int = 4, incident_j: int = 0, incident_k: int = 0) -> list[dict]:
    """
    Streams S-matrix-complex.dat and extracts channel-resolved probabilities for each energy:
      - P_elastic
      - P_inelastic (arr 1, non-elastic)
      - P_react_arr2 (arr 2)
      - P_react_arr3 (arr 3)
      - P_react (arr 2 + 3)
      - P_total_loss (1 - P_elastic)
    """
    ch_in = find_incident_channel(channels, target_v=incident_v, target_j=incident_j, target_k=incident_k)
    if not ch_in:
        return []

    i0 = ch_in['idx']
    ic = (i0 - 1) * max(1, npes) + 1  # 1-based incident electronic channel
    n_spatial = len(channels)
    nchan = n_spatial * max(1, npes)

    results = []
    cur_E_eV = None
    in_block = False
    row_count = 0

    with open(smat_file, 'r') as f:
        for line in f:
            if line.startswith('E(eV)='):
                m = re.search(r'E\(eV\)=\s*([-\d.eE+]+)', line)
                if m:
                    cur_E_eV = float(m.group(1))
                in_block = False
                row_count = 0
            elif line.startswith('S-matrix complex block'):
                in_block = True
                row_count = 0
            elif in_block:
                row_count += 1
                if row_count == ic:
                    parts = [float(x) for x in line.split()]
                    p_elastic = 0.0
                    p_inelastic = 0.0
                    p_arr2 = 0.0
                    p_arr3 = 0.0

                    for jc in range(1, min(len(parts) // 2, nchan) + 1):
                        re_s = parts[2 * (jc - 1)]
                        im_s = parts[2 * (jc - 1) + 1]
                        tp = re_s**2 + im_s**2

                        j_spatial = (jc - 1) // max(1, npes) + 1
                        if j_spatial <= len(channels):
                            arr = channels[j_spatial - 1]['arr']
                        else:
                            arr = 1

                        if jc == ic:
                            p_elastic += tp
                        elif arr == 1:
                            p_inelastic += tp
                        elif arr == 2:
                            p_arr2 += tp
                        elif arr == 3:
                            p_arr3 += tp

                    p_react = p_arr2 + p_arr3
                    p_loss = max(0.0, 1.0 - p_elastic)
                    results.append({
                        'found': True,
                        'E_tot_eV': cur_E_eV,
                        'P_elastic': p_elastic,
                        'P_inelastic': p_inelastic,
                        'P_react_arr2': p_arr2,
                        'P_react_arr3': p_arr3,
                        'P_react': p_react,
                        'P_total_loss': p_loss,
                        'unitarity_sum': p_elastic + p_inelastic + p_react
                    })
                    in_block = False  # Skip remaining rows of this block

    return results


def analyze_run_directory(rundir: str | Path, v_in: int = 4, j_in: int = 0, k_in: int = 0) -> dict | None:
    """Analyzes a full run directory containing output/S-open*.dat and input.d."""
    rundir = Path(rundir)
    outdir = rundir / 'output'
    if not outdir.exists():
        outdir = rundir

    # Parse input.d
    input_file = rundir / 'input.d'
    if not input_file.exists():
        input_file = rundir.parent / 'input.d'

    masses = [1, 1, 1]
    jtot = 0
    enrg = 2.161
    dnrg = 1e-7
    nnrg = 100

    if input_file.exists():
        txt = input_file.read_text()
        m_mass = re.search(r'^\s*mass\s*=\s*([\d\s,]+)', txt, re.M)
        if m_mass:
            masses = [int(x) for x in m_mass.group(1).replace(',', ' ').split()]
        m_jtot = re.search(r'^\s*jtot\s*=\s*(\d+)', txt, re.M)
        if m_jtot:
            jtot = int(m_jtot.group(1))
        m_enrg = re.search(r'^\s*enrg\s*=\s*([-\d.eE+]+)', txt, re.M)
        if m_enrg:
            enrg = float(m_enrg.group(1))
        m_dnrg = re.search(r'^\s*dnrg\s*=\s*([-\d.eE+]+)', txt, re.M)
        if m_dnrg:
            dnrg = float(m_dnrg.group(1))
        m_nnrg = re.search(r'^\s*nnrg\s*=\s*(\d+)', txt, re.M)
        if m_nnrg:
            nnrg = int(m_nnrg.group(1))

    mu = get_mu_au(masses)

    # Check diatomic energy
    diatomic_file = outdir / 'energy_diatomic.dat'
    e_incident = None
    if diatomic_file.exists():
        channels = parse_energy_diatomic(diatomic_file)
        ch_in = find_incident_channel(channels, target_v=v_in, target_j=j_in, target_k=k_in)
        if ch_in:
            e_incident = ch_in['E_eV']

    sopen_files = sorted(glob.glob(str(outdir / 'S-open*.dat')), key=lambda p: int(re.search(r'S-open(\d+)\.dat', p).group(1)) if re.search(r'S-open(\d+)\.dat', p) else 0)
    smat_complex_file = outdir / 'S-matrix-complex.dat'

    results = []
    if sopen_files:
        for fpath in sopen_files:
            idx_match = re.search(r'S-open(\d+)\.dat', fpath)
            idx = int(idx_match.group(1)) if idx_match else 1
            res = extract_from_sopen(fpath, incident_v=v_in, incident_j=j_in, incident_k=k_in)
            if res['found']:
                tot_E = enrg + (idx - 1) * dnrg
                coll_E = tot_E - e_incident if e_incident is not None else (idx - 1) * dnrg
                res['E_tot_eV'] = tot_E
                res['E_coll_eV'] = max(coll_E, 1e-12)
                res['E_coll_K'] = res['E_coll_eV'] * EV_TO_K
                res['k_react'] = float(rate_from_prob(np.array([res['E_coll_eV']]), np.array([res['P_react']]), mu, jtot)[0])
                res['k_inelastic'] = float(rate_from_prob(np.array([res['E_coll_eV']]), np.array([res['P_inelastic']]), mu, jtot)[0])
                res['k_total_loss'] = float(rate_from_prob(np.array([res['E_coll_eV']]), np.array([res['P_total_loss']]), mu, jtot)[0])
                results.append(res)
    elif smat_complex_file.exists() and diatomic_file.exists():
        raw_results = extract_from_smatrix_complex(smat_complex_file, channels, npes=2,
                                                   incident_v=v_in, incident_j=j_in, incident_k=k_in)
        for idx, res in enumerate(raw_results, 1):
            tot_E = res['E_tot_eV'] if res['E_tot_eV'] is not None else (enrg + (idx - 1) * dnrg)
            coll_E = tot_E - e_incident if e_incident is not None else (idx - 1) * dnrg
            res['E_coll_eV'] = max(coll_E, 1e-12)
            res['E_coll_K'] = res['E_coll_eV'] * EV_TO_K
            res['k_react'] = float(rate_from_prob(np.array([res['E_coll_eV']]), np.array([res['P_react']]), mu, jtot)[0])
            res['k_inelastic'] = float(rate_from_prob(np.array([res['E_coll_eV']]), np.array([res['P_inelastic']]), mu, jtot)[0])
            res['k_total_loss'] = float(rate_from_prob(np.array([res['E_coll_eV']]), np.array([res['P_total_loss']]), mu, jtot)[0])
            results.append(res)
    else:
        return None

    return {
        'rundir': str(rundir),
        'masses': masses,
        'jtot': jtot,
        'mu_au': mu,
        'e_incident_eV': e_incident,
        'data': results
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Extract channel-resolved quantum scattering rates.")
    parser.add_argument("target_dir", nargs="?", default=".", help="Run directory to analyze")
    parser.add_argument("--v", type=int, default=4, help="Incident vibrational quantum number v (default: 4)")
    parser.add_argument("--j", type=int, default=0, help="Incident rotational quantum number j (default: 0)")
    args = parser.parse_args()

    res = analyze_run_directory(args.target_dir, v_in=args.v, j_in=args.j)
    if res and res['data']:
        print(f"Successfully analyzed {len(res['data'])} energy points.")
        d0 = res['data'][0]
        print(f"Sample low energy point: E_coll = {d0['E_coll_K']:.4e} K")
        print(f"  P_elastic     = {d0['P_elastic']:.6e}")
        print(f"  P_inelastic   = {d0['P_inelastic']:.6e} (k = {d0['k_inelastic']:.4e} cm^3/s)")
        print(f"  P_react       = {d0['P_react']:.6e} (k = {d0['k_react']:.4e} cm^3/s)")
        print(f"  P_total_loss  = {d0['P_total_loss']:.6e} (k = {d0['k_total_loss']:.4e} cm^3/s)")
        print(f"  Ratio (Total Loss / Pure React) = {d0['P_total_loss'] / max(d0['P_react'], 1e-30):.2f}x")
    else:
        print("No valid S-open data found for the specified state.")
