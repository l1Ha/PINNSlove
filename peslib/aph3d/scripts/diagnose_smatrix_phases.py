#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Diagnose the relative phase / amplitude structure of the APH3D S-matrix for
H + H2(v=4,j=0) -> H + H2(para) against the Kendrick-Hazra-Balakrishnan
JCP 145, 164303 (2016) Table I structure:
  * f_N = S[inc -> (arr1, v'j')]           (non-reactive, arrangement 1)
  * f_R = (S[inc -> (arr2,v'j')] + S[inc -> (arr3,v'j')])/sqrt(2)  (even j')
  * eps = |f_R|/|f_N|  should be ~0.5-0.7 per Table I (|f_R|^2/|f_N|^2 = 0.26-0.45)
  * cos(Delta) = Re(f_N f_R*)/(|f_N||f_R|) should be ~0.99 at 1 uK
Checks also S-matrix symmetry, unitarity, and |S|^2 flux distribution.
"""
import sys
import numpy as np

SM = sys.argv[1] if len(sys.argv) > 1 else \
    "run_h3_kendrick2018/Output_Data/APH3D/APH3D_H3_DMBE/Output/Smatrix.txt"
if len(sys.argv) > 1 and sys.argv[1].endswith(".txt") is False:
    SM = sys.argv[1] + "/Output_Data/APH3D/APH3D_H3_DMBE/Output/Smatrix.txt"

# Parse first energy block only
qn = []
Srows = []
energies = []
with open(SM) as fh:
    cur_S = None
    in_qn = False
    for line in fh:
        s = line.strip()
        if s.startswith("# energy"):
            toks = s.split()
            energies.append(float(toks[2]))
            if len(energies) > 1:
                break
            in_qn = True
            cur_S = []
        elif in_qn and s.startswith("# QN"):
            if "chanl" in s:
                continue
            toks = s.split()
            qn.append(tuple(int(t) for t in toks[2:7]))
        elif s and not s.startswith("#"):
            re_, im_ = s.split()
            cur_S.append(complex(float(re_), float(im_)))
            if len(cur_S) == len(qn) ** 2:
                in_qn = False

S = np.array(cur_S).reshape(len(qn), len(qn))
E1 = energies[0]
nch = len(qn)
print(f"E_tot = {E1:.6e} Ha, Nch = {nch}")

# 1) unitarity per column
dmax = 0.0
for j in range(nch):
    col = np.sum(np.abs(S[:, j]) ** 2)
    dmax = max(dmax, abs(col - 1.0))
print(f"max |sum_j |S_ij|^2 - 1| = {dmax:.3e}")

# 2) symmetry S_ij vs S_ji
asym = np.max(np.abs(S - S.T)) / max(1e-30, np.max(np.abs(S)))
print(f"max asymmetry |S-S^T|/|S| = {asym:.3e}")

# 3) flux distribution: entrance channel candidates
inc = [i for i, q in enumerate(qn) if q == (1, 0, 4, 0, 0)]
print(f"incident channel indices (1,0,4,0,0): {[i+1 for i in inc]}")
a = inc[0]
row = S[a, :]
diag = abs(row[a]) ** 2
print(f"|S_aa|^2 (elastic prob) = {diag:.6f}")
# group flux by (arrangement, parity of j)
groups = {}
for i, q in enumerate(qn):
    arr, _, v, j, l = q
    key = (arr, "even" if j % 2 == 0 else "odd")
    groups.setdefault(key, 0.0)
    groups[key] += abs(row[i]) ** 2
print("flux |S_a->b|^2 by (arr, j-parity):")
for k in sorted(groups):
    print(f"   arr{k[0]} j_{k[1]}: {groups[k]:.6e}")

# 4) per-product-state eps, Delta for even j' (v',j') != (4,0)
print("\nper-state: (v',j')  |f_N|^2     |f_R|^2     eps=|fR/fN|  cosDelta  |f2-f3|/|f2+f3|")
lit_dom = [(0, 8), (0, 6), (0, 10), (1, 4), (1, 8), (2, 4), (2, 8), (3, 4), (3, 8), (3, 0), (2, 0), (1, 0), (0, 0)]
idx = {}
for i, q in enumerate(qn):
    arr, e, v, j, l = q
    idx.setdefault((arr, v, j), i)
tot_p_even = 0.0
cos_list = []
for (v, j) in lit_dom:
    if j % 2:
        continue
    k1 = (1, v, j)
    k2 = (2, v, j)
    k3 = (3, v, j)
    if k1 not in idx or k2 not in idx or k3 not in idx:
        continue
    fN = row[idx[k1]]
    f2 = row[idx[k2]]
    f3 = row[idx[k3]]
    fR = (f2 + f3) / np.sqrt(2)
    if abs(fN) < 1e-30:
        print(f"   ({v},{j:2d})  f_N ~ 0  |f_R|^2={abs(fR)**2:.3e}")
        continue
    eps = abs(fR) / abs(fN)
    cosd = (fN * np.conj(fR)).real / (abs(fN) * abs(fR)) if abs(fR) > 1e-30 else np.nan
    pauli = abs(f2 - f3) / abs(f2 + f3) if abs(f2 + f3) > 1e-30 else np.inf
    print(f"   ({v},{j:2d})  {abs(fN)**2:.4e}  {abs(fR)**2:.4e}  {eps:8.4f}  {cosd:8.4f}  {pauli:10.4f}")
    tot_p_even += abs(fN) ** 2 + abs(fR) ** 2
    if abs(fR) > 1e-30:
        cos_list.append(cosd)

print(f"\nsum even-j' (|fN|^2+|fR|^2) over listed states = {tot_p_even:.4e}")
if cos_list:
    print(f"median cos(Delta) = {np.nanmedian(cos_list):.4f}, mean = {np.nanmean(cos_list):.4f}")

# 5) incident-channel k-vector sanity: energy of each channel relative to Etot
# channel energy = E_ch; open if Etot > E_ch. Order in QN dump should follow
# ascending threshold. Print arrangement-1 thresholds ordering implied by dump.
print("\narrangement-1 QN dump order (v,j):", [(q[2], q[3]) for q in qn if q[0] == 1])
