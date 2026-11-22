#!/usr/bin/env python3
"""
Test S3 permutation covariance of the APH3D labeled-basis S-matrix.

For the NGP single-surface H3 Hamiltonian (fully symmetric under all
permutations of the three identical H nuclei), the labeled S-matrix must
satisfy   S[p(a) <- p(f)] = S[a <- f]  for every permutation p.

Concretely, with incident channels a1=(1,4,0), a2=(2,4,0), a3=(3,4,0):
    S[(1,v'j') <- a1]  must equal  S[(2,v'j') <- a2]  and  S[(3,v'j') <- a3]
    S[(2,v'j') <- a1]  must equal  S[(3,v'j') <- a1]  and  S[(1,v'j') <- a2] ...
(The last line reflects the transposition structure: elements connecting
 different arrangements are all equal to the same exchange amplitude x.)
"""
import sys
import numpy as np

SM = sys.argv[1] if len(sys.argv) > 1 else \
    "run_h3_kendrick2018/Output_Data/APH3D/APH3D_H3_DMBE/Output/Smatrix.txt"
if not SM.endswith(".txt"):
    SM = SM + "/Output_Data/APH3D/APH3D_H3_DMBE/Output/Smatrix.txt"

qn, S = [], None
with open(SM) as fh:
    started = False
    cur = []
    for line in fh:
        s = line.strip()
        if s.startswith("# energy"):
            if started:
                break
            started = True
        elif s.startswith("# QN") and "chanl" not in s:
            t = s.split()
            qn.append(tuple(int(x) for x in t[2:7]))
        elif s and not s.startswith("#") and started:
            re_, im_ = s.split()
            cur.append(complex(float(re_), float(im_)))
            if len(cur) == len(qn) ** 2:
                S = np.array(cur).reshape(len(qn), len(qn))
                break

idx = {q: i for i, q in enumerate(qn)}
n = len(qn)
print(f"Nch={n}")

a1 = idx[(1, 0, 4, 0, 0)]
a2 = idx[(2, 0, 4, 0, 0)]
a3 = idx[(3, 0, 4, 0, 0)]
print(f"incident rows: a1={a1+1}, a2={a2+1}, a3={a3+1} (1-based)")

# S3 covariance tests
print("\n[1] diagonal-covariance: S[(A,v'j')<-aA] vs S[(1,v'j')<-a1]  (f_N copies)")
worst = 0.0
for q in qn:
    if q[0] != 1 or (q[2], q[3]) == (4, 0):
        continue
    v, j, l = q[2], q[3], q[4]
    f1 = idx.get((1, 0, v, j, l)); f2 = idx.get((2, 0, v, j, l)); f3 = idx.get((3, 0, v, j, l))
    if f2 is None or f3 is None:
        continue
    ref = S[f1, a1]
    d2 = abs(S[f2, a2] - ref); d3 = abs(S[f3, a3] - ref)
    scale = max(abs(ref), 1e-30)
    worst = max(worst, d2 / scale, d3 / scale)
    if d2 / scale > 0.05 or d3 / scale > 0.05:
        print(f"   (v={v},j={j}): |d2|/|ref|={d2/scale:.3f}  |d3|/|ref|={d3/scale:.3f}"
              f"   ref={ref:.3e} S22={S[f2,a2]:.3e} S33={S[f3,a3]:.3e}")
print(f"   worst relative violation = {worst:.4f}")

print("\n[2] transposition structure: S[(2,v'j')<-a1] vs S[(3,v'j')<-a1] vs S[(1,v'j')<-a2]")
worst_x = 0.0
for q in qn:
    if q[0] != 1 or (q[2], q[3]) == (4, 0):
        continue
    v, j, l = q[2], q[3], q[4]
    f1 = idx.get((1, 0, v, j, l)); f2 = idx.get((2, 0, v, j, l)); f3 = idx.get((3, 0, v, j, l))
    if f2 is None or f3 is None:
        continue
    x23 = S[f2, a1]; x32 = S[f3, a1]; x12 = S[f1, a2]; x13 = S[f1, a3]
    scale = max(abs(x23), 1e-30)
    viol = max(abs(x23 - x32), abs(x12 - x32), abs(x13 - x23)) / scale
    worst_x = max(worst_x, viol)
print(f"   worst relative violation = {worst_x:.4f}")

# [3] does the (2,4,0) row look like the (1,4,0) row at all?
row1 = S[a1, :]; row2 = S[a2, :]
print(f"\n[3] |S[a1,a1]|={abs(row1[a1]):.6f}  |S[a2,a2]|={abs(row2[a2]):.6f}"
      f"  |S[a1,a2]|={abs(row1[a2]):.3e}")

# [4] the physical S3-symmetrized amplitude f_N + sqrt(2) f_R vs f_N +- f_R
print("\n[4] total rates with different symmetrization formulas (even j', excl (4,0))")
def prefactor(e_coll_Ha):
    AMU_TO_AU = 1822.888486209; M_H = 1.007825032071
    BOHR_CM = 0.529177210903e-8; AUTIME_S = 2.4188843265857e-17
    mu = (2.0/3.0) * M_H * AMU_TO_AU
    return (np.pi / np.sqrt(2.0 * mu**3 * e_coll_Ha)) * BOHR_CM**3 / AUTIME_S

pref = prefactor(3.1668114e-12)
s_minus = s_plus = s_sym3 = 0.0
incoh = 0.0
for q in qn:
    if q[0] != 1 or q[3] % 2 or (q[2], q[3]) == (4, 0):
        continue
    v, j, l = q[2], q[3], q[4]
    f1 = idx.get((1, 0, v, j, l)); f2 = idx.get((2, 0, v, j, l)); f3 = idx.get((3, 0, v, j, l))
    if f2 is None or f3 is None:
        continue
    fN = S[f1, a1]; fR = (S[f2, a1] + S[f3, a1]) / np.sqrt(2.0)
    s_minus += abs(fN - fR) ** 2
    s_plus += abs(fN + fR) ** 2
    s_sym3 += abs(fN + np.sqrt(2.0) * fR) ** 2
    incoh += abs(fN) ** 2 + abs(fR) ** 2
print(f"   |f_N - f_R|^2   -> k = {pref*s_minus:.4e} cm3/s")
print(f"   |f_N + f_R|^2   -> k = {pref*s_plus:.4e} cm3/s")
print(f"   |f_N+sqrt2 f_R|^2 -> k = {pref*s_sym3:.4e} cm3/s   (S3-symmetrized form)")
print(f"   incoherent sum  -> k = {pref*incoh:.4e} cm3/s")
print("   literature (PRL Table I, 1 uK): GP total = 2.48e-11, NGP total = 2.02e-12")
