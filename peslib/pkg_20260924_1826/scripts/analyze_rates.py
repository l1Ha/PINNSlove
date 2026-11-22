#!/usr/bin/env python3
"""Extract P_react(E_coll) for the incident channel of every run under SCRATCH_DIR
and convert to physical rate coefficients k = (pi/(mu*k_au)) * P * 6.1263e-9 cm^3/s.

Fixes vs the pkg_20260921_1407 version:
  * NPES=2 (GP/diabatic) extraction bug fixed: the old version summed wrong
    S-matrix columns (P came out 0 for D+HD). Now P_react = 1 - |S_ii|^2 with
    incident smat index ic = (i0-1)*NPES + 1, i0 = incident base-channel row in
    energy_diatomic.dat (same indexing as the Fortran ic_in=(j_min-1)*NPES+1).
  * Incident channel auto-selected per run from its own energy_diatomic.dat:
    highest-energy arrangement-1 channel with j=0,k=0 open at enrg
    (-> v=4 for the 2.161 eV runs, v=0 for the 0.3 eV production runs).
  * Reduced mass from the input.d mass line (H3 vs D+HD vs D+H2).
  * Also records the total CRP (ReactionProbabilities*.dat, summed over all
    open reactant channels) as a secondary dataset <name>_crptot.

Usage:  python3 analyze_rates.py [SCRATCH_DIR]
Output: <pkg>/analysis/rates.npz and rates_summary.txt

If <pkg>/precomputed/rates_precomputed.npz exists it is merged in (fresh
scratch results override precomputed entries with the same name), so the
comparison figures work before/while the slow GP runs finish.
"""
import numpy as np, os, sys, glob, re, math

PKG = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AU2CC = 6.1263e-9         # a.u. rate -> cm^3/s
EVK = 11604.5180          # eV -> K
AMU_AU = 1822.8885

MMASS = {1: 1.007825, 2: 2.014000, 4: 4.002600, 19: 18.99840,
         39: 38.9637064864, 40: 39.963998166, 41: 40.9618252579,
         85: 84.9117897379, 87: 86.9091805310}


def parse_input(path):
    txt = open(path).read()
    def g(n, dflt=None):
        m = re.search(rf'^\s*{n}\s*=\s*([-\d.eE+]+)', txt, re.M)
        return float(m.group(1)) if m else dflt
    mline = re.search(r'^\s*mass\s*=\s*([\d\s,]+)', txt, re.M)
    masses = [int(x) for x in mline.group(1).replace(',', ' ').split()] if mline else [1, 1, 1]
    return dict(npes=int(g('NPES', 1)), masses=masses,
                enrg=g('enrg', 2.1610001), dnrg=g('dnrg', 1e-7), nnrg=int(g('nnrg', 300)),
                jtot=int(g('jtot', 0)),
                iscale=int(g('iscale', 0)),
                ecol_min=g('ecol_min', 8.617333e-11),
                ecol_max=g('ecol_max', 8.617333e-3))


def mu_au(masses):
    mA = MMASS.get(masses[0], masses[0] * 1.0)
    mBC = MMASS.get(masses[1], masses[1] * 1.0) + MMASS.get(masses[2], masses[2] * 1.0)
    return (mA * mBC / (mA + mBC)) * AMU_AU


def rate_ccs(E_eV, P, mua):
    E = np.maximum(np.asarray(E_eV, float), 1e-30) / 27.211386
    kk = np.sqrt(2 * mua * E)
    return AU2CC * np.pi * np.asarray(P, float) / (mua * kk)


def find_incident(rundir, enrg):
    """Highest-energy arrangement-1 channel open at enrg; prefer j=0,k=0.
    Parity-odd blocks (jpar=-1) contain no j=0 -> fall back to the lowest j
    (server-validated: J0-odd -> v=3,j=1 threshold 1.743 eV).
    Returns (i0, E0, v0, j0, k0)."""
    fn = os.path.join(rundir, 'output', 'energy_diatomic.dat')
    cands = []
    with open(fn) as f:
        for line in f:
            m = re.match(r'\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(-?\d+)\s+(\S+)', line)
            if not m:
                continue
            idx, a, v, j, k = int(m.group(1)), int(m.group(2)), int(m.group(3)), int(m.group(4)), int(m.group(5))
            E = float(m.group(6))
            if a == 1 and E <= enrg * (1 + 1e-6) + 1e-12:
                cands.append((idx, E, v, j, k))
    if not cands:
        return (None, None, None, None, None)
    pref = [c for c in cands if c[3] == 0 and c[4] == 0]
    if not pref:
        jmin = min(c[3] for c in cands)
        pref = [c for c in cands if c[3] == jmin]
    return max(pref, key=lambda c: c[1])


def get_all_channels(rundir):
    fn = os.path.join(rundir, 'output', 'energy_diatomic.dat')
    chans = []
    if os.path.exists(fn):
        with open(fn) as f:
            for line in f:
                m = re.match(r'\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(-?\d+)\s+(\S+)', line)
                if m:
                    chans.append((int(m.group(1)), int(m.group(2)), int(m.group(3)), int(m.group(4)), int(m.group(5)), float(m.group(6))))
    return chans


def parse_real(rundir, enrg, dnrg, nnrg, v, j, k, e0, iscale=0, ecol_min=None, ecol_max=None, is_h3=True):
    """NPES=1: extracts unsymmetrized P_react and coherent symmetrized P_sym according to Kendrick 2018 Eq. 13."""
    out_raw = []
    out_sym_ngp = []
    out_sym_gp = []
    lab = (str(v), str(j), str(k))
    for iN in range(1, nnrg + 1):
        fn = f'{rundir}/output/S-open{iN}.dat'
        p_react = 0.0
        found = False
        states = {}
        try:
            with open(fn) as f:
                for line in f:
                    p = line.split()
                    if len(p) >= 11 and p[0] == '1' and tuple(p[1:4]) == lab:
                        found = True
                        arr, vf, jf = int(p[4]), int(p[5]), int(p[6])
                        sr, si = float(p[8]), float(p[9])
                        s_val = complex(sr, si)
                        key = (vf, jf)
                        if key not in states:
                            states[key] = {}
                        states[key][arr] = s_val
                        if arr in (2, 3):
                            p_react += float(p[10])
        except OSError:
            pass
        if iscale == 1 and ecol_min and ecol_max:
            frac = (iN - 1) / max(1, nnrg - 1)
            e_coll = math.exp(math.log(ecol_min) + frac * (math.log(ecol_max) - math.log(ecol_min)))
        else:
            e_coll = enrg + (iN - 1) * dnrg - e0

        p_ngp = 0.0
        p_gp = 0.0
        if found:
            for (vf, jf), ad in states.items():
                if is_h3:
                    if jf % 2 == 0:  # Para-H2
                        if vf == v and jf == j:
                            continue  # exclude elastic
                        s_N = ad.get(1, 0j)
                        s_R = ad.get(2, 0j)
                        t_ngp = (1.0 / math.sqrt(2.0)) * (s_N - math.sqrt(2.0) * s_R)
                        t_gp = (1.0 / math.sqrt(2.0)) * (s_N + math.sqrt(2.0) * s_R)
                        p_ngp += abs(t_ngp)**2
                        p_gp += abs(t_gp)**2
                else:
                    # D + HD: permutation P_13 between Arr 1 and Arr 3
                    if vf == v and jf == j:
                        continue
                    s1 = ad.get(1, 0j)
                    s3 = ad.get(3, 0j)
                    s_even = (1.0 / math.sqrt(2.0)) * (s1 + s3)
                    s_odd = (1.0 / math.sqrt(2.0)) * (s1 - s3)
                    p_ngp += abs(s_even)**2 * (2.0 / 3.0)
                    p_gp += abs(s_odd)**2 * (1.0 / 3.0)
        else:
            p_react = np.nan
            p_ngp = np.nan
            p_gp = np.nan

        out_raw.append([e_coll, p_react if found else np.nan])
        out_sym_ngp.append([e_coll, p_ngp if found else np.nan])
        out_sym_gp.append([e_coll, p_gp if found else np.nan])

    return np.array(out_raw), np.array(out_sym_ngp), np.array(out_sym_gp)


def parse_complex(rundir, i0, nnrg, npes, e0):
    """NPES=2: P_react = sum_{jc in arr 2,3} |S_ic,jc|^2.
    Streams S-matrix-complex.dat; parses incident row ic and sums reaction channels."""
    fn = f'{rundir}/output/S-matrix-complex.dat'
    chans = get_all_channels(rundir)
    ic = (i0 - 1) * max(1, npes) + 1          # 1-based smat row of incident channel
    out = []
    curE = None; rowidx = 0; inb = False; pval = np.nan
    try:
        with open(fn) as f:
            for line in f:
                if line.startswith('E(eV)='):
                    if curE is not None:
                        out.append([curE - e0, pval])
                        if len(out) >= nnrg:
                            break
                    curE = float(line.split('=')[1].split()[0]); pval = np.nan
                    inb = False; rowidx = 0
                elif line.startswith('S-matrix complex block'):
                    inb = True; rowidx = 0
                elif inb:
                    rowidx += 1
                    if rowidx == ic:
                        p = [float(x) for x in line.split()]
                        if chans:
                            p_react = 0.0
                            for jc in range(1, len(p) // 2 + 1):
                                j_spatial = (jc - 1) // max(1, npes) + 1
                                if j_spatial <= len(chans):
                                    arr = chans[j_spatial - 1][1]
                                    if arr in (2, 3):
                                        re_s = p[2 * (jc - 1)]
                                        im_s = p[2 * (jc - 1) + 1]
                                        p_react += re_s**2 + im_s**2
                            pval = p_react
                        else:
                            col0 = 2 * (ic - 1)
                            s11 = complex(p[col0], p[col0 + 1])
                            pval = max(0.0, 1.0 - abs(s11) ** 2)
                        inb = False
    except OSError:
        pass
    if curE is not None and len(out) < nnrg:
        out.append([curE - e0, pval])
    return np.array(out)


def parse_crptot(rundir, npes, e0):
    for name in ('ReactionProbabilities-complex.dat', 'ReactionProbabilities.dat'):
        fn = os.path.join(rundir, 'output', name)
        if os.path.exists(fn) and os.path.getsize(fn) > 100:
            try:
                d = np.loadtxt(fn)
                if d.ndim == 2 and d.shape[1] >= 2:
                    return np.column_stack([d[:, 0] - e0, d[:, 1]])
            except Exception:
                pass
    return None


if __name__ == '__main__':
    scratch = sys.argv[1] if len(sys.argv) > 1 else os.path.join(PKG, 'runs')
    os.makedirs(f'{PKG}/analysis', exist_ok=True)
    stamp = os.popen('date "+%F %T"').read().strip()
    out = [f'analysis of runs under {scratch}  ({stamp})']
    R = {}
    for rdir in sorted(glob.glob(f'{scratch}/*/')):
        name = os.path.basename(rdir.rstrip('/'))
        idat = os.path.join(rdir, 'input.d')
        if not os.path.exists(idat):
            continue
        inp = parse_input(idat)
        i0, e0, v0, j0, k0 = find_incident(rdir, inp['enrg'])
        if i0 is None:
            out.append(f'{name}: no open arrangement-1 channel in energy_diatomic.dat - skipped')
            continue
        mua = mu_au(inp['masses'])
        is_h3 = (abs(inp['masses'][0] - inp['masses'][1]) < 0.1 and abs(inp['masses'][1] - inp['masses'][2]) < 0.1)
        if inp['npes'] == 1:
            v_raw, v_sngp, v_sgp = parse_real(rdir, inp['enrg'], inp['dnrg'], inp['nnrg'], v0, j0, k0, e0,
                                              inp.get('iscale', 0), inp.get('ecol_min'), inp.get('ecol_max'), is_h3=is_h3)
            v = v_raw
        else:
            v = parse_complex(rdir, i0, inp['nnrg'], inp['npes'], e0)
            v_raw = v; v_sngp = v; v_sgp = v
        vc = v[~np.isnan(v[:, 1])] if v.size else v
        ctot = parse_crptot(rdir, inp['npes'], e0)
        if ctot is not None:
            R[name + '_crptot'] = ctot
        if v.size == 0 or vc.shape[0] == 0:
            out.append(f'{name}: no incident-channel data extracted')
            continue
        R[name] = np.column_stack([vc[:, 0], vc[:, 1], rate_ccs(vc[:, 0], vc[:, 1], mua)])   # E_coll(eV), P, k(cm3/s)
        EcK = vc[:, 0] * EVK
        kk = rate_ccs(vc[:, 0], vc[:, 1], mua)
        out.append(f'{name}: npes={inp["npes"]} mua={mua:.1f}au n={len(vc)} '
                   f'E_c {EcK[0]:.4f}-{EcK[-1]:.2f} K  P(thr)={vc[0, 1]:.4f}  '
                   f'k({EcK[0]:.3f}K)={kk[0]:.3e}  k(10K)={np.interp(10, EcK, kk):.3e} cm3/s')

        # Symmetrized observables
        if inp['npes'] == 1:
            if is_h3:
                vc_sngp = v_sngp[~np.isnan(v_sngp[:, 1])]
                vc_sgp  = v_sgp[~np.isnan(v_sgp[:, 1])]
                if vc_sngp.size:
                    R[name + '_sym_ngp'] = np.column_stack([vc_sngp[:, 0], vc_sngp[:, 1], rate_ccs(vc_sngp[:, 0], vc_sngp[:, 1], mua)])
                if vc_sgp.size:
                    R[name + '_sym_gp']  = np.column_stack([vc_sgp[:, 0], vc_sgp[:, 1], rate_ccs(vc_sgp[:, 0], vc_sgp[:, 1], mua)])
            else:
                vc_even = v_sngp[~np.isnan(v_sngp[:, 1])]
                vc_odd  = v_sgp[~np.isnan(v_sgp[:, 1])]
                if vc_even.size:
                    R[name + '_even'] = np.column_stack([vc_even[:, 0], vc_even[:, 1], rate_ccs(vc_even[:, 0], vc_even[:, 1], mua)])
                if vc_odd.size:
                    R[name + '_odd']  = np.column_stack([vc_odd[:, 0], vc_odd[:, 1], rate_ccs(vc_odd[:, 0], vc_odd[:, 1], mua)])

    # Compute full partial wave sum if J=0, 1, 2 are present
    if 'h3_J0_even_ngp_log' in R and 'h3_J1_pneg_log' in R:
        v0 = R['h3_J0_even_ngp_log']
        v1 = R['h3_J1_pneg_log']
        v2 = R.get('h3_J2_ppos_log')
        n_pts = len(v0)
        k0 = v0[:, 2]
        k1 = 3.0 * v1[:, 2]
        k2 = np.zeros(n_pts)
        if v2 is not None:
            k2[-len(v2):] = 5.0 * v2[:, 2]
        ktot = k0 + k1 + k2
        R['h3_full_Jsum'] = np.column_stack([v0[:, 0], v0[:, 1], ktot])
        out.append(f'h3_full_Jsum: J=0+3*J=1+5*J=2, k(1uK)={ktot[0]:.3e}, k(10K)={np.interp(10, v0[:,0]*EVK, ktot):.3e}, k(100K)={ktot[-1]:.3e} cm3/s')

    # Strictly save only freshly calculated datasets from runs/ (NO precomputed archive merge)
    np.savez(f'{PKG}/analysis/rates.npz', **R)
    open(f'{PKG}/analysis/rates_summary.txt', 'w').write('\n'.join(out))
    print('\n'.join(out))
