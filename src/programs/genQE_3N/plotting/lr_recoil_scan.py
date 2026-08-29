#!/usr/bin/env python3
"""L/R region fractions of the (theta12, theta23) plane vs. recoil-momentum
threshold, for detection-based `events` trees (make_events_from_genT.py /
SRC_analysis_2N schema).

Applies the unified paper (e,e'pp) selection of plot_fig5_theta_heatmap_pair.py
(exactly 2 detected protons, lead cuts, pmiss window, xB, reconstructed third
above kF) and then scans an additional |p_recoil| > p_rec requirement. FSI
ejecta that fake the second proton are soft, true SRC partners are not, so
this threshold controls how much of the 2N background is "pn pair + fake
proton" (left-skewed) versus "pp pair" (right-skewed, back-to-back).

Usage: lr_recoil_scan.py LABEL=path [LABEL=path ...] [--kF 0.22] [--prec 0.22,0.3,...]
"""
import argparse
import numpy as np
import uproot


def angdeg(u, v):
    return np.degrees(np.arccos(np.clip(
        np.sum(u * v, axis=1)
        / (np.linalg.norm(u, axis=1) * np.linalg.norm(v, axis=1) + 1e-30), -1, 1)))


def regions(t12, t23, A=135., K=4.):
    Ar = np.radians(A)
    B = np.degrees(np.arctan(np.sin(Ar) / (K + np.cos(Ar)))) / (180. - A)
    R = ((t23 <= B * (t12 - 180.) + 180.) & (t23 >= (1. / B) * (t12 - 180.) + 180.)
         & (t23 >= -(t12 - A) + 180. - (180. - A) * B))
    m = 360. - t12 - t23
    L = ((t23 <= B * (m - 180.) + 180.) & (t23 >= (1. / B) * (m - 180.) + 180.)
         & (t23 >= -(m - A) + 180. - (180. - A) * B))
    return L, R


def scan(path, label, kf, precs, theta_e_max=45., Q2_min=1., xB_max=1.2,
         theta_pq_max=8., lead_over_q_min=0.75, pmiss_lo=0.25, pmiss_hi=0.9):
    t = uproot.open(path)['events']
    a = t.arrays(['weight', 'lead', 'recoil', 'q', 'Q2', 'xB', 'pmiss',
                  'scattering_angle', 'nProtonAboveKF'], library='np')
    w = a['weight']; lead = a['lead']; rec = a['recoil']; q = a['q']
    lm = np.linalg.norm(lead[:, :3], axis=1); qm = np.linalg.norm(q[:, :3], axis=1)
    ang = angdeg(lead[:, :3], q[:, :3])
    p1 = lead[:, :3] - q[:, :3]; p3 = -(p1 + rec[:, :3]); p3m = np.linalg.norm(p3, axis=1)
    rm = np.linalg.norm(rec[:, :3], axis=1)
    m0 = (np.isfinite(w) & (w > 0) & (a['nProtonAboveKF'] == 2)
          & (a['scattering_angle'] < theta_e_max) & (a['Q2'] > Q2_min)
          & (a['xB'] < xB_max) & (ang < theta_pq_max) & (lm / qm > lead_over_q_min)
          & (a['pmiss'] > pmiss_lo) & (a['pmiss'] < pmiss_hi) & (p3m > kf))
    rows = []
    for prec in precs:
        m = m0 & (rm > prec)
        t12 = angdeg(p1[m], rec[m, :3]); t23 = angdeg(rec[m, :3], p3[m])
        L, R = regions(t12, t23); ww = w[m]; wt = ww.sum() + 1e-30
        fL = ww[L].sum() / wt; fR = ww[R].sum() / wt
        rows.append((prec, int(m.sum()), fL, fR, fL / max(fR, 1e-30)))
        print(f"  {label:8s} p_rec>{prec:.2f}: pass={m.sum():8d}  "
              f"f_L={100*fL:5.2f}%  f_R={100*fR:5.2f}%  L/R={fL/max(fR,1e-30):6.3f}")
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('samples', nargs='+', help='LABEL=path/to/events.root')
    ap.add_argument('--kF', type=float, default=0.22)
    ap.add_argument('--prec', default='0.22,0.30,0.35,0.40,0.50')
    args = ap.parse_args()
    precs = [float(x) for x in args.prec.split(',')]
    for s in args.samples:
        label, path = s.split('=', 1)
        print(f"{label}: {path}")
        scan(path, label, args.kF, precs)


if __name__ == '__main__':
    main()
