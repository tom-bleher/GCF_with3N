#!/usr/bin/env python3
"""L/R region fractions with the DETECTION-BASED 3N definitions, from the
post-FSI particle list (part_*), for each wavefunction table.

Definitions (conventions from plotting/plot_detection_heatmaps.py):
  detected nucleon = final-state p or n with |p| > kF = 0.25 GeV/c
  lead             = detected nucleon closest in angle to q
  p1               = lead - q          (reconstructed initial struck momentum)
  p2               = highest-|p| other detected nucleon
  p3               = -(p1 + p2)        (reconstructed third, zero initial CM)
  theta12 = angle(p1,p2)   theta23 = angle(p2,p3)
  "=2"  : exactly 2 detected nucleons
  ">=2" : at least 2 detected (third still reconstructed, real one ignored)
Cuts: 8 < theta_e < 45 deg, Q2 > 1, xB < 1.2 (same as the truth-level table).
Errors from sum(w^2).
"""
import os, sys
import numpy as np, awkward as ak, uproot
sys.path.insert(0, "/home/ilayws/Ilay-Generators/src/programs/genQE_3N/plotting")
from plot_fig5_theta_heatmap_pair import in_region_L_array, in_region_R_array

EV = "/home/ilayws/Ilay-Generators/src/programs/genQE_3N/events/wf_test"
KF, MN, EB = 0.25, 0.93892, 6.0
LAB = {"u1": "AV18 (u=1)", "u4": "CD-Bonn (u=4)", "u6": "pD AV18 (u=6)", "u9": "pD CD-Bonn (u=9)"}


def ang(a, b):
    return np.degrees(np.arccos(np.clip(np.sum(a*b, axis=1) /
                     (np.linalg.norm(a, axis=1)*np.linalg.norm(b, axis=1) + 1e-30), -1, 1)))


def analyse(path):
    t = uproot.open(path)["genT"]
    a = t.arrays(["weight", "pe", "part_pdg", "part_px", "part_py", "part_pz"], library="ak")
    w = ak.to_numpy(a.weight); pe = ak.to_numpy(a.pe)
    q = np.column_stack([-pe[:, 0], -pe[:, 1], EB - pe[:, 2]])
    qm = np.linalg.norm(q, axis=1); nu = EB - np.linalg.norm(pe, axis=1)
    Q2 = qm**2 - nu**2; xB = Q2/(2*MN*nu + 1e-30)
    th_e = np.degrees(np.arctan2(np.hypot(pe[:, 0], pe[:, 1]), pe[:, 2]))

    isN = (a.part_pdg == 2212) | (a.part_pdg == 2112)
    pm = np.sqrt(a.part_px**2 + a.part_py**2 + a.part_pz**2)
    det = isN & (pm > KF)
    ndet = ak.to_numpy(ak.sum(det, axis=1))

    # lead = detected nucleon closest in angle to q
    px, py, pz = (ak.where(det, v, 0.0) for v in (a.part_px, a.part_py, a.part_pz))
    pmag = ak.where(det, pm, 1e-30)
    cosq = (px*q[:, 0] + py*q[:, 1] + pz*q[:, 2]) / (pmag*qm)
    cosq = ak.where(det, cosq, -2.0)
    il = ak.argmax(cosq, axis=1, keepdims=True)
    lead = np.column_stack([ak.to_numpy(ak.firsts(v[il])) for v in (a.part_px, a.part_py, a.part_pz)])
    # p2 = highest-|p| other detected nucleon
    idx = ak.local_index(a.part_pdg)
    other = det & (idx != ak.firsts(il))
    pm2 = ak.where(other, pm, -1.0)
    i2 = ak.argmax(pm2, axis=1, keepdims=True)
    has2 = ak.to_numpy(ak.sum(other, axis=1)) > 0
    p2 = np.column_stack([ak.to_numpy(ak.firsts(v[i2])) for v in (a.part_px, a.part_py, a.part_pz)])
    p2 = np.nan_to_num(p2)

    p1 = lead - q
    p3 = -(p1 + p2)
    t12, t23 = ang(p1, p2), ang(p2, p3)
    base = (np.isfinite(w) & (w > 0) & (th_e > 8.) & (th_e < 45.) & (Q2 > 1.)
            & (xB < 1.2) & has2 & (ndet >= 2))
    return dict(w=w, t12=t12, t23=t23, base=base, ndet=ndet)


print(f"{'sample':22s} {'def':>5s} {'events':>9s} {'L%':>7s} {'R%':>7s} {'L/R':>8s} {'+-':>7s}")
for mode in ("fsi", "pwia"):
    for tag in ("u1", "u4", "u6", "u9"):
        p = f"{EV}/{mode}_{tag}.root"
        if not os.path.exists(p):
            continue
        d = analyse(p)
        for name, sel in (("=2", d["base"] & (d["ndet"] == 2)), (">=2", d["base"])):
            L = sel & in_region_L_array(d["t12"], d["t23"])
            R = sel & in_region_R_array(d["t12"], d["t23"])
            w = d["w"]; tot = w[sel].sum()
            sL, sR = w[L].sum(), w[R].sum()
            if sL <= 0 or sR <= 0:
                print(f"{mode+' '+LAB[tag]:22s} {name:>5s}  (empty)"); continue
            r = sL/sR
            err = r*np.hypot(np.sqrt((w[L]**2).sum())/sL, np.sqrt((w[R]**2).sum())/sR)
            print(f"{mode+' '+LAB[tag]:22s} {name:>5s} {int(sel.sum()):9d} "
                  f"{100*sL/tot:7.2f} {100*sR/tot:7.2f} {r:8.3f} {err:7.3f}")
    print()
