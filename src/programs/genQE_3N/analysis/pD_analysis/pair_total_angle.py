#!/usr/bin/env python3
"""Angle between the odd proton and the TOTAL momentum of the remaining (p,n)
pair:   a = p_odd proton ,  b = p_proton + p_neutron ,  theta = angle(a,b).

This carries no wavefunction information.  In the 3N rest frame
p1 + p2 + p3 = 0, so b = -a exactly and theta = 180 deg for every configuration,
whatever the potential.  Verified: with use_CM off (events/wf_test/
pwia_u1_noCM.root) the angle is 180.0000 deg in every event (min = max = 180).

What it does measure:
  PWIA  -- the motion of the 3N cluster inside the nucleus.  The target is 12C
           (SetTargetNucleus(12,6), QEGenerator_3N.cc:49-50), so the triplet has
           a CM momentum v_cm (sigCM = 0.15 GeV/c per component, the literature
           12C value) and recoils against the A-3 = 9Be remnant.  Then
           p1+p2+p3 = v_cm != 0 and b = v_cm - a, with the width set by
           |v_cm| / |a|.  This is physical; the no-CM sample is the artificial
           one.
  FSI   -- genuine physics: the measured final state does not sum to zero,
           because the cascade redistributes momentum, pushes nucleons below the
           detection threshold and knocks out extra ones.

Definitions follow opening_angle_from_events.py exactly; only b changes
(total momentum of the pair instead of its relative momentum).

Writes pair_total_angle.png and pair_total_angle_summary.txt .
"""
import os
import sys

import numpy as np
import awkward as ak
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from opening_angle_from_events import EV, KF, EB, P, N, SAMPLES, ang, kinematics

HERE = os.path.dirname(os.path.abspath(__file__))
BINS = np.linspace(0, 180, 73)        # 2.5 deg


def pwia(path):
    br = ["weight", "pe", "p1_gen", "p2_gen", "p3_gen", "N1_type", "N2_type", "N3_type"]
    o = {k: [] for k in ("w", "t1", "t2")}
    for a in uproot.iterate(f"{path}:genT", br, library="np", step_size=1_000_000):
        q, qm, Q2, the = kinematics(a["pe"])
        p1 = a["p1_gen"] - q
        p2, p3 = a["p2_gen"], a["p3_gen"]
        n2p = a["N2_type"] == P
        rp = np.where(n2p[:, None], p2, p3)          # recoil proton
        rn = np.where(n2p[:, None], p3, p2)          # recoil neutron
        ok = ((a["N1_type"] == P) & (a["N2_type"] != a["N3_type"])
              & np.isfinite(a["weight"]) & (a["weight"] > 0)
              & (the > 8) & (the < 45) & (Q2 > 1))
        o["w"].append(a["weight"][ok])
        o["t1"].append(ang(p1[ok], rp[ok] + rn[ok]))     # odd = struck proton
        o["t2"].append(ang(rp[ok], p1[ok] + rn[ok]))     # odd = recoil proton
    return {k: np.concatenate(v) for k, v in o.items()}


def fsi(path):
    br = ["weight", "pe", "part_pdg", "part_px", "part_py", "part_pz"]
    o = {k: [] for k in ("w", "t1", "t2")}
    for a in uproot.iterate(f"{path}:genT", br, library="ak", step_size=400_000):
        pe = ak.to_numpy(a.pe)
        q, qm, Q2, the = kinematics(pe)
        pm = np.sqrt(a.part_px ** 2 + a.part_py ** 2 + a.part_pz ** 2)
        isp, isn = a.part_pdg == P, a.part_pdg == N
        det = (isp | isn) & (pm > KF)
        cosq = ((a.part_px * q[:, 0] + a.part_py * q[:, 1] + a.part_pz * q[:, 2])
                / np.maximum(pm * qm, 1e-300))
        il = ak.argmax(ak.where(det, cosq, -2.0), axis=1, keepdims=True)
        lead = np.column_stack([ak.to_numpy(ak.firsts(v[il])) for v in
                                (a.part_px, a.part_py, a.part_pz)])
        lead_is_p = ak.to_numpy(ak.fill_none(ak.firsts(a.part_pdg[il]), 0)) == P
        rest = det & (ak.local_index(a.part_pdg) != ak.firsts(il))
        got, vecs = [], []
        for sel in (rest & isp, rest & isn):
            i = ak.argmax(ak.where(sel, pm, -1.0), axis=1, keepdims=True)
            got.append(ak.to_numpy(ak.sum(sel, axis=1)) > 0)
            vecs.append(np.nan_to_num(np.column_stack(
                [ak.to_numpy(ak.firsts(v[i])) for v in (a.part_px, a.part_py, a.part_pz)])))
        rp, rn = vecs
        w = ak.to_numpy(a.weight)
        p1 = lead - q
        ok = (lead_is_p & got[0] & got[1] & np.isfinite(w) & (w > 0)
              & (the > 8) & (the < 45) & (Q2 > 1))
        o["w"].append(w[ok])
        o["t1"].append(ang(p1[ok], rp[ok] + rn[ok]))
        o["t2"].append(ang(rp[ok], p1[ok] + rn[ok]))
    return {k: np.concatenate(v) for k, v in o.items()}


def shape(t, w):
    h, _ = np.histogram(t, bins=BINS, weights=w)
    s = max(h.sum(), 1e-300)
    return h / s / np.diff(BINS), float(np.sum(0.5 * (BINS[1:] + BINS[:-1]) * h) / s)


print("reading ...")
D = {}
for lab, u, _ in SAMPLES:
    D[("PWIA", u)] = pwia(f"{EV}/pwia_{u}.root")
    D[("FSI", u)] = fsi(f"{EV}/fsi_{u}.root")
    print(f"  {lab}")
D[("PWIA", "noCM")] = pwia(f"{EV}/pwia_u1_noCM.root")

out = []
fig, axes = plt.subplots(2, 4, figsize=(19, 8.6), sharex=True)
for r, src in enumerate(("PWIA", "FSI")):
    for c, (lab, u, _) in enumerate(SAMPLES):
        ax, d = axes[r][c], D[(src, u)]
        for key, col, ls, nm in (("t1", "#2b6cb0", "-", r"odd = struck $p$"),
                                 ("t2", "#dd8800", "--", r"odd = recoil $p$")):
            h, mu = shape(d[key], d["w"])
            ax.stairs(h, BINS, color=col, lw=2.0, ls=ls, label=f"{nm}   $\\langle\\theta\\rangle$={mu:.1f}°")
            out.append(f"{src:4s} {lab:11s} {nm:18s} <theta> = {mu:6.2f} deg")
        if src == "PWIA" and u == "u1":
            ax.axvline(180, color="#c53030", lw=2.5, alpha=.8)
            ax.text(176, ax.get_ylim()[1] * .55, "CM off (artificial):\nexactly 180° always",
                    ha="right", fontsize=8.5, color="#c53030")
        ax.set_title(f"{src} — {lab}", fontsize=11)
        ax.set_xlim(0, 180); ax.set_ylim(bottom=0)
        if r == 1: ax.set_xlabel(r"$\theta$ [deg]")
        if c == 0: ax.set_ylabel("normalised")
        ax.legend(frameon=False, fontsize=8, loc="upper left")

fig.suptitle(r"$^{12}$C: angle between the odd proton and the TOTAL momentum of the remaining $(p,n)$ pair"
             "\n180° exactly if the triplet were at rest — PWIA measures the 3N cluster's CM motion in "
             r"$^{12}$C ($\sigma_{CM}=0.15$), FSI adds the cascade on top", fontsize=12)
fig.tight_layout(rect=(0, 0, 1, 0.91))
png = os.path.join(HERE, "pair_total_angle.png")
fig.savefig(png, dpi=150)
# unbinned, so the number is not clipped to the last bin centre (178.75)
dn = D[("PWIA", "noCM")]
out.append(f"PWIA AV18 no-CM control (unbinned): <theta> = "
           f"{np.average(dn['t1'], weights=dn['w']):.4f} deg, "
           f"min {dn['t1'].min():.4f}, max {dn['t1'].max():.4f}  "
           f"-- exactly 180 by momentum conservation")
open(os.path.join(HERE, "pair_total_angle_summary.txt"), "w").write("\n".join(out) + "\n")
print("\n".join(out)); print("saved", png)
