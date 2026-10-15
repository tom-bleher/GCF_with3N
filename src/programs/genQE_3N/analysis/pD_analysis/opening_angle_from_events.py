#!/usr/bin/env python3
"""Opening angle  theta_p,(pn)  from generated events: PWIA, and the full
generator with FSI.  Companion to opening_angle_from_density.py, same observable.

Target: 12C (QEGenerator_3N.cc:49-50); ejecting a ppn triplet leaves 9Be.  The
density table is a 3He 3N density used as the universal 3N SRC density, the way
GCF uses the deuteron for 2N pairs.

The ejected 3N triplet is ppn, so picking one PROTON as the odd nucleon always leaves a (p,n)
pair.  The pair is not identical-particle, so its relative momentum has a fixed
sign; we use (proton - neutron)/2 everywhere, as in the density script.

  case 1   odd nucleon = the STRUCK proton
           a = p_miss                      b = (p_recoil_p - p_n)/2
  case 2   odd nucleon = the RECOIL proton
           a = p_recoil_p                  b = (p_struck_p - p_n)/2
  theta = angle(a, b)

PWIA  -- generated (initial) momenta, exactly as they would be measured.
p_miss = p1_gen - q is the initial momentum of the struck nucleon (p1_gen is the
OUTGOING lead, so q must be removed).  Events are required to have
N1_type = proton.  No CM correction is applied: these samples were generated
with use_CM on (sigCM = 0.15 GeV/c per component, QEGenerator_3N.cc:26 -- a 12C
2N-pair width, so it is spurious for 3He, where there is no A-3 remnant to carry
momentum).  Subtracting P_tot/3 moves <theta> by only 0.3 deg, and it cannot be
done on the FSI side at all (P_tot is not measurable there), so both rows use
the raw, measurable momenta.

FSI   -- the measured final state, from the post-cascade particle list:
  detected    = final-state p or n with |p| > kF = 0.25 GeV/c
  lead        = detected nucleon closest in angle to q; required to be a proton
  recoil p/n  = highest-|p| detected proton / neutron among the rest
  p_miss      = p_lead - q
No CM removal is possible here (P_tot is not measurable and the cascade
redistributes it), so the FSI panels use the raw momenta.

Cuts on every sample: 8 < theta_e < 45 deg, Q2 > 1 GeV^2.  Produced twice, with
and without |a| > 0.25 GeV/c on the odd proton (the row-2 cut of the density
figure).

Writes opening_angle_events[_pcut].png and opening_angle_events_summary.txt .
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
import pd_jacobi as J

EV = "/home/ilayws/Ilay-Generators/src/programs/genQE_3N/events/wf_test"
ARR = "/home/ilayws/Ilay-Generators/src/programs/genQE_3N/array"
HERE = os.path.dirname(os.path.abspath(__file__))
KF, MN, EB, HBARC = 0.25, 0.93892, 6.0, 0.1973269804
P, N = 2212, 2112
NT, NK = 61, 101
BINS = np.linspace(0, 180, 37)        # 5 deg; the FSI samples are weight-limited
SAMPLES = [("AV18", "u1", "he3_3Ndensity_AV18.dnst"),
           ("CD-Bonn", "u4", "he3_3Ndensity_CDBonn.dnst"),
           ("pD AV18", "u6", "he3_3Ndensity_pD_AV18.dnst"),
           ("pD CD-Bonn", "u9", "he3_3Ndensity_pD_CDBonn.dnst")]


def ang(a, b):
    c = np.sum(a * b, 1) / np.maximum(np.linalg.norm(a, axis=1) * np.linalg.norm(b, axis=1), 1e-300)
    return np.degrees(np.arccos(np.clip(c, -1, 1)))


def kinematics(pe):
    q = np.column_stack([-pe[:, 0], -pe[:, 1], EB - pe[:, 2]])
    qm = np.linalg.norm(q, axis=1)
    nu = EB - np.linalg.norm(pe, axis=1)
    Q2 = qm ** 2 - nu ** 2
    the = np.degrees(np.arctan2(np.hypot(pe[:, 0], pe[:, 1]), pe[:, 2]))
    return q, qm, Q2, the


def from_pwia(path):
    """generated initial momenta, struck proton, CM removed"""
    br = ["weight", "pe", "p1_gen", "p2_gen", "p3_gen", "N1_type", "N2_type", "N3_type"]
    o = {k: [] for k in ("w", "t1", "t2", "a1", "a2")}
    for a in uproot.iterate(f"{path}:genT", br, library="np", step_size=1_000_000):
        q, qm, Q2, the = kinematics(a["pe"])
        p1 = a["p1_gen"] - q                       # initial struck-nucleon momentum
        p2, p3 = a["p2_gen"], a["p3_gen"]
        n2p = a["N2_type"] == P
        rp = np.where(n2p[:, None], p2, p3)        # the recoil PROTON
        rn = np.where(n2p[:, None], p3, p2)        # the recoil NEUTRON
        ok = ((a["N1_type"] == P) & (a["N2_type"] != a["N3_type"])
              & np.isfinite(a["weight"]) & (a["weight"] > 0)
              & (the > 8) & (the < 45) & (Q2 > 1))
        o["w"].append(a["weight"][ok])
        o["t1"].append(ang(p1[ok], 0.5 * (rp[ok] - rn[ok])))
        o["t2"].append(ang(rp[ok], 0.5 * (p1[ok] - rn[ok])))
        o["a1"].append(np.linalg.norm(p1[ok], axis=1))
        o["a2"].append(np.linalg.norm(rp[ok], axis=1))
    return {k: np.concatenate(v) for k, v in o.items()}


def from_fsi(path):
    """measured final state after the cascade"""
    br = ["weight", "pe", "part_pdg", "part_px", "part_py", "part_pz"]
    o = {k: [] for k in ("w", "t1", "t2", "a1", "a2")}
    for a in uproot.iterate(f"{path}:genT", br, library="ak", step_size=400_000):
        pe = ak.to_numpy(a.pe)
        q, qm, Q2, the = kinematics(pe)
        pm = np.sqrt(a.part_px ** 2 + a.part_py ** 2 + a.part_pz ** 2)
        isp, isn = a.part_pdg == P, a.part_pdg == N
        det = (isp | isn) & (pm > KF)
        # lead = detected nucleon closest in angle to q
        cosq = ((a.part_px * q[:, 0] + a.part_py * q[:, 1] + a.part_pz * q[:, 2])
                / np.maximum(pm * qm, 1e-300))
        il = ak.argmax(ak.where(det, cosq, -2.0), axis=1, keepdims=True)
        lead = np.column_stack([ak.to_numpy(ak.firsts(v[il])) for v in
                                (a.part_px, a.part_py, a.part_pz)])
        lead_is_p = ak.to_numpy(ak.fill_none(ak.firsts(a.part_pdg[il]), 0)) == P
        rest = det & (ak.local_index(a.part_pdg) != ak.firsts(il))
        got, vecs = [], []
        for sel in (rest & isp, rest & isn):       # recoil proton, recoil neutron
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
        o["t1"].append(ang(p1[ok], 0.5 * (rp[ok] - rn[ok])))
        o["t2"].append(ang(rp[ok], 0.5 * (p1[ok] - rn[ok])))
        o["a1"].append(np.linalg.norm(p1[ok], axis=1))
        o["a2"].append(np.linalg.norm(rp[ok], axis=1))
    return {k: np.concatenate(v) for k, v in o.items()}


_DCACHE = {}


def density_curve(fn, support=True):
    """Same observable straight from the table (choice A), for the reference
    curve.  support=True restricts to the generator's own support -- all three
    internal momenta above fkMinInternal = 0.25 GeV/c (QEGenerator_3N.cc:27)."""
    if fn not in _DCACHE:
        _DCACHE[fn] = J.histograms(np.loadtxt(os.path.join(ARR, fn)).reshape(NT, NK, NK), BINS)
    return _DCACHE[fn]["thA" + ("_sup" if support else "")]["h"]


def shape(t, w, m):
    """normalised density, its mean, and the effective statistics
    N_eff = (sum w)^2 / sum w^2 -- the weighted samples are far from their
    raw entry count, which is what makes the FSI panels noisy"""
    h, _ = np.histogram(t[m], bins=BINS, weights=w[m])
    s = max(h.sum(), 1e-300)
    ww = w[m]
    neff = float(ww.sum() ** 2 / max((ww ** 2).sum(), 1e-300))
    return h / s / np.diff(BINS), float(np.sum(0.5 * (BINS[1:] + BINS[:-1]) * h) / s), neff


if __name__ == "__main__":
    print("reading samples ...")
    D = {}
    for lab, u, _ in SAMPLES:
        D[("PWIA", u)] = from_pwia(f"{EV}/pwia_{u}.root"); print(f"  pwia_{u}: {len(D[('PWIA',u)]['w'])} events")
        D[("FSI", u)] = from_fsi(f"{EV}/fsi_{u}.root");    print(f"  fsi_{u} : {len(D[('FSI',u)]['w'])} events")

    out = []
    for support in (True,):
        fig, axes = plt.subplots(2, 4, figsize=(19, 8.6), sharex=True)
        for r, src in enumerate(("PWIA", "FSI")):
            for c, (lab, u, fn) in enumerate(SAMPLES):
                ax, d = axes[r][c], D[(src, u)]
                m1 = m2 = np.ones(len(d["w"]), bool)
                h1, mu1, n1 = shape(d["t1"], d["w"], m1)
                h2, mu2, n2 = shape(d["t2"], d["w"], m2)
                ax.stairs(density_curve(fn, support), BINS, color="0.55", lw=1.3, ls=":",
                          label="wavefunction")
                ax.stairs(h1, BINS, color="#2b6cb0", lw=2.2, label=r"odd = struck $p$")
                ax.stairs(h2, BINS, color="#dd8800", lw=1.6, ls="--", label=r"odd = recoil $p$")
                ax.set_title(f"{src} — {lab}", fontsize=11)
                ax.set_xlim(0, 180); ax.set_ylim(bottom=0)
                if r == 1: ax.set_xlabel(r"$\theta$ [deg]")
                if c == 0: ax.set_ylabel("normalised")
                if r == 0 and c == 0: ax.legend(frameon=False, fontsize=8)
                out.append(f"{src:4s} {lab:11s}  <theta> odd=struck {mu1:6.2f}  odd=recoil {mu2:6.2f}"
                           f"   |  raw {len(d['w']):8d} ev   N_eff {n1:9.0f}")
        fig.suptitle(r"$^3$He: angle between the odd proton and $(\vec p_p-\vec p_n)/2$ of the remaining pair"
                     "\ndotted grey = same observable from the density table, on the generator's own "
                     "support (all internal $|p_i| > 0.25$ GeV/c)", fontsize=13)
        fig.tight_layout(rect=(0, 0, 1, 0.92))
        png = os.path.join(HERE, "opening_angle_events.png")
        fig.savefig(png, dpi=150); plt.close(fig); print("saved", png)

    open(os.path.join(HERE, "opening_angle_events_summary.txt"), "w").write("\n".join(out) + "\n")
    print("\n".join(out))
