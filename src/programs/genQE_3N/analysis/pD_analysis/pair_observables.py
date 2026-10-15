#!/usr/bin/env python3
"""Three more (p,n)-pair observables, for every source and potential.

  1  prel      |(p_proton - p_neutron)/2|        relative momentum of the pair
  2  ang_pn    angle(p_proton, p_neutron)        opening angle of the two nucleons
  3  Mmiss_pd  missing mass of A(e,e'pd)X        treating the pair as a deuteron

Target is 12C (QEGenerator_3N.cc:49-50, SetTargetNucleus(12,6)); ejecting a ppn
triplet leaves A-3 = 9Be. Verified from the samples: E_lead+E2+E3-nu = 2.78172
GeV against M(12C)-M(9Be) = 2.78214 (the 0.4 MeV is E*). The 3N density table is
a 3He density used as the universal 3N SRC density, exactly as GCF uses the
deuteron for 2N pairs -- the nucleus is 12C, not 3He.

Observable 3, following the definition asked for:
    p_miss_pd = p_lead + (p_proton + p_neutron) - q
    E_miss_pd = nu + M(12C) - E_lead - E_d ,  E_d = sqrt(M_d^2 + |p_pair|^2)
    M_miss_pd = sqrt(E_miss_pd^2 - |p_miss_pd|^2)
The undetected system is the 9Be residual, so M_miss_pd has a sharp lower edge
at M(9Be) = 8.3928 GeV. We plot M_miss_pd - M(9Be): that excess is the energy
cost of treating the free (p,n) pair as a bound deuteron, i.e. how deuteron-like
the pair is -- 0 for a perfect deuteron. It is the direct pD observable.
Only one pairing makes sense here (the pair is the two recoils), so unlike
observables 1-2 there is no struck/recoil choice.

Observables 1-2 are computed for both choices of the odd proton, and overlaid
with the same quantity straight from the density table (generator support).

All curves are normalised to unit AREA over the plotted range (h/sum(h)/binwidth,
weights = event weights), so the y axis is a probability density.  Entries
outside the plotted range are excluded from the normalisation; the summary file
reports that fraction per curve.

Writes pair_relmom.png, pair_opening_angle.png, missing_mass_pd.png and
pair_observables_summary.txt .
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
from opening_angle_from_events import EV, KF, EB, P, N, SAMPLES, ang, kinematics

HERE = os.path.dirname(os.path.abspath(__file__))
MP, MN, MD, MC12, MBE9 = 0.93827, 0.93957, 1.87561, 11.17493, 8.39279
B_PREL = np.linspace(0, 1.2, 49)
B_ANG = np.linspace(0, 180, 37)
B_MM = np.linspace(0, 1.5, 46)


def _fill(p1, rp, rn, lead, q, nu, w):
    """the three observables for one choice of pairing"""
    pair = rp + rn
    pm = np.linalg.norm(lead + pair - q, axis=1)
    Ed = np.sqrt(MD ** 2 + np.sum(pair ** 2, axis=1))
    El = np.sqrt(MP ** 2 + np.sum(lead ** 2, axis=1))
    Em = nu + MC12 - El - Ed
    return dict(w=w,
                prel=np.linalg.norm(0.5 * (rp - rn), axis=1),
                ang_pn=ang(rp, rn),
                mm=np.sqrt(np.maximum(Em ** 2 - pm ** 2, 0.0)) - MBE9)


def from_pwia(path):
    br = ["weight", "pe", "p1_gen", "p2_gen", "p3_gen", "N1_type", "N2_type", "N3_type"]
    o = {k: [] for k in ("w", "prel1", "ang1", "mm", "prel2", "ang2")}
    for a in uproot.iterate(f"{path}:genT", br, library="np", step_size=1_000_000):
        q, qm, Q2, the = kinematics(a["pe"])
        nu = EB - np.linalg.norm(a["pe"], axis=1)
        lead = a["p1_gen"]
        p1 = lead - q
        n2p = a["N2_type"] == P
        rp = np.where(n2p[:, None], a["p2_gen"], a["p3_gen"])      # recoil proton
        rn = np.where(n2p[:, None], a["p3_gen"], a["p2_gen"])      # recoil neutron
        ok = ((a["N1_type"] == P) & (a["N2_type"] != a["N3_type"])
              & np.isfinite(a["weight"]) & (a["weight"] > 0)
              & (the > 8) & (the < 45) & (Q2 > 1))
        # case 1: odd = struck proton  -> pair = (recoil p, n)   [also the (e,e'pd) pairing]
        c1 = _fill(p1[ok], rp[ok], rn[ok], lead[ok], q[ok], nu[ok], a["weight"][ok])
        # case 2: odd = recoil proton  -> pair = (struck p, n)
        o["w"].append(c1["w"]); o["prel1"].append(c1["prel"])
        o["ang1"].append(c1["ang_pn"]); o["mm"].append(c1["mm"])
        o["prel2"].append(np.linalg.norm(0.5 * (p1[ok] - rn[ok]), axis=1))
        o["ang2"].append(ang(p1[ok], rn[ok]))
    return {k: np.concatenate(v) for k, v in o.items()}


def from_fsi(path):
    br = ["weight", "pe", "part_pdg", "part_px", "part_py", "part_pz"]
    o = {k: [] for k in ("w", "prel1", "ang1", "mm", "prel2", "ang2")}
    for a in uproot.iterate(f"{path}:genT", br, library="ak", step_size=400_000):
        pe = ak.to_numpy(a.pe)
        q, qm, Q2, the = kinematics(pe)
        nu = EB - np.linalg.norm(pe, axis=1)
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
        c1 = _fill(p1[ok], rp[ok], rn[ok], lead[ok], q[ok], nu[ok], w[ok])
        o["w"].append(c1["w"]); o["prel1"].append(c1["prel"])
        o["ang1"].append(c1["ang_pn"]); o["mm"].append(c1["mm"])
        o["prel2"].append(np.linalg.norm(0.5 * (p1[ok] - rn[ok]), axis=1))
        o["ang2"].append(ang(p1[ok], rn[ok]))
    return {k: np.concatenate(v) for k, v in o.items()}


def shape(v, w, bins):
    """Unit-AREA normalisation: h / sum(h) / bin_width, so the curve integrates
    to 1 over the plotted range.  Weights are the event weights, so sum(h) is
    sum(w), not a count.  np.histogram drops out-of-range entries from both the
    numerator and sum(h), so the normalisation is over the plotted window only --
    the returned `over` is the weight fraction that falls outside it.  The
    returned mean is over the FULL sample, so it is not biased by that
    truncation (it can therefore sit slightly off the drawn curve's own mean)."""
    h, _ = np.histogram(v, bins=bins, weights=w)
    s = max(h.sum(), 1e-300)
    tot = max(w.sum(), 1e-300)
    over = 1.0 - s / tot
    # mean over the FULL sample, not the window: with a few % of the FSI weight
    # beyond the axis a windowed mean is biased low
    mean = float(np.average(v, weights=w))
    return h / s / np.diff(bins), mean, over


print("reading ...")
D, out = {}, []
for lab, u, _ in SAMPLES:
    D[("PWIA", u)] = from_pwia(f"{EV}/pwia_{u}.root")
    D[("FSI", u)] = from_fsi(f"{EV}/fsi_{u}.root")
    print(f"  {lab}")

DENS = {}
for lab, u, fn in SAMPLES:
    DENS[u] = J.pair_histograms(np.loadtxt(os.path.join(J.ARR, fn)).reshape(J.NT, J.NK, J.NK),
                                B_PREL, B_ANG)

PANELS = [("prel", "prel1", "prel2", B_PREL,
           r"$|(\vec p_p-\vec p_n)/2|$  [GeV/c]", "pair_relmom.png",
           r"relative momentum of the $(p,n)$ pair"),
          ("ang_pn", "ang1", "ang2", B_ANG,
           r"$\theta(\vec p_p,\vec p_n)$  [deg]", "pair_opening_angle.png",
           r"opening angle between the two nucleons of the $(p,n)$ pair")]

for dkey, k1, k2, bins, xl, png, title in PANELS:
    fig, axes = plt.subplots(2, 4, figsize=(19, 8.6), sharex=True)
    for r, src in enumerate(("PWIA", "FSI")):
        for c, (lab, u, _) in enumerate(SAMPLES):
            ax, d = axes[r][c], D[(src, u)]
            ax.stairs(DENS[u][dkey + "_sup"]["h"], bins, color="0.55", lw=1.4, ls=":",
                      label="wavefunction")
            for key, col, ls, nm in ((k1, "#2b6cb0", "-", r"odd = struck $p$"),
                                     (k2, "#dd8800", "--", r"odd = recoil $p$")):
                h, mu, over = shape(d[key], d["w"], bins)
                ax.stairs(h, bins, color=col, lw=2.0, ls=ls, label=f"{nm}  $\\langle\\cdot\\rangle$={mu:.3g}")
                out.append(f"{dkey:7s} {src:4s} {lab:11s} {nm:18s} mean = {mu:8.4f}"
                           f"   weight outside plotted range = {100*over:5.2f} %")
            ax.set_title(f"{src} — {lab}", fontsize=11)
            ax.set_xlim(bins[0], bins[-1]); ax.set_ylim(bottom=0)
            if r == 1: ax.set_xlabel(xl)
            if c == 0: ax.set_ylabel("normalised")
            ax.legend(frameon=False, fontsize=8)
    fig.suptitle(r"$^{12}$C(e,e'ppn): " + title +
                 "\ndotted grey = same quantity from the density table, on the generator support",
                 fontsize=12.5)
    fig.tight_layout(rect=(0, 0, 1, 0.91))
    fig.savefig(os.path.join(HERE, png), dpi=150); plt.close(fig)
    print("saved", png)
    out.append(f"{dkey:7s} wavefunction (gen support) mean = "
               + ", ".join(f"{l} {DENS[u][dkey+'_sup']['mean']:.4f}" for l, u, _ in SAMPLES))

# ---- observable 3: missing mass of (e,e'pd)
fig, axes = plt.subplots(1, 2, figsize=(13, 5.2))
for ax, src in zip(axes, ("PWIA", "FSI")):
    for (lab, u, _), col in zip(SAMPLES, ("#2b6cb0", "#c53030", "#2f855a", "#6b46c1")):
        d = D[(src, u)]
        h, mu, over = shape(d["mm"], d["w"], B_MM)
        ax.stairs(h, B_MM, color=col, lw=2.0, label=f"{lab}   $\\langle\\cdot\\rangle$={mu:.3f} GeV")
        out.append(f"Mmiss_pd {src:4s} {lab:11s} <M_miss - M(9Be)> = {mu:7.4f} GeV   "
                   f"frac < 50 MeV = {np.average(d['mm'] < 0.05, weights=d['w']):.4f}"
                   f"   weight above {B_MM[-1]:.1f} GeV = {100*over:5.2f} %")
    ax.axvline(0, color="0.3", lw=1.2, ls="--")
    ax.text(0.012, ax.get_ylim()[1] * .95, r"$M(^{9}$Be$)$ : perfect deuteron",
            fontsize=9, color="0.3", va="top")
    ax.set_title(src, fontsize=12)
    ax.set_xlabel(r"$M_{miss}(e,e'pd) - M(^{9}\mathrm{Be})$   [GeV]")
    ax.set_ylabel("normalised"); ax.set_xlim(0, 1.5); ax.set_ylim(bottom=0)
    ax.legend(frameon=False, fontsize=9)
fig.suptitle(r"$^{12}$C(e,e'pd)X missing mass:  $p_{miss}=p_{lead}+(p_p+p_n)-q$,  "
             r"$E_{miss}=\nu+M(^{12}\mathrm{C})-E_{lead}-E_d$"
             "\nthe undetected system is $^{9}$Be, so the excess above $M(^{9}$Be$)$ is how far the "
             r"free $(p,n)$ pair is from a bound deuteron", fontsize=12)
fig.tight_layout(rect=(0, 0, 1, 0.88))
fig.savefig(os.path.join(HERE, "missing_mass_pd.png"), dpi=150)
print("saved missing_mass_pd.png")

open(os.path.join(HERE, "pair_observables_summary.txt"), "w").write("\n".join(out) + "\n")
print("\n".join(out))
