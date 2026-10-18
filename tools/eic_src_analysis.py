#!/usr/bin/env python3
"""Acceptance analysis of run_eic_sim.sh runs: e + A -> e' N_lead N_rec (A-2), cf. PRC 105, 034001.

Two steps:

  extract RUN_DIR [RUN_DIR ...] [-j N]
      Reduce each run to RUN_DIR/analysis_events.npz, one row per simulated event:
      generator kinematics (rest frame), cross-section weight, and per nucleon (leading,
      recoil) the far-forward hit-level and reco-level detection flags.

  plots --setting LABEL=RUN_DIR[,RUN_DIR...] [--setting ...] -o OUTDIR
      Figures and a summary table. Q2 slices of one beam setting are listed together;
      each run's events are normalized by that run's own attempts, so slices add.

Event weight: w = (generator weight) / (attempts of its run), in nb, so sum(w) over a
selection is its cross section; yields for 10 fb^-1 are sigma[nb] * 1e7, times 0.5 for
nucleon transparency as in the paper (applied only where a figure says so).

Detection, hit level (the paper's definition: enough hits to reconstruct a track, or a
deposit in a calorimeter's volume):
  B0 tracker >= 3 stations, off-momentum (OMD) >= 2 stations, Roman pots (RP) >= 2 stations,
  counting hits made by the nucleon itself; ZDC: energy deposited by the nucleon or its
  shower (hits attributed through calorimeter contributions to the nucleon's descendants).
Detection, reco level:
  B0: a B0 CKF track associated to the nucleon; central/B0 charged: an associated
  ReconstructedChargedParticle; ZDC and central/B0 calorimeters: an associated cluster of
  the nucleon's descendants. RP and OMD have no reconstruction for ion beams in EICrecon,
  so they stay at hit level in the reco-level flags.
Angles are measured with respect to the ion beam (rest-frame kinematics boosted head-on,
before the crossing angle), as in the paper's Figs. 3 and 4.
"""
import argparse
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import awkward as ak
import numpy as np
import uproot

TRACKERS = {"B0TrackerHits": 3, "ForwardOffMTrackerHits": 2, "ForwardRomanPotHits": 2}
ZDC_HITS = ["HcalFarForwardZDCHits", "EcalFarForwardZDCHits"]
ZDC_CLUSTERS = ["HcalFarForwardZDCClusterAssociations", "EcalFarForwardZDCClusterAssociations"]
CENTRAL_CLUSTERS = ["LFHCALClusterAssociations", "EcalEndcapPClusterAssociations",
                    "HcalEndcapPInsertClusterAssociations", "B0ECalClusterAssociations",
                    "EcalBarrelClusterAssociations", "HcalBarrelClusterAssociations"]
STATION_GAP_MM = 30.0
P_IRF_MIN = 0.3  # GeV/c, the paper's cut against evaporation nucleons
ROLES = ("lead", "rec")


# ----------------------------------------------------------------------------- extract

def ancestors(st, pb, pe, pidx):
    """Index of the generator-level (status 1) ancestor of every MCParticle, -1 if none."""
    n = len(st)
    parent = np.full(n, -1)
    has = pe > pb
    parent[has] = pidx[pb[has]]
    anc = np.arange(n)
    for _ in range(64):
        step = (st[anc] != 1) & (parent[anc] >= 0)
        if not step.any():
            break
        anc[step] = parent[anc[step]]
    anc[st[anc] != 1] = -1
    return anc


def stations(z):
    if len(z) == 0:
        return 0
    z = np.sort(z)
    return 1 + int(np.sum(np.diff(z) > STATION_GAP_MM))


def extract_file(sim_path):
    """Per-event detection flags of the leading and recoil nucleon for one sim/reco chunk."""
    sim_path = Path(sim_path)
    reco_path = sim_path.with_name(sim_path.name.replace("sim_", "reco_").replace(".edm4hep.", ".edm4eic."))
    mc = ["MCParticles/MCParticles." + b for b in ("generatorStatus", "PDG", "parents_begin", "parents_end")]
    mc.append("_MCParticles_parents/_MCParticles_parents.index")
    br = ["EventHeader/EventHeader.eventNumber"] + mc
    for c in TRACKERS:
        br += [f"{c}/{c}.position.z", f"_{c}_particle/_{c}_particle.index"]
    for c in ZDC_HITS:
        br += [f"{c}Contributions/{c}Contributions.energy", f"_{c}Contributions_particle/_{c}Contributions_particle.index"]
    s = uproot.open(sim_path)["events"].arrays(br, library="ak")

    r = None
    if reco_path.exists():
        rt = uproot.open(reco_path)["events"]
        keys = set(rt.keys())
        rbr = ["EventHeader/EventHeader.eventNumber"] + mc
        assoc = ["B0TrackerCKFTrackAssociations", "ReconstructedChargedParticleAssociations"]
        assoc += [c for c in ZDC_CLUSTERS + CENTRAL_CLUSTERS if f"_{c}_sim/_{c}_sim.index" in keys]
        rbr += [f"_{c}_sim/_{c}_sim.index" for c in assoc]
        rbr += ["ScatteredElectronsTruth_objIdx/ScatteredElectronsTruth_objIdx.index",
                "ReconstructedElectrons_objIdx/ReconstructedElectrons_objIdx.index"]
        r = rt.arrays(rbr, library="ak")
        r_index = {int(n[0]): k for k, n in enumerate(r["EventHeader/EventHeader.eventNumber"])}

    rows = []
    for k in range(len(s)):
        num = int(s["EventHeader/EventHeader.eventNumber"][k][0])
        st = ak.to_numpy(s[mc[0]][k])
        prim = np.flatnonzero(st == 1)
        anc = ancestors(st, ak.to_numpy(s[mc[2]][k]), ak.to_numpy(s[mc[3]][k]), ak.to_numpy(s[mc[4]][k]))
        row = {"num": num, "n_primaries": len(prim)}
        idx = {"lead": prim[1], "rec": prim[2]}  # converter order: e', lead, recoil, [remnant]
        for role, j in idx.items():
            for c in TRACKERS:
                z = ak.to_numpy(s[f"{c}/{c}.position.z"][k])[ak.to_numpy(s[f"_{c}_particle/_{c}_particle.index"][k]) == j]
                row[f"{role}_{c}"] = stations(z)
            e = 0.0
            for c in ZDC_HITS:
                ce = ak.to_numpy(s[f"{c}Contributions/{c}Contributions.energy"][k])
                cp = ak.to_numpy(s[f"_{c}Contributions_particle/_{c}Contributions_particle.index"][k])
                e += float(ce[anc[cp] == j].sum()) if len(cp) else 0.0
            row[f"{role}_zdc_edep"] = e
        if r is not None and num in r_index:
            q = r_index[num]
            rst = ak.to_numpy(r[mc[0]][q])
            ranc = ancestors(rst, ak.to_numpy(r[mc[2]][q]), ak.to_numpy(r[mc[3]][q]), ak.to_numpy(r[mc[4]][q]))
            rprim = np.flatnonzero(rst == 1)
            ridx = {"lead": rprim[1], "rec": rprim[2]} if len(rprim) == len(prim) else None
            row["reco"] = ridx is not None
            row["e_truth_matched"] = len(r["ScatteredElectronsTruth_objIdx/ScatteredElectronsTruth_objIdx.index"][q]) > 0
            row["e_reco"] = len(r["ReconstructedElectrons_objIdx/ReconstructedElectrons_objIdx.index"][q]) > 0
            for role in ROLES:
                j = ridx[role] if ridx else -2
                sim_of = lambda c: ak.to_numpy(r[f"_{c}_sim/_{c}_sim.index"][q])  # noqa: E731
                row[f"{role}_b0track"] = bool(np.any(sim_of("B0TrackerCKFTrackAssociations") == j))
                row[f"{role}_charged"] = bool(np.any(sim_of("ReconstructedChargedParticleAssociations") == j))
                for group, cols in (("zdc_cluster", ZDC_CLUSTERS), ("central_cluster", CENTRAL_CLUSTERS)):
                    hit = False
                    for c in cols:
                        if f"_{c}_sim/_{c}_sim.index" in r.fields:
                            m = sim_of(c)
                            hit |= bool(np.any(ranc[m[m >= 0]] == j)) if len(m) else False
                    row[f"{role}_{group}"] = hit
        else:
            row["reco"] = False
        rows.append(row)
    return rows


def boost_angles(p_irf, a, p_nucleon, m_ion):
    """Lab polar angle w.r.t. the ion beam [mrad] and lab momentum of rest-frame 3-momenta."""
    bg = a * p_nucleon / m_ion
    g = np.hypot(1.0, bg)
    px, py, pz = -p_irf[:, 0], p_irf[:, 1], -p_irf[:, 2]
    e = np.sqrt(0.93892**2 + (p_irf**2).sum(1))
    pz_lab = g * pz + bg * e
    return np.arctan2(np.hypot(px, py), pz_lab) * 1e3, np.sqrt(px**2 + py**2 + pz_lab**2)


def extract(run, jobs):
    run = Path(run)
    conv = json.loads((run / "gen.hepmc3.json").read_text())
    beams = conv["beams"]
    sims = sorted(run.glob("sim_*.edm4hep.root"))
    if not sims:
        sys.exit(f"{run}: no sim files")
    with ProcessPoolExecutor(jobs) as ex:
        rows = [r for part in ex.map(extract_file, sims) for r in part]
    num = np.array([r["num"] for r in rows])
    gen = uproot.open(run / "gen.root")["genT"].arrays(["pe", "pLead", "pRec", "lead_type", "rec_type", "weight"], library="np")
    e_irf = beams["E_IRF_run_GeV"]
    pe = gen["pe"][num]
    ee = np.linalg.norm(pe, axis=1)
    q2 = 2 * e_irf * (ee - pe[:, 2])
    xb = q2 / (2 * 0.938272 * (e_irf - ee))
    out = {"num": num, "w_nb": gen["weight"][num] / conv["attempts"], "q2": q2, "xb": xb,
           "lead_type": gen["lead_type"][num], "rec_type": gen["rec_type"][num]}
    z, a, m = {"C12": (6, 12, 11.174863), "d": (1, 2, 1.875613), "He4": (2, 4, 3.727379)}[beams["ion"]]
    for role, b in (("lead", "pLead"), ("rec", "pRec")):
        p = gen[b][num]
        out[f"{role}_p_irf"] = np.linalg.norm(p, axis=1)
        out[f"{role}_theta_mrad"], out[f"{role}_p_lab"] = boost_angles(p, a, beams["ion_momentum_per_nucleon_GeV"], m)
    keys = [k for k in rows[0] if k not in ("num",)]
    for k in keys:
        out[k] = np.array([r.get(k, False) for r in rows])
    meta = {"run": str(run), "beams": beams, "attempts": conv["attempts"], "cuts": conv["cuts"],
            "events": len(num), "q2_range": json.loads((run / "run.json").read_text()).get("q2_range_GeV2")
            if (run / "run.json").exists() else None}
    np.savez_compressed(run / "analysis_events.npz", meta=json.dumps(meta), **out)
    bad = int(np.sum(out["n_primaries"] < 3))
    print(f"{run.name}: {len(num)} events -> analysis_events.npz (reco matched {out['reco'].mean():.3f}"
          f"{f', {bad} events with < 3 primaries' if bad else ''})")


# ----------------------------------------------------------------------------- selections

def detected(d, role, level):
    """Far-forward detection of a nucleon (role: lead|rec) at 'hit' or 'reco' level."""
    is_p = d[f"{role}_type"] == 2212
    if level == "hit":
        p_ok = (d[f"{role}_B0TrackerHits"] >= 3) | (d[f"{role}_ForwardOffMTrackerHits"] >= 2) | (d[f"{role}_ForwardRomanPotHits"] >= 2)
        n_ok = d[f"{role}_zdc_edep"] > 0
    else:
        p_ok = d[f"{role}_b0track"] | (d[f"{role}_ForwardOffMTrackerHits"] >= 2) | (d[f"{role}_ForwardRomanPotHits"] >= 2)
        n_ok = d[f"{role}_zdc_cluster"]
    return np.where(is_p, p_ok, n_ok)


def detected_anywhere(d, role):
    """Far-forward (hit level) or central/B0 reconstruction (charged particle or cluster)."""
    return detected(d, role, "hit") | d[f"{role}_charged"] | d[f"{role}_central_cluster"]


def load(dirs):
    parts = [np.load(Path(r) / "analysis_events.npz") for r in dirs]
    keys = [k for k in parts[0].files if k != "meta"]
    d = {k: np.concatenate([p[k] for p in parts]) for k in keys}
    metas = [json.loads(str(p["meta"])) for p in parts]
    return d, metas


# ----------------------------------------------------------------------------- plots

BLUE, ORANGE, AQUA, INK, MUTED = "#2a78d6", "#eb6834", "#1baf7a", "#0b0b0b", "#52514e"


def style(ax):
    ax.grid(alpha=0.25, lw=0.5)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def wstats(x, w, ok, bins):
    gen, _ = np.histogram(x, bins, weights=w)
    acc, _ = np.histogram(x[ok], bins, weights=w[ok])
    s1, _ = np.histogram(x[ok], bins, weights=w[ok] ** 2)
    s0, _ = np.histogram(x[~ok], bins, weights=w[~ok] ** 2)
    with np.errstate(invalid="ignore", divide="ignore"):
        eff = acc / gen
        err = np.sqrt(s1 * (1 - eff) ** 2 + s0 * eff**2) / gen
    return gen, acc, eff, err


def plots(settings, outdir, lumi_nb_inv, transparency):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    data = {lab: load(dirs) for lab, dirs in settings.items()}
    scale = lumi_nb_inv * transparency
    summary = {}

    # The figures follow the style of PRC 105, 034001 (Figs. 3-6). The paper's evaporation nucleons
    # (from BeAGLE) are not simulated here, so Figs. 3 and 4 show recoil and leading nucleons only.
    green, red, orange, magenta, brown = "#45c045", "#f0281e", "#ff9800", "#ee5fbf", "#9a5a1c"
    low_energy = {lab: metas[0]["beams"]["ion_momentum_per_nucleon_GeV"] < 75 for lab, (_, metas) in data.items()}
    for lab in data:
        summary.setdefault(lab, {})

    def place_label(ax, text, color, pts, anchor, taken, bw=0.18, bh=0.13, **kw):
        """Put text at the free spot nearest anchor (axes fractions): no data point inside its box."""
        best = None
        for xf in np.linspace(0.1, 0.9, 41):
            for yf in np.linspace(0.1, 0.78, 35):
                if np.any((np.abs(pts[:, 0] - xf) < bw / 2) & (np.abs(pts[:, 1] - yf) < bh / 2)):
                    continue
                if any(abs(xf - a) < (bw + aw) / 2 and abs(yf - b) < (bh + ah) / 2 for a, b, aw, ah in taken):
                    continue
                dd = (xf - anchor[0]) ** 2 + (yf - anchor[1]) ** 2
                if best is None or dd < best[0]:
                    best = (dd, xf, yf)
        xf, yf = best[1:] if best else anchor
        taken.append((xf, yf, bw, bh))
        ax.text(xf, yf, text, transform=ax.transAxes, color=color, ha="center", va="center", **kw)

    def frame(ax, minor=True):
        ax.tick_params(direction="in", which="both", top=True, right=True, labelsize=12)
        if minor:
            ax.minorticks_on()

    # Fig. 3: angle w.r.t. the ion beam of recoil (green) and leading (red) nucleons, 10 fb^-1, no p_IRF cut
    fig, axes = plt.subplots(len(data), 1, figsize=(6.6, 3.9 * len(data)), squeeze=False)
    bins = np.arange(0, 51, 1.0)
    cen, hw3 = 0.5 * (bins[1:] + bins[:-1]), 0.5 * np.diff(bins)
    for ax, (lab, (d, _)) in zip(axes[:, 0], data.items()):
        w = d["w_nb"] * lumi_nb_inv
        curves, errs = {}, {}
        for role, c, mk in (("rec", green, "o"), ("lead", red, "^")):
            x = d[f"{role}_theta_mrad"]
            h = np.histogram(x, bins, weights=w)[0] / 1e3
            e = np.sqrt(np.histogram(x, bins, weights=w**2)[0]) / 1e3
            g = h > 0
            ax.errorbar(cen[g], h[g], xerr=hw3[g], yerr=e[g], fmt=mk, ms=4, color=c, elinewidth=1, capsize=0)
            curves[role], errs[role] = h, e
        ymax = 1.15 * max(float(np.max(v)) for v in curves.values())
        ax.set_xlim(0, 50)
        ax.set_ylim(0, ymax)
        frame(ax)
        for x0 in (20, 30):
            ax.axvline(x0, color=orange if x0 == 20 else magenta, lw=1, ls=":")
        y_arrow = 0.965 * ymax
        ax.annotate("", xy=(0.3, y_arrow), xytext=(19.7, y_arrow), arrowprops=dict(arrowstyle="<->", color=orange, lw=1.8))
        ax.annotate("", xy=(49.5, y_arrow), xytext=(30.3, y_arrow), arrowprops=dict(arrowstyle="->", color=magenta, lw=1.8))
        ax.text(4.2, 1.01, "far-forward region", transform=ax.get_xaxis_transform(), ha="left", va="bottom",
                color=orange, fontsize=12, fontweight="bold")
        ax.text(40, 1.01, "central hadron endcap", transform=ax.get_xaxis_transform(), ha="center", va="bottom",
                color=magenta, fontsize=12, fontweight="bold")
        ax.text(0.0, 1.01, r"$\times10^3$", transform=ax.transAxes, ha="left", va="bottom", fontsize=12)
        zdc = 4.0
        ax.plot([zdc, zdc], [-0.07, 0.04], transform=ax.get_xaxis_transform(), color=brown, lw=2, clip_on=False)
        ax.text(zdc, -0.2, "ZDC limit", transform=ax.get_xaxis_transform(), ha="center", va="bottom",
                color=brown, fontsize=12, fontweight="bold")
        ax.text(0.92, 0.84, lab, transform=ax.transAxes, ha="right", fontsize=12, fontweight="bold")
        # direct labels in free space next to each curve (data points incl. error bars kept clear)
        pts = []
        for role in ("rec", "lead"):
            h, e = curves[role], errs[role]
            for y in (h, h + e, np.maximum(h - e, 0)):
                g = h > 0
                pts.append(np.column_stack([cen[g] / 50, y[g] / ymax]))
        pts = np.vstack(pts + [np.array([[0.92, 0.86], [0.85, 0.86]])])  # keep the setting tag clear
        taken = []
        for role, c, text, anchor in (("rec", green, "Recoil\nNucleon", lambda k: (cen[k] / 50 + 0.14, 0.62 * curves["rec"][k] / ymax)),
                                       ("lead", red, "Leading\nNucleon", lambda k: (cen[k] / 50, 0.25))):
            k = int(np.argmax(curves[role]))
            place_label(ax, text, c, pts, anchor(k), taken, fontsize=12, fontweight="bold")
        ax.set_ylabel("Events per bin", fontsize=13)
        ax.set_xlabel(r"$\theta$ [mrad]", fontsize=13, loc="right")
    fig.tight_layout(h_pad=2.5)
    fig.savefig(outdir / "fig3_angles.pdf")
    plt.close(fig)

    # Fig. 4: recoil momentum in the ion rest frame vs opening angle to the ion beam (jet map, empty bins white)
    fig, axes = plt.subplots(len(data), 1, figsize=(4.6, 2.9 * len(data)), squeeze=False)
    for ax, (lab, (d, _)) in zip(axes[:, 0], data.items()):
        xmax = 15 if low_energy[lab] else 8
        H, xe, ye = np.histogram2d(d["rec_theta_mrad"], d["rec_p_irf"], bins=[np.linspace(0, xmax, 61), np.linspace(0, 0.8, 41)],
                                   weights=d["w_nb"] * lumi_nb_inv)
        ax.pcolormesh(xe, ye, np.ma.masked_less_equal(H, 0).T, cmap="jet", vmin=0, shading="flat", rasterized=True)
        ax.axhline(P_IRF_MIN, color="red", lw=2, ls=(0, (5, 3)))
        ax.text(0.05, 0.9, lab, transform=ax.transAxes, fontsize=15, va="top", color="#333333",
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=2))
        ax.text(0.95, 0.08, "Recoil\nNucleons", transform=ax.transAxes, fontsize=13, ha="right", va="bottom", color="#333333")
        ax.set_xticks(np.arange(0, xmax + 1, 5 if xmax >= 10 else 2))
        ax.set_xlim(0, xmax)
        ax.set_ylim(0, 0.8)
        ax.set_xlabel("Opening Angle [mrad]", fontsize=12, color="#333333")
        ax.set_ylabel("Momentum (IRF) [GeV]", fontsize=12, color="#333333")
        ax.tick_params(labelsize=11, colors="#333333")
    fig.tight_layout()
    fig.savefig(outdir / "fig4_pirf_vs_theta.pdf")
    plt.close(fig)

    # Fig. 5: Q2 of events with a recoil detected far forward, 10 fb^-1 x transparency, one point per GeV^2
    fig, ax = plt.subplots(figsize=(6.6, 3.8))
    bins = np.arange(3, 25.5, 1.0)
    cen5 = 0.5 * (bins[1:] + bins[:-1])
    curves5 = {}
    for lab, (d, _) in data.items():
        m = (d["rec_p_irf"] > P_IRF_MIN) & detected(d, "rec", "hit")
        h, _ = np.histogram(d["q2"][m], bins, weights=d["w_nb"][m] * scale)
        summary[lab]["q2_reach_100_per_GeV2"] = float(bins[1:][h >= 100].max()) if (h >= 100).any() else None
        c, mk = (red, "o") if low_energy[lab] else ("#1029ff", "s")
        g = h > 0
        ax.plot(cen5[g], h[g], ls="none", marker=mk, ms=5, color=c)
        curves5[lab] = (h, c)
    ax.set_yscale("log")
    ax.set_xlim(1, 25)
    lo = min(float(np.min(h[h > 0])) for h, _ in curves5.values())
    hi = max(float(np.max(h)) for h, _ in curves5.values())
    ax.set_ylim(lo / 2, hi * 2)
    # direct labels: the higher curve labeled above itself, the lower below
    order = sorted(curves5, key=lambda l: -float(np.sum(curves5[l][0])))
    for i, lab in enumerate(order):
        h, c = curves5[lab]
        k = 12 if i == 0 else 3
        ax.text(cen5[k], h[k] * (9 if i == 0 else 0.2), lab, color=c, fontsize=15, ha="center", va="center")
    frame(ax)
    ax.set_ylabel("Events per bin", fontsize=13, loc="top")
    ax.set_xlabel(r"Q$^2$ [GeV$^2$/c$^4$]", fontsize=13, loc="right")
    fig.tight_layout()
    fig.savefig(outdir / "fig5_q2.pdf")
    plt.close(fig)

    # Fig. 6 in the style of PRC 105, 034001: rows are beam settings, columns recoil p and n, panels
    # share edges. Generated recoils over the full p_IRF range (black crosses), accepted at hit level
    # above the p_IRF cut (filled), accepted at reco level (blue crosses, in the slot of the paper's
    # smeared points), all for 10 fb^-1 x transparency. fig6b: acceptance vs p_IRF.
    pink, pblue = "#f4c4c4", "#1f3fe0"
    lbins = np.linspace(0.0, 1.3, 40)
    lc, hw = 0.5 * (lbins[1:] + lbins[:-1]), 0.5 * np.diff(lbins)
    hists = {}
    for lab, (d, _) in data.items():
        for code, name in ((2212, "proton"), (2112, "neutron")):
            m = d["rec_type"] == code
            x, w = d["rec_p_irf"][m], d["w_nb"][m] * scale
            cut = x > P_IRF_MIN
            sel = {"gen": np.ones(len(x), bool), "hit": cut & detected(d, "rec", "hit")[m],
                   "reco": cut & detected(d, "rec", "reco")[m], "cut": cut}
            hists[lab, name] = {k: (np.histogram(x[v], lbins, weights=w[v])[0],
                                    np.sqrt(np.histogram(x[v], lbins, weights=w[v] ** 2)[0])) for k, v in sel.items()}
            hists[lab, name]["w2"] = {k: np.histogram(x[cut & ~v], lbins, weights=w[cut & ~v] ** 2)[0] for k, v in sel.items()
                                      if k in ("hit", "reco")}
            summary[lab][f"recoil {name[0]}"] = {k: float(w[sel[k]].sum() / w[cut].sum()) for k in ("hit", "reco")}

    def paper_grid(ylabel):
        n = len(data)
        fig, axes = plt.subplots(n, 2, figsize=(8.4, 2.8 * n + 0.9), sharey="row", squeeze=False,
                                 gridspec_kw={"wspace": 0, "hspace": 0})
        for row in range(n):
            for col in range(2):
                ax = axes[row, col]
                ax.tick_params(direction="in", which="both", top=True, right=True, labelsize=11)
                ax.set_xlim(0, 1.3)
                ax.set_xticks([0.0, 0.5, 1.0] if col == 0 else [0.5, 1.0])
                if row < n - 1:
                    ax.tick_params(labelbottom=False)
            axes[row, 0].set_ylabel(ylabel, fontsize=12)
        fig.subplots_adjust(left=0.13, right=0.98, top=1 - 0.45 / fig.get_figheight(), bottom=0.55 / fig.get_figheight())
        fig.supxlabel(r"Spectator Momentum (Ion Rest Frame) [GeV/$c$]", fontsize=12, y=0.1 / fig.get_figheight())
        return fig, axes

    def labels(ax, letter, text):
        ax.text(0.05, 0.93, f"({letter})", transform=ax.transAxes, va="top", fontsize=13)
        ax.text(0.95, 0.93, text, transform=ax.transAxes, va="top", ha="right", fontsize=12)

    def spectra(log, letters, fname):
        fig, axes = paper_grid("Events per bin")
        top = max(float(np.max(hists[k]["gen"][0] + hists[k]["gen"][1])) for k in hists)
        k10 = 0 if log else int(np.floor(np.log10(top)))
        f = 10.0 ** k10
        it = iter(letters)
        for row, lab in enumerate(data):
            ymax = 0.0
            for col, name in enumerate(("proton", "neutron")):
                ax, H = axes[row, col], hists[lab, name]
                ax.stairs(H["hit"][0] / f, lbins, fill=True, color=pink, lw=0, label="Accepted", zorder=1)
                for key, c, lbl in (("gen", INK, "Generated"), ("reco", pblue, "Reconstructed")):
                    v, e = H[key]
                    g = v > 0
                    ax.errorbar(lc[g], v[g] / f, xerr=hw[g], yerr=e[g] / f, fmt="none", ecolor=c, elinewidth=1.4,
                                label=lbl, zorder=2)
                    if g.any():
                        ymax = max(ymax, float(np.max(v[g] + e[g])) / f)
                labels(ax, next(it), f"{name}, {lab}")
            if log:
                axes[row, 0].set_yscale("log")
                axes[row, 0].set_ylim(1, ymax * 8)
            else:
                axes[row, 0].set_ylim(0, ymax * 1.3)
        if not log and k10:
            axes[0, 0].text(0.0, 1.02, rf"$\times10^{{{k10}}}$", transform=axes[0, 0].transAxes, fontsize=12, va="bottom")
        h, l_ = axes[0, 1].get_legend_handles_labels()
        order = [l_.index(x) for x in ("Generated", "Accepted", "Reconstructed")]
        h, l_ = [h[i] for i in order], [l_[i] for i in order]
        if log:  # the paper's log panels carry no legend; put it in one row above the grid
            fig.legend(h, l_, frameon=False, fontsize=11, ncol=3, loc="upper center", bbox_to_anchor=(0.555, 1.0))
        else:
            axes[0, 1].legend(h, l_, frameon=False, fontsize=11, loc="upper right", bbox_to_anchor=(1.0, 0.85))
        fig.savefig(outdir / fname)
        plt.close(fig)

    spectra(True, "abcdefgh", "fig6_recoil_acceptance.pdf")
    spectra(False, "efghijkl", "fig6_recoil_acceptance_linear.pdf")

    # fig6b: acceptance vs p_IRF (accepted / generated above the p_IRF cut), hit and reco level
    fig, axes = paper_grid("Acceptance")
    it = iter("abcdefgh")
    for row, lab in enumerate(data):
        for col, name in enumerate(("proton", "neutron")):
            ax, H = axes[row, col], hists[lab, name]
            gen = H["cut"][0]
            ok = gen > 0
            for j, (key, c, lbl) in enumerate((("hit", INK, "hit level"), ("reco", pblue, "reco level"))):
                acc = H[key][0]
                eff = np.divide(acc, gen, out=np.zeros_like(acc), where=ok)
                err = np.sqrt(H[key][1] ** 2 * (1 - eff) ** 2 + H["w2"][key] * eff ** 2) / np.where(ok, gen, 1)
                ax.errorbar(lc[ok] + (j - 0.5) * 0.25 * hw[ok], eff[ok], yerr=err[ok], fmt="o", ms=3.5, color=c,
                            elinewidth=1.2, label=lbl)
            ax.axvline(P_IRF_MIN, color=MUTED, lw=0.8, ls=":")
            labels(ax, next(it), f"{name}, {lab}: {100 * summary[lab]['recoil ' + name[0]]['hit']:.0f}%")
        axes[row, 0].set_ylim(0, 1.15)
    axes[0, 1].legend(frameon=False, fontsize=11, loc="upper right", bbox_to_anchor=(1.0, 0.85))
    fig.savefig(outdir / "fig6b_acceptance_vs_pirf.pdf")
    plt.close(fig)

    # Pair-type table: recoil-tagged and triple-coincidence fractions
    lines = []
    for lab, (d, metas) in data.items():
        m0 = d["rec_p_irf"] > P_IRF_MIN
        w = d["w_nb"]
        rec_ok = detected(d, "rec", "hit")
        lead_ok = detected_anywhere(d, "lead")
        e_ok = d["e_truth_matched"] if "e_truth_matched" in d else np.ones(len(w), bool)
        sigma = float(w[m0].sum())
        summary[lab]["sigma_nb_cuts_pIRF"] = sigma
        summary[lab]["electron_truth_matched"] = float(w[m0 & e_ok].sum() / w[m0].sum())
        summary[lab]["runs"] = [m_["run"] for m_ in metas]
        for lt, rt, pair in ((2212, 2212, "pp"), (2212, 2112, "pn"), (2112, 2212, "np"), (2112, 2112, "nn")):
            m = m0 & (d["lead_type"] == lt) & (d["rec_type"] == rt)
            if not m.any():
                continue
            f = lambda ok: float(w[m & ok].sum() / w[m].sum())  # noqa: E731
            summary[lab][f"pair_{pair}"] = {"weight_fraction": float(w[m].sum() / sigma), "recoil": f(rec_ok),
                                            "lead": f(lead_ok), "triple_e_lead_rec": f(rec_ok & lead_ok & e_ok)}
            lines.append(f"{lab:8s} {pair}  share {w[m].sum() / sigma:5.3f}  recoil {f(rec_ok):5.2f}  "
                         f"lead {f(lead_ok):5.2f}  e+lead+recoil {f(rec_ok & lead_ok & e_ok):5.2f}")
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (outdir / "pairs.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    for lab in data:
        print(f"{lab}: recoil p {summary[lab]['recoil p']}, recoil n {summary[lab]['recoil n']}, "
              f"Q2 reach (>=100/GeV2) {summary[lab]['q2_reach_100_per_GeV2']}, "
              f"electron {summary[lab]['electron_truth_matched']:.2f}")
    print(f"wrote {outdir}/fig3_angles.pdf fig4_pirf_vs_theta.pdf fig5_q2.pdf fig6_recoil_acceptance.pdf fig6_recoil_acceptance_linear.pdf fig6b_acceptance_vs_pirf.pdf summary.json pairs.txt")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    e = sub.add_parser("extract", help="reduce run directories to analysis_events.npz")
    e.add_argument("runs", nargs="+")
    e.add_argument("-j", "--jobs", type=int, default=16)
    p = sub.add_parser("plots", help="figures and summary from extracted runs")
    p.add_argument("--setting", action="append", required=True, help="LABEL=RUN_DIR[,RUN_DIR...] (Q2 slices of one setting)")
    p.add_argument("-o", "--outdir", required=True)
    p.add_argument("--lumi-fb", type=float, default=10.0, help="integrated luminosity per nucleus [fb^-1] (default 10)")
    p.add_argument("--transparency", type=float, default=0.5, help="flat nucleon transparency factor (default 0.5)")
    a = ap.parse_args()
    if a.cmd == "extract":
        for r in a.runs:
            extract(r, a.jobs)
    else:
        settings = {}
        for s in a.setting:
            lab, dirs = s.split("=", 1)
            settings[lab] = dirs.split(",")
        plots(settings, a.outdir, a.lumi_fb * 1e6, a.transparency)


if __name__ == "__main__":
    main()
