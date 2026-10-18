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
    from matplotlib.colors import LogNorm

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    data = {lab: load(dirs) for lab, dirs in settings.items()}
    colors = dict(zip(data, (BLUE, ORANGE, AQUA)))
    scale = lumi_nb_inv * transparency
    summary = {}

    # Fig. 3: angles of leading and recoil nucleons w.r.t. the ion beam
    fig, axes = plt.subplots(1, len(data), figsize=(5.2 * len(data), 4.2), squeeze=False)
    for ax, (lab, (d, _)) in zip(axes[0], data.items()):
        bins = np.linspace(0, 50, 101)
        for role, c, name in (("rec", BLUE, "Recoil nucleon"), ("lead", ORANGE, "Leading nucleon")):
            m = d["rec_p_irf"] > P_IRF_MIN
            h, _ = np.histogram(d[f"{role}_theta_mrad"][m], bins, weights=d["w_nb"][m] * lumi_nb_inv)
            ax.stairs(h, bins, color=c, lw=2, label=name)
        ax.axvspan(0, 20, color=MUTED, alpha=0.08, lw=0)
        ax.axvline(20, color=MUTED, lw=0.8, ls=":")
        ax.axvline(30, color=MUTED, lw=0.8, ls=":")
        ax.text(10, 1.01, "far forward", transform=ax.get_xaxis_transform(), ha="center", va="bottom", color=MUTED, fontsize=8)
        ax.text(40, 1.01, "central endcap", transform=ax.get_xaxis_transform(), ha="center", va="bottom", color=MUTED, fontsize=8)
        ax.set_yscale("log")
        ax.set_ylim(bottom=1)
        ax.set_xlabel(r"$\theta$ w.r.t. ion beam [mrad]")
        ax.set_ylabel(r"events / 0.5 mrad (10 fb$^{-1}$)")
        ax.set_title(f"$^{{12}}$C$(e,e'NN)$, {lab}", fontsize=11, pad=16)
        ax.legend(frameon=False, fontsize=9, loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=2)
        style(ax)
    fig.tight_layout()
    fig.savefig(outdir / "fig3_angles.pdf")
    plt.close(fig)

    # Fig. 4: p_IRF vs angle of the recoil nucleon
    fig, axes = plt.subplots(1, len(data), figsize=(5.2 * len(data), 3.8), squeeze=False)
    for ax, (lab, (d, _)) in zip(axes[0], data.items()):
        _, _, _, im = ax.hist2d(d["rec_theta_mrad"], d["rec_p_irf"], bins=[np.linspace(0, 30, 61), np.linspace(0, 1.2, 49)],
                                weights=d["w_nb"] * lumi_nb_inv, cmap="Blues", norm=LogNorm(vmin=1))
        fig.colorbar(im, ax=ax, label=r"events / bin (10 fb$^{-1}$)")
        ax.axhline(P_IRF_MIN, color=ORANGE, lw=1.5, ls="--")
        ax.set_xlabel(r"recoil $\theta$ w.r.t. ion beam [mrad]")
        ax.set_ylabel(r"recoil $p_\mathrm{IRF}$ [GeV/$c$]")
        ax.set_title(f"recoil nucleons, {lab}", fontsize=11)
        style(ax)
    fig.tight_layout()
    fig.savefig(outdir / "fig4_pirf_vs_theta.pdf")
    plt.close(fig)

    # Fig. 5: Q2 of events with a detected recoil (hit level), 10 fb^-1 x transparency
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    bins = np.arange(3, 25.5, 1.0)
    for lab, (d, _) in data.items():
        m = (d["rec_p_irf"] > P_IRF_MIN) & detected(d, "rec", "hit")
        h, _ = np.histogram(d["q2"][m], bins, weights=d["w_nb"][m] * scale)
        ax.stairs(h, bins, color=colors[lab], lw=2, label=lab)
        summary.setdefault(lab, {})["q2_reach_100_per_GeV2"] = float(bins[1:][h >= 100].max()) if (h >= 100).any() else None
    ax.axhline(100, color=MUTED, lw=0.8, ls=":")
    ax.set_yscale("log")
    ax.set_ylim(bottom=0.1)
    ax.set_xlabel(r"$Q^2$ [GeV$^2$]")
    ax.set_ylabel(rf"events / GeV$^2$ ({lumi_nb_inv / 1e6:g} fb$^{{-1}}$ $\times$ {transparency:g})")
    ax.set_title(r"$^{12}$C$(e,e'N_\mathrm{rec})$, recoil detected far forward, $x_B>1.2$, $p_\mathrm{IRF}>0.3$", fontsize=10)
    ax.legend(frameon=False)
    style(ax)
    fig.tight_layout()
    fig.savefig(outdir / "fig5_q2.pdf")
    plt.close(fig)

    # Fig. 6: recoil momentum (rest frame), generated / accepted (hit) / accepted (reco), and acceptance
    pbins = np.linspace(0.3, 1.2, 19)
    centers = 0.5 * (pbins[1:] + pbins[:-1])
    fig, axes = plt.subplots(2, 2 * len(data), figsize=(4.2 * 2 * len(data), 6.0), sharex=True,
                             gridspec_kw={"height_ratios": [3, 2]}, squeeze=False)
    col = 0
    for lab, (d, _) in data.items():
        for code, name in ((2212, "recoil p"), (2112, "recoil n")):
            m = (d["rec_type"] == code) & (d["rec_p_irf"] > P_IRF_MIN)
            x, w = d["rec_p_irf"][m], d["w_nb"][m] * scale
            top, bot = axes[0, col], axes[1, col]
            res = {}
            for level, c, ls in (("hit", BLUE, "-"), ("reco", ORANGE, "--")):
                ok = detected(d, "rec", level)[m]
                gen, acc, eff, err = wstats(x, w, ok, pbins)
                if level == "hit":
                    top.stairs(gen, pbins, color=INK, lw=1.5, label="generated")
                top.stairs(acc, pbins, color=c, lw=2, ls=ls, label=f"accepted ({level} level)")
                bot.errorbar(centers, eff, yerr=err, fmt="o", ms=4, color=c, label=level)
                res[level] = float(w[ok].sum() / w.sum())
            summary[lab][name] = res
            top.set_yscale("log")
            top.set_ylim(bottom=1)
            top.set_title(f"{name}, {lab}: {100 * res['hit']:.0f}% (hit), {100 * res['reco']:.0f}% (reco)", fontsize=10)
            top.set_ylabel(r"events / bin (10 fb$^{-1}$ $\times$ 0.5)" if col == 0 else "")
            bot.set_ylim(0, 1.05)
            bot.set_xlabel(r"$p_\mathrm{IRF}$ [GeV/$c$]")
            bot.set_ylabel("acceptance" if col == 0 else "")
            for a_ in (top, bot):
                style(a_)
            col += 1
    axes[0, 0].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "fig6_recoil_acceptance.pdf")
    plt.close(fig)

    # Fig. 6, linear, laid out like PRC 105, 034001 Fig. 6(e)-(h): generated over the full p_IRF
    # range (points), accepted at hit level (filled), accepted at reco level (points), 10 fb^-1 x 0.5
    lbins = np.linspace(0.0, 1.3, 40)
    lc = 0.5 * (lbins[1:] + lbins[:-1])
    hw = 0.5 * np.diff(lbins)
    fig, axes = plt.subplots(len(data), 2, figsize=(9.6, 3.4 * len(data)), sharex=True, sharey="row", squeeze=False)
    letters = iter("efghijklmnop")
    for row, (lab, (d, _)) in enumerate(data.items()):
        ymax = 0.0
        for col, (code, name) in enumerate(((2212, "proton"), (2112, "neutron"))):
            ax = axes[row, col]
            m = d["rec_type"] == code
            x, w = d["rec_p_irf"][m], d["w_nb"][m] * scale
            gen, _ = np.histogram(x, lbins, weights=w)
            gerr = np.sqrt(np.histogram(x, lbins, weights=w**2)[0])
            cut = x > P_IRF_MIN
            hit = cut & detected(d, "rec", "hit")[m]
            rec = cut & detected(d, "rec", "reco")[m]
            acc, _ = np.histogram(x[hit], lbins, weights=w[hit])
            rco, _ = np.histogram(x[rec], lbins, weights=w[rec])
            rerr = np.sqrt(np.histogram(x[rec], lbins, weights=w[rec] ** 2)[0])
            ax.stairs(acc, lbins, fill=True, color=BLUE, alpha=0.35, lw=0, label="Accepted (hit level)")
            g = gen > 0
            ax.errorbar(lc[g], gen[g], xerr=hw[g], yerr=gerr[g], fmt="none", ecolor=INK, elinewidth=1.2, label="Generated")
            r = rco > 0
            ax.errorbar(lc[r], rco[r], xerr=hw[r], yerr=rerr[r], fmt="none", ecolor=ORANGE, elinewidth=1.2,
                        label="Accepted (reco level)")
            ymax = max(ymax, float(np.max(gen + gerr)), float(np.max(acc)))
            ax.text(0.04, 0.92, f"({next(letters)})", transform=ax.transAxes, va="top", fontsize=12)
            ax.text(0.96, 0.92, f"{name}, {lab}", transform=ax.transAxes, va="top", ha="right", fontsize=11)
            ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0), useMathText=True)
            if col == 0:
                ax.set_ylabel(r"Events per bin (10 fb$^{-1}$ $\times$ 0.5)")
            style(ax)
        axes[row, 0].set_ylim(0, 1.12 * ymax)  # rows share y; set after both panels are drawn
    for ax in axes[-1]:
        ax.set_xlim(0, 1.3)
        ax.set_xlabel(r"Spectator momentum (ion rest frame) [GeV/$c$]")
    axes[0, 1].legend(frameon=False, fontsize=9, loc="upper right", bbox_to_anchor=(1.0, 0.84))
    fig.tight_layout()
    fig.savefig(outdir / "fig6_recoil_acceptance_linear.pdf")
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
    print(f"wrote {outdir}/fig3_angles.pdf fig4_pirf_vs_theta.pdf fig5_q2.pdf fig6_recoil_acceptance.pdf fig6_recoil_acceptance_linear.pdf summary.json pairs.txt")


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
