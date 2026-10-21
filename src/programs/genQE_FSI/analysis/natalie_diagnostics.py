#!/usr/bin/env python3
"""Diagnostics requested by N. Wright (2026-09) for the (e,e'pp) transparency > 1
and the mean-field c.m. puzzle:

  D1  multiplicity of final-state nucleons with p > 350 MeV/c
  D2  theta(p_miss, recoil) and theta(p_lead, recoil) for
        (1) actual GCF pairs   - the true SRC partner, untouched by the cascade
        (2) all MF candidates  - every MF recoil candidate (several fills per event)
        (3) selected MF pairs  - the candidate our (e,e'pp) selection keeps
  D3  which FSI process produces the SRC-like MF "pairs" (fate of the recoil's ancestor)

Truth level (no CLAS acceptance, no smearing): this is a sanity check of the GENIE
hA cascade, not a detector comparison. Lead = highest-momentum final-state proton.
SRC cuts: xB>1.2, theta_pq<25 deg, 0.62<|p_lead|/|q|<0.96, 0.4<p_miss<1.0, m_miss<1.1.

Usage: natalie_diagnostics.py [--max-events N] [--outdir DIR]
"""
import argparse, os, sys
import numpy as np
import awkward as ak
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "/home/ilayws/Ilay-Generators/src/programs/genQE_FSI"
SRC_FILE = os.path.join(BASE, "events/natalie/fsi_hA_SRC_rho1.root")
MF_FILE = os.path.join(BASE, "events/natalie/fsi_hA_MF_localFG.root")
M_AVG, M_P = 0.9389, 0.93827
PREC_MIN = 0.35
# Intranuke hA2018 fate codes (INukeHadroFates2018.h; Elas commented out)
FATES = {-1: "no interaction", 0: "undefined", 1: "NoInteraction", 2: "ChargeExchange",
         3: "Inelastic", 4: "Absorption", 5: "KnockOut", 6: "Compound", 7: "PiProd"}
BR = ["weight", "q", "lead", "nu", "xB", "nPart", "part_pdg", "part_px", "part_py",
      "part_pz", "part_E", "part_generation", "part_parentRole", "part_rescatterCode"]


def vec(a):
    return ak.zip({"x": a.part_px, "y": a.part_py, "z": a.part_pz})


def analyse(path, label, max_events, nbins_ang=45):
    """Stream one sample, return accumulated histograms and counters."""
    ang_bins = np.linspace(0, 180, nbins_ang + 1)
    out = dict(
        label=label,
        mult_p=np.zeros(12), mult_p_src=np.zeros(12),          # protons > 350 MeV/c
        mult_N=np.zeros(12), mult_N_src=np.zeros(12),          # nucleons > 350 MeV/c
        ang_miss=dict(), ang_lead=dict(),                       # per category
        fates=dict(), gens=dict(), roles=dict(),
        n_ev=0.0, n_src=0.0, n_epp=0.0, w_src=0.0, w_epp=0.0)
    for cat in ("gcf_true", "all_cand", "selected"):
        out["ang_miss"][cat] = np.zeros(nbins_ang)
        out["ang_lead"][cat] = np.zeros(nbins_ang)

    seen = 0
    for a in uproot.iterate(f"{path}:events", BR, step_size=200_000, library="ak"):
        if max_events and seen >= max_events:
            break
        if max_events and seen + len(a) > max_events:
            a = a[: max_events - seen]
        seen += len(a)
        w = ak.to_numpy(a.weight)
        out["n_ev"] += len(a)

        p3 = vec(a)
        pmag = np.sqrt(p3.x**2 + p3.y**2 + p3.z**2)
        is_p = a.part_pdg == 2212
        is_N = (a.part_pdg == 2212) | (a.part_pdg == 2112)
        hard = pmag > PREC_MIN

        # --- D1 multiplicities ---
        np_hard = ak.to_numpy(ak.sum(is_p & hard, axis=1))
        nN_hard = ak.to_numpy(ak.sum(is_N & hard, axis=1))
        for arr, key in ((np_hard, "mult_p"), (nN_hard, "mult_N")):
            h, _ = np.histogram(np.clip(arr, 0, 11), bins=np.arange(13), weights=w)
            out[key] += h

        # --- lead = highest-momentum final-state proton ---
        pmag_p = ak.where(is_p, pmag, -1.0)
        has_p = ak.to_numpy(ak.sum(is_p, axis=1)) > 0
        ilead = ak.argmax(pmag_p, axis=1, keepdims=True)
        lead = ak.firsts(p3[ilead]);  lead_p = ak.to_numpy(ak.firsts(pmag[ilead]))
        lx, ly, lz = (ak.to_numpy(ak.firsts(p3[ilead][f])) for f in ("x", "y", "z"))
        q = ak.to_numpy(a.q)  # (px,py,pz,E)
        qx, qy, qz = q[:, 0], q[:, 1], q[:, 2]
        qmag = np.sqrt(qx**2 + qy**2 + qz**2)
        nu = ak.to_numpy(a.nu); xB = ak.to_numpy(a.xB)

        mx, my, mz = lx - qx, ly - qy, lz - qz          # p_miss vector
        pmiss = np.sqrt(mx**2 + my**2 + mz**2)
        cos_pq = (lx * qx + ly * qy + lz * qz) / np.maximum(lead_p * qmag, 1e-12)
        th_pq = np.degrees(np.arccos(np.clip(cos_pq, -1, 1)))
        E_lead = np.sqrt(lead_p**2 + M_P**2)
        mm2 = (2 * M_AVG + nu - E_lead) ** 2 - pmiss**2
        mmiss = np.where(mm2 > 0, np.sqrt(np.maximum(mm2, 0.0)), -1.0)

        src = (has_p & (xB > 1.2) & (th_pq < 25) & (lead_p / np.maximum(qmag, 1e-12) > 0.62)
               & (lead_p / np.maximum(qmag, 1e-12) < 0.96) & (pmiss > 0.4) & (pmiss < 1.0)
               & (mmiss < 1.1) & (mmiss > 0))
        out["n_src"] += src.sum(); out["w_src"] += w[src].sum()
        for arr, key in ((np_hard, "mult_p_src"), (nN_hard, "mult_N_src")):
            h, _ = np.histogram(np.clip(arr[src], 0, 11), bins=np.arange(13), weights=w[src])
            out[key] += h

        # --- recoil candidates: protons > 350 MeV/c that are not the lead ---
        idx = ak.local_index(a.part_pdg)
        lead_idx = ak.firsts(ilead)
        cand = is_p & hard & (idx != lead_idx)
        ncand = ak.to_numpy(ak.sum(cand, axis=1))

        def angles(cx, cy, cz, ev):
            """angles of candidate vectors vs p_miss and p_lead, for events `ev`."""
            cm = np.sqrt(cx**2 + cy**2 + cz**2)
            a1 = np.degrees(np.arccos(np.clip((cx * mx[ev] + cy * my[ev] + cz * mz[ev])
                                              / np.maximum(cm * pmiss[ev], 1e-12), -1, 1)))
            a2 = np.degrees(np.arccos(np.clip((cx * lx[ev] + cy * ly[ev] + cz * lz[ev])
                                              / np.maximum(cm * lead_p[ev], 1e-12), -1, 1)))
            return a1, a2

        # (2) all candidates, several fills per event  [SRC-cut events only]
        c3s = p3[cand][src]
        counts = ak.to_numpy(ak.num(c3s))
        ev = np.repeat(np.where(src)[0], counts)
        cx = ak.to_numpy(ak.flatten(c3s.x)); cy = ak.to_numpy(ak.flatten(c3s.y)); cz = ak.to_numpy(ak.flatten(c3s.z))
        if len(cx):
            a1, a2 = angles(cx, cy, cz, ev)
            out["ang_miss"]["all_cand"] += np.histogram(a1, bins=ang_bins, weights=w[ev])[0]
            out["ang_lead"]["all_cand"] += np.histogram(a2, bins=ang_bins, weights=w[ev])[0]

        # (3) selected pair: highest-momentum candidate
        pm_c = ak.where(cand, pmag, -1.0)
        isel = ak.argmax(pm_c, axis=1, keepdims=True)
        sel_ok = src & (ncand > 0)
        sx = ak.to_numpy(ak.firsts(p3[isel].x)); sy = ak.to_numpy(ak.firsts(p3[isel].y)); sz = ak.to_numpy(ak.firsts(p3[isel].z))
        ev = np.where(sel_ok)[0]
        if len(ev):
            a1, a2 = angles(sx[ev], sy[ev], sz[ev], ev)
            out["ang_miss"]["selected"] += np.histogram(a1, bins=ang_bins, weights=w[ev])[0]
            out["ang_lead"]["selected"] += np.histogram(a2, bins=ang_bins, weights=w[ev])[0]
            out["n_epp"] += len(ev); out["w_epp"] += w[ev].sum()
            # --- D3 provenance of the selected recoil ---
            gsel = ak.to_numpy(ak.firsts(a.part_generation[isel]))[ev]
            fsel = ak.to_numpy(ak.firsts(a.part_rescatterCode[isel]))[ev]
            rsel = ak.to_numpy(ak.firsts(a.part_parentRole[isel]))[ev]
            wev = w[ev]
            lfate = ak.to_numpy(ak.firsts(a.part_rescatterCode[ilead]))[ev]
            lgen = ak.to_numpy(ak.firsts(a.part_generation[ilead]))[ev]
            for arr, key in ((gsel, "gens"), (fsel, "fates"), (rsel, "roles")):
                for v in np.unique(arr):
                    out[key][int(v)] = out[key].get(int(v), 0.0) + wev[arr == v].sum()
            for f, g, ww in ((fsel, gsel, wev),):
                for pair in set(zip(f.tolist(), g.tolist())):
                    m = (f == pair[0]) & (g == pair[1])
                    out.setdefault("rec_fate_gen", {})[pair] = out.setdefault("rec_fate_gen", {}).get(pair, 0.0) + ww[m].sum()
            for pair in set(zip(lfate.tolist(), lgen.tolist())):
                m = (lfate == pair[0]) & (lgen == pair[1])
                out.setdefault("lead_fate_gen", {})[pair] = out.setdefault("lead_fate_gen", {}).get(pair, 0.0) + wev[m].sum()

        if "SRC" in label.upper() or "GCF" in label.upper():
            lead_role = ak.to_numpy(ak.firsts(a.part_parentRole[ilead]))
            out["lead_is_partner"] = out.get("lead_is_partner", 0.0) + w[src & (lead_role == 1)].sum()
            out["lead_is_lead"] = out.get("lead_is_lead", 0.0) + w[src & (lead_role == 0)].sum()

        # (1) actual GCF pairs: the true SRC partner, when it is not the lead
        if "SRC" in label.upper() or "GCF" in label.upper():
            # exclude the lead itself: in ~1/6 of events the SRC partner IS the
            # highest-momentum proton, i.e. it was detected as the lead.
            # generation == 1 does NOT mean "untouched": every cascade output
            # is at depth >= 1, and in hA a rescattered recoil or a nucleon it
            # knocked out is also at depth 1 with role 1. So "gcf_true" also
            # holds those; rescatterCode == 1 would isolate the untouched
            # partner.
            true_rec = is_p & (a.part_generation == 1) & (a.part_parentRole == 1) & (idx != lead_idx)
            tr = p3[true_rec]
            ev2 = np.where(src & (ak.to_numpy(ak.sum(true_rec, axis=1)) > 0))[0]
            if len(ev2):
                tx = ak.to_numpy(ak.firsts(tr.x))[ev2]; ty = ak.to_numpy(ak.firsts(tr.y))[ev2]; tz = ak.to_numpy(ak.firsts(tr.z))[ev2]
                keep = np.sqrt(tx**2 + ty**2 + tz**2) > PREC_MIN
                ev2, tx, ty, tz = ev2[keep], tx[keep], ty[keep], tz[keep]
                if len(ev2):
                    a1, a2 = angles(tx, ty, tz, ev2)
                    out["ang_miss"]["gcf_true"] += np.histogram(a1, bins=ang_bins, weights=w[ev2])[0]
                    out["ang_lead"]["gcf_true"] += np.histogram(a2, bins=ang_bins, weights=w[ev2])[0]
    out["ang_bins"] = ang_bins
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-events", type=int, default=5_000_000)
    ap.add_argument("--outdir", default=os.path.join(BASE, "analysis/plots/natalie_diagnostics"))
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    res = {}
    for path, lab in ((SRC_FILE, "GCF (SRC)"), (MF_FILE, "MF")):
        print(f"[{lab}] streaming {os.path.basename(path)} ...", flush=True)
        res[lab] = analyse(path, lab, args.max_events)
        r = res[lab]
        print(f"  events read {r['n_ev']:.0f}   pass SRC cuts (weighted) {r['w_src']:.4g}"
              f"   with a >0.35 recoil {r['w_epp']:.4g}   ratio {r['w_epp']/max(r['w_src'],1e-30):.4f}")
    np.save(os.path.join(args.outdir, "diagnostics.npy"), res, allow_pickle=True)
    make_plots(res, args.outdir)
    return res


def make_plots(res, outdir):
    plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 11,
                         "axes.labelsize": 13, "xtick.direction": "in", "ytick.direction": "in"})
    G, M = res["GCF (SRC)"], res["MF"]
    lines = []

    # --- D1 multiplicity ---
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
    x = np.arange(12)
    for a, keyp, keyN, title in ((ax[0], "mult_p", "mult_N", "all generated events"),
                                 (ax[1], "mult_p_src", "mult_N_src", "after SRC cuts on the lead")):
        for r, col in ((G, "#3366bb"), (M, "#cc3355")):
            h = r[keyp] / max(r[keyp].sum(), 1e-30)
            a.step(x, h, where="mid", color=col, lw=2, label=f"{r['label']}  protons")
            hN = r[keyN] / max(r[keyN].sum(), 1e-30)
            a.step(x, hN, where="mid", color=col, lw=1.2, ls="--", label=f"{r['label']}  p+n")
        a.set_xlabel("final-state nucleons with p > 350 MeV/c"); a.set_xlim(-0.5, 7.5)
        a.set_ylabel("fraction of events"); a.set_title(title, fontsize=11)
    ax[0].legend(frameon=False, fontsize=9)
    fig.tight_layout(); fig.savefig(os.path.join(outdir, "D1_multiplicity.png"), dpi=180); plt.close(fig)

    for r in (G, M):
        for key, lab in (("mult_p_src", "protons"), ("mult_N_src", "nucleons")):
            h = r[key]; tot = max(h.sum(), 1e-30)
            mean = (h * np.arange(12)).sum() / tot
            lines.append(f"{r['label']:10s} {lab:8s} >350 MeV/c after SRC cuts: mean {mean:.2f}, "
                         f"P(>=2) = {100*h[2:].sum()/tot:5.1f} %, P(>=3) = {100*h[3:].sum()/tot:4.1f} %")

    # --- D2 angles ---
    c = 0.5 * (G["ang_bins"][1:] + G["ang_bins"][:-1])
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
    series = [(G, "gcf_true", "#3366bb", "-", "actual GCF pairs (true SRC partner)"),
              (G, "selected", "#88aadd", "-.", "selected GCF pairs"),
              (M, "all_cand", "#cc3355", ":", "all MF candidates (multi-fill)"),
              (M, "selected", "#118844", "--", "selected MF pairs")]
    for a, key, xl in ((ax[0], "ang_miss", r"$\theta(p_{\rm miss},\, p_{\rm recoil})$ [deg]"),
                       (ax[1], "ang_lead", r"$\theta(p_{\rm lead},\, p_{\rm recoil})$ [deg]")):
        for r, cat, col, ls, lab in series:
            h = r[key][cat]
            if h.sum() <= 0:
                continue
            a.step(c, h / h.sum(), where="mid", color=col, ls=ls, lw=2, label=lab)
        a.set_xlabel(xl); a.set_ylabel("normalised"); a.set_xlim(0, 180)
    ax[0].legend(frameon=False, fontsize=9)
    fig.suptitle("Recoil direction, truth level, after SRC cuts on the lead", fontsize=11)
    fig.tight_layout(); fig.savefig(os.path.join(outdir, "D2_recoil_angles.png"), dpi=180); plt.close(fig)

    for r, cat, _, _, lab in series:
        h = r["ang_miss"][cat]
        if h.sum() > 0:
            lines.append(f"{lab:42s} <theta(pmiss,rec)> = {(h*c).sum()/h.sum():6.1f} deg   "
                         f"frac > 120 deg = {100*h[c>120].sum()/h.sum():5.1f} %")

    # --- D3 provenance of the selected MF recoil ---
    fg = M.get("rec_fate_gen", {})
    if fg:
        tot = sum(fg.values())
        items = sorted(fg.items(), key=lambda kv: -kv[1])
        labs = [f"{FATES.get(k[0], k[0])}\ngen {k[1]}" for k, _ in items]
        vals = [100 * v / tot for _, v in items]
        fig, a = plt.subplots(figsize=(8.5, 4.4))
        a.bar(range(len(vals)), vals, color="#cc3355")
        a.set_xticks(range(len(vals))); a.set_xticklabels(labs, fontsize=8, rotation=20)
        a.set_ylabel("% of selected MF (e,e'pp) pairs")
        a.set_title("How the MF recoil was made (ancestor fate, cascade depth)", fontsize=11)
        fig.tight_layout(); fig.savefig(os.path.join(outdir, "D3_mf_recoil_origin.png"), dpi=180); plt.close(fig)
        lines.append("")
        lines.append("MF selected recoil provenance (ancestor fate, generation):")
        for k, v in items:
            lines.append(f"   {FATES.get(k[0], k[0]):16s} gen {k[1]}   {100*v/tot:5.1f} %")
    lf = M.get("lead_fate_gen", {})
    if lf:
        tot = sum(lf.values())
        lines.append("MF lead provenance:")
        for k, v in sorted(lf.items(), key=lambda kv: -kv[1]):
            lines.append(f"   {FATES.get(k[0], k[0]):16s} gen {k[1]}   {100*v/tot:5.1f} %")
    if "lead_is_partner" in G:
        tot = G["lead_is_partner"] + G["lead_is_lead"]
        lines.append(f"GCF: the detected lead is the SRC *partner* (not the struck nucleon) in "
                     f"{100*G['lead_is_partner']/max(tot,1e-30):.1f} % of SRC-cut events")
    for r in (G, M):
        lines.append(f"{r['label']:10s} (e,e'pp)/(e,e'p) at truth level = "
                     f"{r['w_epp']/max(r['w_src'],1e-30):.4f}")
    txt = os.path.join(outdir, "summary.txt")
    open(txt, "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))
    print("\nsaved plots + summary in", outdir)


if __name__ == "__main__":
    main()
