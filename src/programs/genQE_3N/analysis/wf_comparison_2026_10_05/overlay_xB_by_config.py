#!/usr/bin/env python3
"""Overlay the xB distribution of each 3N angular configuration (Lead Rocket,
Recoil Rocket, Star) across wavefunction tables.

Same regions, same cuts and the same weight clipping as
plotting/plot_fig4_xB_by_config.py (whose region code is imported directly).
Top row: yields per generation attempt (so cross-section differences show).
Bottom row: area-normalised shapes.
"""
import os, sys
import numpy as np, uproot
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

G3N = "/home/ilayws/Ilay-Generators/src/programs/genQE_3N"
sys.path.insert(0, os.path.join(G3N, "plotting"))
from regions import in_region

EV = os.path.join(G3N, "events/wf_test")
OUT = os.path.dirname(os.path.abspath(__file__))
REGIONS = [("Lead Rocket", 0), ("Recoil Rocket", 2), ("Star", 6)]
SAMPLES = [("AV18 (u=1)", "u1", "#3366bb", "-"), ("CD-Bonn (u=4)", "u4", "#118844", "-"),
           ("pD AV18 (u=6)", "u6", "#cc3355", "--"), ("pD CD-Bonn (u=9)", "u9", "#ee9944", "--"),
           ("AV8 old", "av8", "#888888", ":")]
THETA_E_MAX, Q2_MIN, CLIP_PCT = 45.0, 1.0, 99.9
EDGES = np.linspace(0.0, 3.5, 71)
CTR = 0.5*(EDGES[1:] + EDGES[:-1])


def ang(a, b):
    return np.degrees(np.arccos(np.clip(np.sum(a*b, axis=1) /
                                        (np.linalg.norm(a, axis=1)*np.linalg.norm(b, axis=1) + 1e-30), -1, 1)))


def load(mode, tag):
    p = f"{EV}/events/{mode}_{tag}_ev3N.root"
    if not os.path.exists(p):
        return None
    a = uproot.open(p)["events"].arrays(
        ["weight", "lead", "recoil2", "recoil3", "q", "Q2", "xB", "scattering_angle"], library="np")
    att = float(uproot.open(f"{EV}/{mode}_{tag}.root")["hAttempts"].values()[0])
    w, Q2, xB, th_e = a["weight"], a["Q2"], a["xB"], a["scattering_angle"]
    p1 = a["lead"][:, :3] - a["q"][:, :3]
    t12 = ang(p1, a["recoil2"][:, :3])
    t23 = ang(a["recoil2"][:, :3], a["recoil3"][:, :3])
    base = np.isfinite(w) & (w > 0) & (th_e < THETA_E_MAX) & (Q2 > Q2_MIN)
    w = np.minimum(w, float(np.percentile(w[base], CLIP_PCT)))
    return dict(w=w, xB=xB, t12=t12, t23=t23, base=base, att=att)


for mode in ("pwia", "fsi"):
    data = {tag: load(mode, tag) for _, tag, _, _ in SAMPLES}
    data = {k: v for k, v in data.items() if v is not None}
    if not data:
        continue
    fig, axes = plt.subplots(2, 3, figsize=(15, 7.5))
    rows = []
    for col, (rname, ireg) in enumerate(REGIONS):
        for lab, tag, c, ls in SAMPLES:
            d = data.get(tag)
            if d is None:
                continue
            m = d["base"] & in_region(d["t12"], d["t23"], ireg)
            h, _ = np.histogram(d["xB"][m], bins=EDGES, weights=d["w"][m])
            # AV8's weight scale is arbitrary (unknown table convention) -> shape only
            if tag != "av8":
                axes[0, col].plot(CTR, h/d["att"], color=c, ls=ls, lw=1.6, label=lab)
            if h.sum() > 0:
                axes[1, col].plot(CTR, h/h.sum(), color=c, ls=ls, lw=1.6, label=lab)
                mean = (h*CTR).sum()/h.sum()
                frac = h[CTR > 1.0].sum()/h.sum()
                rows.append(f"  {mode:4s} {rname:14s} {lab:18s} <xB> = {mean:5.3f}   "
                            f"frac(xB>1) = {100*frac:5.1f} %   sum(w)/att = {h.sum()/d['att']:.4g}")
        for r in (0, 1):
            axes[r, col].set_xlim(0.5, 3.0)
            axes[r, col].set_xlabel(r"$x_B$")
            axes[r, col].set_yscale("log")      # the high-xB tail is 1e-3 of the peak
        axes[0, col].set_title(rname, fontsize=12)
    axes[0, 0].set_ylabel("weight / generation attempt")
    axes[1, 0].set_ylabel("area-normalised")
    axes[0, 0].legend(frameon=False, fontsize=8)
    axes[1, 0].legend(frameon=False, fontsize=8)
    fig.suptitle(f"$x_B$ by 3N configuration, {mode.upper()} — wavefunction comparison "
                 f"($\\theta_e<45^\\circ$, $Q^2>1$, 99.9% weight clip)", fontsize=12)
    fig.tight_layout()
    p = os.path.join(OUT, f"xB_by_config_overlay_{mode}.png")
    fig.savefig(p, dpi=170); plt.close(fig)
    print("saved", p)
    print("\n".join(rows))
    with open(os.path.join(OUT, "xB_by_config_summary.txt"), "a") as fh:
        fh.write("\n".join(rows) + "\n")

    # --- dominance vs xB threshold (this is where the literature statement lives) ---
    thr = np.linspace(0.4, 2.4, 41)
    fig, axes = plt.subplots(1, len(data), figsize=(4.6*len(data), 4.2), squeeze=False)
    for ax, (tag, d) in zip(axes[0], data.items()):
        lab = [l for l, t, _, _ in SAMPLES if t == tag][0]
        tot = np.array([d["w"][d["base"] & (d["xB"] > x)].sum() for x in thr])
        for rname, ireg in REGIONS:
            m = d["base"] & in_region(d["t12"], d["t23"], ireg)
            y = np.array([d["w"][m & (d["xB"] > x)].sum() for x in thr])
            ax.plot(thr, 100*np.divide(y, tot, out=np.zeros_like(y), where=tot > 0), lw=1.8, label=rname)
        ax.axhline(50, color="k", lw=0.6, ls=":")
        ax.set_xlabel(r"$x_B$ threshold"); ax.set_title(lab, fontsize=11); ax.set_ylim(0, 100)
    axes[0][0].set_ylabel("% of weight above threshold")
    axes[0][0].legend(frameon=False, fontsize=9)
    fig.suptitle(f"Which 3N configuration dominates above an $x_B$ cut — {mode.upper()}", fontsize=12)
    fig.tight_layout()
    p2 = os.path.join(OUT, f"xB_dominance_vs_threshold_{mode}.png")
    fig.savefig(p2, dpi=170); plt.close(fig)
    print("saved", p2)
