#!/usr/bin/env python3
"""Fake-3N opening-angle heatmap from the single-nucleon mean-field generator.

Physics question: can a mean-field (MF) event that ejects a second energetic
nucleon via FSI fake a 3N-SRC signal?

We take MF events (SRC_analysis_2N wf_mode=4: one struck nucleon from a Fermi
gas, GENIE FSI in the A-1 residual, no correlated recoil) in the (e,e'pp)
channel — exactly TWO detected protons above kF, the SAME detection-based
rule as the 3N and 2N paper figures — reconstruct a hypothetical third
nucleon assuming the initial three-nucleon system had zero center-of-mass
momentum, and plot the theta12-theta23 opening-angle heatmap.

Conventions (unified events schema from SRC_analysis_2N):
  * lead = detected proton closest in angle to q; recoil = the other one
  * p1 = lead - q               (reconstructed initial struck-nucleon momentum)
  * p2 = recoil                 (final momentum of the second detected proton)
  * p3 = -(p1 + p2)             (zero initial CM: p1 + p2 + p3 = 0)
  * theta12 = angle(p1, p2),  theta23 = angle(p2, p3)
Region L/R/BR triangles (A=135 deg, K=4.0) and the equal-area center circle are
overlaid, and the weighted fraction of events in each region is written out.

Usage:
    python plot_mf_fake3n_heatmap.py [--input FILE] [--output FILE]
        [--src-cuts] [--interplane-max 20] [--symmetrize]
"""
import argparse
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Polygon
import uproot


# Paper styling shared with the genQE_3N theta-heatmap plots.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "genQE_3N", "plotting"))
from paper_style import apply_style, figure_size, CMAP, OVERLAY_COLOR

# ──────────── Constants ────────────
KF_DEFAULT = 0.25   # Fermi momentum [GeV/c]
P_CODE, N_CODE = 2212, 2112


# ──────────── Region definitions (A=135 deg, K=4.0) ────────────
def region_params(A=135.0, K=4.0):
    A_rad = np.radians(A)
    B = np.degrees(np.arctan(np.sin(A_rad) / (K + np.cos(A_rad)))) / (180.0 - A)
    return A, B

def in_region_R(t12, t23, A=135.0, K=4.0):
    _, B = region_params(A, K)
    line1 = B * (t12 - 180.0) + 180.0
    line2 = (1.0 / B) * (t12 - 180.0) + 180.0
    line3 = -(t12 - A) + 180.0 - (180.0 - A) * B
    return (t23 <= line1) & (t23 >= line2) & (t23 >= line3)

def in_region_L(t12, t23, A=135.0, K=4.0):
    return in_region_R(360.0 - t12 - t23, t23, A, K)

def in_region_BR(t12, t23, A=135.0, K=4.0):
    return in_region_R(t12, 360.0 - t12 - t23, A, K)

def in_region_center(t12, t23, radius):
    return (t12 - 120.0) ** 2 + (t23 - 120.0) ** 2 < radius ** 2

def _triangle_area(p1, p2, p3):
    return 0.5 * abs((p2[0] - p1[0]) * (p3[1] - p1[1]) -
                     (p3[0] - p1[0]) * (p2[1] - p1[1]))

def get_triangle_vertices():
    """Vertices of the R and L region triangles (from plot_all_heatmaps.py)."""
    A, B = region_params()
    def _li(m1, b1, m2, b2):
        d = m1 - m2
        if abs(d) < 1e-15:
            return None
        x = (b2 - b1) / d
        return (float(x), float(m1 * x + b1))
    m1_R, b1_R = B, 180.0 * (1.0 - B)
    m2_R, b2_R = (1.0 / B), 180.0 * (1.0 - (1.0 / B))
    m3_R, b3_R = -1.0, (A + 180.0 - (180.0 - A) * B)
    p12_R = _li(m1_R, b1_R, m2_R, b2_R)
    p13_R = _li(m1_R, b1_R, m3_R, b3_R)
    p23_R = _li(m2_R, b2_R, m3_R, b3_R)

    s1_L = -(m1_R / (1.0 + m1_R))
    s2_L = -(m2_R / (1.0 + m2_R))
    x_vert_L = 360.0 - b3_R
    p12_L = (0.0, 180.0)
    p13_L = (x_vert_L, s1_L * x_vert_L + 180.0)
    p23_L = (x_vert_L, s2_L * x_vert_L + 180.0)
    return (p12_R, p13_R, p23_R), (p12_L, p13_L, p23_L)

def add_triangle(ax, points, *, edgecolor, label):
    pts = [p for p in points if p is not None and np.isfinite(p[0]) and np.isfinite(p[1])]
    if len(pts) != 3:
        return
    cx = sum(p[0] for p in pts) / 3.0
    cy = sum(p[1] for p in pts) / 3.0
    pts = sorted(pts, key=lambda p: np.arctan2(p[1] - cy, p[0] - cx))
    ax.add_patch(Polygon(pts, closed=True, fill=False,
                         edgecolor=edgecolor, linewidth=2.0, label=label))


# ──────────── Angle helpers ────────────
def angle_deg(a, b):
    """Opening angle [deg] between arrays of 3-vectors a, b (shape (N,3))."""
    dot = np.sum(a * b, axis=1)
    am = np.linalg.norm(a, axis=1)
    bm = np.linalg.norm(b, axis=1)
    return np.degrees(np.arccos(np.clip(dot / (am * bm + 1e-30), -1, 1)))

def interplane_angle(shared, a, b):
    """Angle [deg] between planes {shared,a} and {shared,b}."""
    n1 = np.cross(shared, a)
    n2 = np.cross(shared, b)
    n1m = np.linalg.norm(n1, axis=1, keepdims=True)
    n2m = np.linalg.norm(n2, axis=1, keepdims=True)
    cos_a = np.clip(np.abs(np.sum(n1 / (n1m + 1e-30) * n2 / (n2m + 1e-30), axis=1)), 0, 1)
    return np.degrees(np.arccos(cos_a))


# ──────────── Load + select + reconstruct ────────────
def build_fake3n(path, ebeam, kf, require_pp=True):
    """Load detection-based MF `events` (SRC_analysis_2N unified schema),
    select the (e,e'pp) channel — exactly 2 detected protons above kF, the
    SAME rule as the 3N and 2N figures — reconstruct the zero-CM third
    nucleon, and return a dict of per-event arrays for the selected events
    (plus derived SRC-cut quantities).

    The tree already contains the detection-based lead (detected proton
    closest to q) and recoil (highest-momentum other detected proton);
    require_pp=False relaxes the channel to any 2 detected nucleons
    (nNucleonAboveKF == 2)."""
    tree = uproot.open(path)["events"]
    flat = tree.arrays(["weight", "lead", "recoil", "q",
                        "Q2", "scattering_angle", "xB",
                        "nProtonAboveKF", "nNucleonAboveKF"], library="np")
    w_all = flat["weight"]
    w_all_file = float(w_all[np.isfinite(w_all) & (w_all > 0)].sum())
    n_all_file = len(w_all)

    n_tot = n_all_file
    lead = flat["lead"][:, :3]
    rec = flat["recoil"][:, :3]
    qv = flat["q"][:, :3]

    if require_pp:
        sel = flat["nProtonAboveKF"] == 2
        chan = "(e,e'pp): exactly 2 detected protons above kF"
    else:
        sel = flat["nNucleonAboveKF"] == 2
        chan = "pN: exactly 2 detected nucleons above kF"

    print(f"  {n_tot} events")
    print(f"  selection [{chan}]: {int(sel.sum())} events")

    p1 = lead - qv                  # reconstructed initial struck-nucleon momentum
    p2 = rec                        # second detected proton (final momentum)
    p3 = -(p1 + p2)                 # zero initial CM:  p1 + p2 + p3 = 0

    d = {}
    d["weight"] = flat["weight"][sel]
    d["p1"] = p1[sel]
    d["p2"] = p2[sel]
    d["p3"] = p3[sel]
    d["theta12"] = angle_deg(d["p1"], d["p2"])
    d["theta23"] = angle_deg(d["p2"], d["p3"])
    d["theta13"] = angle_deg(d["p1"], d["p3"])
    d["scattering_angle"] = flat["scattering_angle"][sel]
    d["Q2"] = flat["Q2"][sel]
    d["xB"] = flat["xB"][sel]
    # SRC-cut quantities
    d["pmiss"] = np.linalg.norm(d["p1"], axis=1)
    pN = np.linalg.norm(lead[sel], axis=1)
    q3 = np.linalg.norm(qv[sel], axis=1)
    d["pN_over_q"] = pN / (q3 + 1e-30)
    d["theta_pq"] = angle_deg(lead[sel], qv[sel])
    d["interplane"] = interplane_angle(d["p1"], d["p2"], d["p3"])
    # zero-CM sanity
    resid = np.linalg.norm(d["p1"] + d["p2"] + d["p3"], axis=1)
    print(f"  zero-CM check: max|p1+p2+p3| = {resid.max():.2e} GeV/c")
    d["_w_all_file"] = w_all_file
    d["_n_all_file"] = n_all_file
    return d


# ──────────── Plot ────────────
def plot(d, args):
    w = d["weight"]
    mask_base = (np.isfinite(w) & (w > 0)
                 & (d["scattering_angle"] > args.theta_e_min)
                 & (d["scattering_angle"] < args.theta_e_max)
                 & (d["Q2"] > args.Q2_min))
    cut_desc = f"{args.theta_e_min}<scattering<{args.theta_e_max}, Q2>{args.Q2_min}"
    if args.xB_max is not None:
        mask_base &= (d["xB"] < args.xB_max)
        cut_desc += f", xB<{args.xB_max}"
    mask_src = (mask_base
                & (d["theta_pq"] < 8.0) & (d["pN_over_q"] > 0.75)
                & (d["pmiss"] > 0.25) & (d["pmiss"] < 0.9) & (d["xB"] < 1.2))

    mask = mask_base.copy()
    if args.src_cuts:
        mask = mask_src.copy()
        cut_desc += ", SRC(thpq<8,pN/q>0.75,0.25<pmiss<0.9,xB<1.2)"
    if args.interplane_max is not None:
        mask &= (d["interplane"] < args.interplane_max)
        cut_desc += f", interplane<{args.interplane_max}"

    t12, t23 = d["theta12"][mask], d["theta23"][mask]
    ww = w[mask]
    print(f"  events after cuts: {int(mask.sum())}  ({cut_desc})")

    if args.symmetrize:
        # 2<->3 swap (cf. plot_all_heatmaps '3detected'): also fill (theta13, theta23).
        t12 = np.concatenate([t12, d["theta13"][mask]])
        t23 = np.concatenate([t23, d["theta23"][mask]])
        ww = np.concatenate([ww, ww])

    h, xe, ye = np.histogram2d(t12, t23, bins=args.bins,
                               range=[[0, 180], [0, 180]], weights=ww)
    h_norm = h / (h.sum() + 1e-30)
    h_plot = np.where(h_norm > 0, h_norm, np.nan)
    pos = h_norm[h_norm > 0]
    vmax = args.vmax if args.vmax is not None else (pos.max() if pos.size else 1e-2)
    vmin = args.vmin if args.vmin is not None else vmax * args.vmin_frac

    # Region geometry (integration only; nothing is drawn on the plot).
    tri_R, tri_L = get_triangle_vertices()
    p12_R, p13_R, p23_R = tri_R
    tri_area = (_triangle_area(p12_R, p13_R, p23_R)
                if (p12_R and p13_R and p23_R) else 0.0)
    center_radius = np.sqrt(tri_area / np.pi) if tri_area > 0 else 10.0

    apply_style()
    fig, ax = plt.subplots(figsize=figure_size(cols=1, ratio=0.92))
    im = ax.pcolormesh(xe, ye, h_plot.T, cmap=CMAP,
                       norm=LogNorm(vmin=vmin, vmax=vmax),
                       shading="auto", rasterized=True)
    cbar = fig.colorbar(im, ax=ax, label="Normalized weight",
                        pad=0.02, fraction=0.046)
    cbar.ax.tick_params(labelsize=8)
    cbar.ax.yaxis.label.set_size(9)

    x = np.linspace(0, 180, 200)
    line_kw = dict(color=OVERLAY_COLOR, linestyle="--", lw=0.7, alpha=0.85)
    ax.plot(x, 180 - x / 2, **line_kw)
    ax.plot(180 - x / 2, x, **line_kw)
    ax.plot(x, x, **line_kw)

    ax.set_xlim(0, 180)
    ax.set_ylim(0, 180)
    ax.set_aspect("equal")
    ax.set_xticks(np.arange(0, 181, 30))
    ax.set_yticks(np.arange(0, 181, 30))
    ax.set_xlabel(r"$\theta_{12}$ (deg)")
    ax.set_ylabel(r"$\theta_{23}$ (deg)")
    ax.set_title(r"MF+FSI on $^{12}$C ($n_{>k_F}=2$, $N=%d$)" % int(mask.sum()),
                 fontsize=8)

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    fig.tight_layout()
    fig.savefig(args.output, dpi=args.dpi)
    plt.close(fig)

    # Region fractions (weighted), on the un-symmetrized selection.
    t12c, t23c, wc = d["theta12"][mask], d["theta23"][mask], w[mask]
    wtot = wc.sum() + 1e-30
    fr = dict(
        L=wc[in_region_L(t12c, t23c)].sum() / wtot,
        R=wc[in_region_R(t12c, t23c)].sum() / wtot,
        BR=wc[in_region_BR(t12c, t23c)].sum() / wtot,
        center=wc[in_region_center(t12c, t23c, center_radius)].sum() / wtot,
    )
    LR_ratio = fr["L"] / fr["R"] if fr["R"] > 0 else float("inf")
    # Pass fractions relative to ALL generated events in the file (weighted);
    # the eq2 detector tag from build_fake3n is folded into both masks.
    f_base = w[mask_base].sum() / (d["_w_all_file"] + 1e-30)
    f_src = w[mask_src].sum() / (d["_w_all_file"] + 1e-30)
    frac_path = os.path.splitext(os.path.abspath(args.output))[0] + "_region_fractions.txt"
    with open(frac_path, "w") as fh:
        fh.write(f"# MF+FSI fake-3N region fractions (weighted)\n")
        fh.write(f"# cuts: {cut_desc}\n")
        fh.write(f"# n_events={int(mask.sum())}  center_radius={center_radius:.3f} deg\n")
        for k, v in fr.items():
            fh.write(f"{k:8s} {v:.6f}\n")
        fh.write(f"L/R      {LR_ratio:.6f}\n")
        fh.write(f"# pass fractions vs all {d['_n_all_file']} generated events (weighted)\n")
        fh.write(f"base     {f_base:.6e}\n")
        fh.write(f"src      {f_src:.6e}\n")
    print(f"  region fractions: L={fr['L']:.4f} R={fr['R']:.4f} "
          f"BR={fr['BR']:.4f} center={fr['center']:.4f}  L/R={LR_ratio:.3f}")
    print(f"  base-cut pass fraction (weighted, incl. eq2 tag): {f_base:.4e}  "
          f"({int(mask_base.sum())}/{d['_n_all_file']} events)")
    print(f"  SRC-cut  pass fraction (weighted, incl. eq2 tag): {f_src:.4e}  "
          f"({int(mask_src.sum())}/{d['_n_all_file']} events)")
    print(f"Saved {args.output}")
    print(f"Saved {frac_path}")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", default="events/hA/events_2N_fsi_hA_501_MF.root",
                   help="MF FSI ROOT file (SRC_analysis_2N wf_mode=4)")
    p.add_argument("--output",
                   default="analysis/plots/mf_fake3n/mf_fake3n_heatmap.pdf")
    p.add_argument("--ebeam", type=float, default=5.01)
    p.add_argument("--kf", type=float, default=KF_DEFAULT)
    p.add_argument("--bins", type=int, default=120)
    p.add_argument("--dpi", type=int, default=300)
    p.add_argument("--theta-e-max", type=float, default=45.0)
    p.add_argument("--theta-e-min", type=float, default=0.0)
    p.add_argument("--Q2-min", type=float, default=1.0)
    p.add_argument("--xB-max", type=float, default=None,
                   help="if set, apply xB < this as a BASE cut (it is always "
                        "part of the SRC cuts)")
    p.add_argument("--vmin", type=float, default=None)
    p.add_argument("--vmin-frac", type=float, default=1e-4)
    p.add_argument("--vmax", type=float, default=None)
    p.add_argument("--src-cuts", action="store_true",
                   help="apply 3N SRC cuts (thpq<8, pN/q>0.75, 0.25<pmiss<0.9, xB<1.2)")
    p.add_argument("--interplane-max", type=float, default=None,
                   help="if set, require interplane angle < this (deg), e.g. 20")
    p.add_argument("--symmetrize", action="store_true",
                   help="also fill the 2<->3 swapped (theta13,theta23) entries")
    p.add_argument("--all-nucleons", action="store_true",
                   help="count neutrons too (default: require final pp — both "
                        "nucleons above kF are protons)")
    args = p.parse_args()

    print(f"Reading {args.input}")
    d = build_fake3n(args.input, args.ebeam, args.kf,
                     require_pp=not args.all_nucleons)
    plot(d, args)
    print("Done.")


if __name__ == "__main__":
    main()
