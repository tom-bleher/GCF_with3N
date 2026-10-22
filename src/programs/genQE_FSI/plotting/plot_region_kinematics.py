#!/usr/bin/env python3
"""
Plot Region L/R kinematic distributions: overlay 2N+FSI vs 3N SRC for each variable.
Reads directly from ROOT TTree files via uproot.

Regions are triangular, defined by parameters A=135 deg, K=4 (as in Fig. 5).
Region R: triangle near (180,180) bounded by lines with slopes B and 1/B through
          (180,180), plus a diagonal with slope -1.
Region L: (x,y) in L iff (360-x-y, y) in R  (theta12 <-> theta13 swap).

For Region L: compares N=2 FSI events (genQE_FSI) with pure 3N SRC (genQE_3N)
For Region R: compares N=3 FSI events (genQE_FSI) with pure 3N SRC (genQE_3N)
"""
import os
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt

# Try to import uproot
try:
    import uproot
except ImportError:
    print("ERROR: uproot not installed. Install with: pip install uproot")
    sys.exit(1)

# ──────────── Physical constants ────────────
mN = 0.93892
kF = 0.25
pCode = 2212
nCode = 2112

# ──────────── Region definitions ────────────
# The paper's regions (A=135 deg, K=4), shared with Fig. 5.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "genQE_3N", "plotting"))
from plot_fig5_theta_heatmap_pair import (region_params,
                                          in_region_R_array as in_region_R,
                                          in_region_L_array as in_region_L)

# ──────────── Vector math ────────────
def angle_deg(ax, ay, az, bx, by, bz):
    dot = ax*bx + ay*by + az*bz
    a_mag = np.sqrt(ax**2 + ay**2 + az**2)
    b_mag = np.sqrt(bx**2 + by**2 + bz**2)
    return np.degrees(np.arccos(np.clip(dot / (a_mag * b_mag + 1e-30), -1, 1)))

# ──────────── Load 2N+FSI data ────────────
def load_fsi(filepath):
    """Load and preprocess 2N+FSI events."""
    f = uproot.open(filepath)
    d = f["events"].arrays(library="np")

    q_px, q_py, q_pz = d["q"][:, 0], d["q"][:, 1], d["q"][:, 2]
    q_mag = np.sqrt(q_px**2 + q_py**2 + q_pz**2)
    lp_px, lp_py, lp_pz = d["lead_post"][:, 0], d["lead_post"][:, 1], d["lead_post"][:, 2]
    rp_px, rp_py, rp_pz = d["recoil_post"][:, 0], d["recoil_post"][:, 1], d["recoil_post"][:, 2]
    lpre_px, lpre_py, lpre_pz = d["lead_pre"][:, 0], d["lead_pre"][:, 1], d["lead_pre"][:, 2]
    rpre_px, rpre_py, rpre_pz = d["recoil_pre"][:, 0], d["recoil_pre"][:, 1], d["recoil_pre"][:, 2]
    e_px, e_py, e_pz, e_E = d["electron"][:, 0], d["electron"][:, 1], d["electron"][:, 2], d["electron"][:, 3]

    # p1 = pmiss
    p1_x, p1_y, p1_z = lp_px - q_px, lp_py - q_py, lp_pz - q_pz
    p1_mag = np.sqrt(p1_x**2 + p1_y**2 + p1_z**2)
    p1_after_mag = np.sqrt(lp_px**2 + lp_py**2 + lp_pz**2)
    p2_mag = np.sqrt(rp_px**2 + rp_py**2 + rp_pz**2)

    # p3_hyp = -(p1 + p2)
    p3h_x, p3h_y, p3h_z = -(p1_x + rp_px), -(p1_y + rp_py), -(p1_z + rp_pz)

    d["theta12"] = angle_deg(p1_x, p1_y, p1_z, rp_px, rp_py, rp_pz)
    d["theta23"] = angle_deg(rp_px, rp_py, rp_pz, p3h_x, p3h_y, p3h_z)
    d["p1_mag"] = p1_mag
    d["p1_after_mag"] = p1_after_mag
    d["p2_mag"] = p2_mag
    d["p1_after_angle_q"] = angle_deg(lp_px, lp_py, lp_pz, q_px, q_py, q_pz)
    d["lead_over_q"] = p1_after_mag / (q_mag + 1e-30)
    d["e_angle_q"] = angle_deg(e_px, e_py, e_pz, q_px, q_py, q_pz)
    d["e_mom"] = e_E
    d["p2_angle_q"] = angle_deg(rp_px, rp_py, rp_pz, q_px, q_py, q_pz)
    d["p1_angle_q"] = angle_deg(p1_x, p1_y, p1_z, q_px, q_py, q_pz)

    # Pre-FSI
    d["p1_pre_angle_q"] = angle_deg(lpre_px - q_px, lpre_py - q_py, lpre_pz - q_pz, q_px, q_py, q_pz)
    d["p2_pre_angle_q"] = angle_deg(rpre_px, rpre_py, rpre_pz, q_px, q_py, q_pz)
    p3pre_x = -((lpre_px - q_px) + rpre_px)
    p3pre_y = -((lpre_py - q_py) + rpre_py)
    p3pre_z = -((lpre_pz - q_pz) + rpre_pz)
    d["p3_pre_angle_q"] = angle_deg(p3pre_x, p3pre_y, p3pre_z, q_px, q_py, q_pz)

    # FSI secondaries: find p3_fsi for N=3
    n = len(d["weight"])
    p3f_px, p3f_py, p3f_pz = np.zeros(n), np.zeros(n), np.zeros(n)
    has_p3 = np.zeros(n, dtype=bool)
    nSec, sec_pdg = d["nSec"], d["sec_pdg"]
    sec_px, sec_py, sec_pz = d["sec_px"], d["sec_py"], d["sec_pz"]
    for i in range(n):
        best = -1.0
        for j in range(nSec[i]):
            pdg = sec_pdg[i][j]
            if pdg == pCode or pdg == nCode:
                mag = np.sqrt(sec_px[i][j]**2 + sec_py[i][j]**2 + sec_pz[i][j]**2)
                if mag > kF and mag > best:
                    best = mag
                    p3f_px[i], p3f_py[i], p3f_pz[i] = sec_px[i][j], sec_py[i][j], sec_pz[i][j]
                    has_p3[i] = True
    d["p3_fsi_mag"] = np.sqrt(p3f_px**2 + p3f_py**2 + p3f_pz**2)
    d["has_p3_fsi"] = has_p3
    d["theta12_n3"] = angle_deg(p1_x, p1_y, p1_z, rp_px, rp_py, rp_pz)
    d["theta23_n3"] = angle_deg(rp_px, rp_py, rp_pz, p3f_px, p3f_py, p3f_pz)
    d["p3_fsi_angle_q"] = angle_deg(p3f_px, p3f_py, p3f_pz, q_px, q_py, q_pz)

    return d


def load_3n(filepath):
    """Load and preprocess 3N SRC events."""
    f = uproot.open(filepath)
    d = f["events"].arrays(library="np")

    q_px, q_py, q_pz = d["q"][:, 0], d["q"][:, 1], d["q"][:, 2]
    q_mag = np.sqrt(q_px**2 + q_py**2 + q_pz**2)
    l_px, l_py, l_pz = d["lead"][:, 0], d["lead"][:, 1], d["lead"][:, 2]
    r2_px, r2_py, r2_pz = d["recoil2"][:, 0], d["recoil2"][:, 1], d["recoil2"][:, 2]
    r3_px, r3_py, r3_pz = d["recoil3"][:, 0], d["recoil3"][:, 1], d["recoil3"][:, 2]
    e_px, e_py, e_pz, e_E = d["electron"][:, 0], d["electron"][:, 1], d["electron"][:, 2], d["electron"][:, 3]

    p1_x, p1_y, p1_z = l_px - q_px, l_py - q_py, l_pz - q_pz
    d["p1_mag"] = np.sqrt(p1_x**2 + p1_y**2 + p1_z**2)
    d["p1_after_mag"] = np.sqrt(l_px**2 + l_py**2 + l_pz**2)
    d["p2_mag"] = np.sqrt(r2_px**2 + r2_py**2 + r2_pz**2)
    d["p3_mag"] = np.sqrt(r3_px**2 + r3_py**2 + r3_pz**2)

    d["theta12"] = angle_deg(p1_x, p1_y, p1_z, r2_px, r2_py, r2_pz)
    d["theta23"] = angle_deg(r2_px, r2_py, r2_pz, r3_px, r3_py, r3_pz)
    d["p1_after_angle_q"] = angle_deg(l_px, l_py, l_pz, q_px, q_py, q_pz)
    d["lead_over_q"] = d["p1_after_mag"] / (q_mag + 1e-30)
    d["e_angle_q"] = angle_deg(e_px, e_py, e_pz, q_px, q_py, q_pz)
    d["e_mom"] = e_E
    d["p1_angle_q"] = angle_deg(p1_x, p1_y, p1_z, q_px, q_py, q_pz)
    d["p2_angle_q"] = angle_deg(r2_px, r2_py, r2_pz, q_px, q_py, q_pz)
    d["p3_angle_q"] = angle_deg(r3_px, r3_py, r3_pz, q_px, q_py, q_pz)

    return d


# ──────────── Overlay plotting ────────────
VARIABLES = [
    ("init_lead_angle_q",  r"Initial lead angle with $\vec{q}$", "deg", (0, 180)),
    ("init_rec1_angle_q",  r"Initial recoil 1 angle with $\vec{q}$", "deg", (0, 180)),
    ("init_rec2_angle_q",  r"Initial recoil 2 angle with $\vec{q}$", "deg", (0, 180)),
    ("final_lead_angle_q", r"Final lead angle with $\vec{q}$", "deg", (0, 50)),
    ("final_rec1_angle_q", r"Final recoil 1 angle with $\vec{q}$", "deg", (0, 180)),
    ("final_rec2_angle_q", r"Final recoil 2 angle with $\vec{q}$", "deg", (0, 180)),
    ("pmiss_LR",           r"$p_{\mathrm{miss}}$", "GeV/c", (0, 1.5)),
    ("pmiss_angle_q",      r"$p_{\mathrm{miss}}$ angle with $\vec{q}$", "deg", (0, 180)),
    ("p2_final_mom",       r"$|p_2|$ (final)", "GeV/c", (0, 1.5)),
    ("p3_final_mom",       r"$|p_3|$ (final)", "GeV/c", (0, 1.5)),
    ("e_angle_q_LR",       r"Electron angle with $\vec{q}$", "deg", (0, 90)),
    ("e_mom_LR",           r"Electron momentum", "GeV/c", (0, 8)),
    ("xB_LR",              r"$x_B$", "", (0, 4)),
    ("Q2_LR",              r"$Q^2$", r"GeV$^2$", (2, 10)),
    ("lead_mom_over_q",    r"$|p_{\mathrm{lead}}|/|\vec{q}|$", "", (0, 2)),
    ("p1_final_mom",       r"$|p_1|$ (final)", "GeV/c", (0, 8)),
]

# Map variable names to data keys
def get_fsi_n2_val(name, d):
    m = {
        "init_lead_angle_q": d["p1_pre_angle_q"], "init_rec1_angle_q": d["p2_pre_angle_q"],
        "init_rec2_angle_q": d["p3_pre_angle_q"],
        "final_lead_angle_q": d["p1_after_angle_q"], "final_rec1_angle_q": d["p2_angle_q"],
        "pmiss_LR": d["pmiss"], "pmiss_angle_q": d["p1_angle_q"],
        "p2_final_mom": d["p2_mag"],
        "e_angle_q_LR": d["e_angle_q"], "e_mom_LR": d["e_mom"],
        "xB_LR": d["xB"], "Q2_LR": d["Q2"],
        "lead_mom_over_q": d["lead_over_q"], "p1_final_mom": d["p1_after_mag"],
    }
    return m.get(name)

def get_fsi_n3_val(name, d):
    m = {
        "init_lead_angle_q": d["p1_pre_angle_q"], "init_rec1_angle_q": d["p2_pre_angle_q"],
        "final_lead_angle_q": d["p1_after_angle_q"], "final_rec1_angle_q": d["p2_angle_q"],
        "final_rec2_angle_q": d["p3_fsi_angle_q"],
        "pmiss_LR": d["pmiss"], "pmiss_angle_q": d["p1_angle_q"],
        "p2_final_mom": d["p2_mag"], "p3_final_mom": d["p3_fsi_mag"],
        "e_angle_q_LR": d["e_angle_q"], "e_mom_LR": d["e_mom"],
        "xB_LR": d["xB"], "Q2_LR": d["Q2"],
        "lead_mom_over_q": d["lead_over_q"], "p1_final_mom": d["p1_after_mag"],
    }
    return m.get(name)

def get_3n_val(name, d):
    m = {
        "init_lead_angle_q": d["p1_angle_q"], "init_rec1_angle_q": d["p2_angle_q"],
        "init_rec2_angle_q": d["p3_angle_q"],
        "final_lead_angle_q": d["p1_after_angle_q"], "final_rec1_angle_q": d["p2_angle_q"],
        "final_rec2_angle_q": d["p3_angle_q"],
        "pmiss_LR": d["pmiss"], "pmiss_angle_q": d["p1_angle_q"],
        "p2_final_mom": d["p2_mag"], "p3_final_mom": d["p3_mag"],
        "e_angle_q_LR": d["e_angle_q"], "e_mom_LR": d["e_mom"],
        "xB_LR": d["xB"], "Q2_LR": d["Q2"],
        "lead_mom_over_q": d["lead_over_q"], "p1_final_mom": d["p1_after_mag"],
    }
    return m.get(name)


def plot_1d_norm(ax, vals, weights, bins, range_, **kwargs):
    counts, edges = np.histogram(vals, bins=bins, range=range_, weights=weights)
    centers = 0.5 * (edges[:-1] + edges[1:])
    s = np.sum(counts)
    if s > 0:
        counts = counts / s
    ax.step(centers, counts, where='mid', **kwargs)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fsi", default=None, help="Path to 2N events ROOT file")
    parser.add_argument("--src3n", default=None, help="Path to events_3N.root")
    parser.add_argument("-o", "--output-dir", default="region_overlays", help="Output directory")
    args = parser.parse_args()

    # Auto-detect paths relative to script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    if args.fsi is None:
        args.fsi = os.path.join(script_dir, '..', 'events', 'hN', 'events_2N_fsi_501.root')
    if args.src3n is None:
        args.src3n = os.path.join(script_dir, '..', '..', 'genQE_3N', 'events_3N.root')

    os.makedirs(args.output_dir, exist_ok=True)

    have_fsi = os.path.exists(args.fsi)
    have_3n = os.path.exists(args.src3n)

    if have_fsi:
        print(f"Loading FSI: {args.fsi}")
        fsi = load_fsi(args.fsi)
        print(f"  {len(fsi['weight'])} events")
    else:
        print(f"FSI file not found: {args.fsi}")
        fsi = None

    if have_3n:
        print(f"Loading 3N:  {args.src3n}")
        src3n = load_3n(args.src3n)
        print(f"  {len(src3n['weight'])} events")
    else:
        print(f"3N file not found: {args.src3n}")
        src3n = None

    if fsi is None and src3n is None:
        print("No data to plot!")
        return

    # Build masks
    if fsi is not None:
        # N=2 cuts
        m_n2 = fsi["nAboveKF"] == 2
        m_n2_cuts = m_n2 & (fsi["scattering_angle"] < 45) & (fsi["p1_after_angle_q"] < 10) & (fsi["lead_over_q"] > 0.7)
        m_n2_L = m_n2_cuts & in_region_L(fsi["theta12"], fsi["theta23"])
        m_n2_R = m_n2_cuts & in_region_R(fsi["theta12"], fsi["theta23"])

        # N=3 cuts
        m_n3 = (fsi["nAboveKF"] >= 3) & fsi["has_p3_fsi"]
        m_n3_cuts = m_n3 & (fsi["scattering_angle"] < 45) & (fsi["p1_after_angle_q"] < 10) & (fsi["lead_over_q"] > 0.7) & (fsi["p2_mag"] > 0.5)
        m_n3_L = m_n3_cuts & in_region_L(fsi["theta12_n3"], fsi["theta23_n3"])
        m_n3_R = m_n3_cuts & in_region_R(fsi["theta12_n3"], fsi["theta23_n3"])

    if src3n is not None:
        m_3n_cuts = (src3n["p1_mag"] > kF) & (src3n["p2_mag"] > kF) & (src3n["p3_mag"] > kF) \
                    & (src3n["scattering_angle"] < 45) & (src3n["p1_after_angle_q"] < 10) \
                    & (src3n["lead_over_q"] > 0.7) & (src3n["p2_mag"] > 0.5)
        m_3n_L = m_3n_cuts & in_region_L(src3n["theta12"], src3n["theta23"])
        m_3n_R = m_3n_cuts & in_region_R(src3n["theta12"], src3n["theta23"])

    _, B = region_params()
    print(f"\nRegion R: triangular near (180,180), A=135 deg, K=4, B={B:.4f}")
    print(f"Region L: (x,y) in L iff (360-x-y, y) in R\n")

    # ── Region L: N=2 FSI vs 3N SRC ──
    print("=== Region L ===")
    for var_name, var_label, unit, range_ in VARIABLES:
        fig, ax = plt.subplots(figsize=(7, 5))
        plotted = False

        if fsi is not None and np.sum(m_n2_L) > 0:
            vals = get_fsi_n2_val(var_name, fsi)
            if vals is not None:
                plot_1d_norm(ax, vals[m_n2_L], fsi["weight"][m_n2_L], 45, range_,
                             label="2N+FSI N=2", linewidth=1.5, color='tab:blue')
                plotted = True

        if src3n is not None and np.sum(m_3n_L) > 0:
            vals = get_3n_val(var_name, src3n)
            if vals is not None:
                plot_1d_norm(ax, vals[m_3n_L], src3n["weight"][m_3n_L], 45, range_,
                             label="3N SRC", linewidth=1.5, color='tab:red')
                plotted = True

        if not plotted:
            plt.close(fig)
            continue

        region_desc = r"Region L (triangular, $A=135°$, $K=4$)"
        ax.set_title(f"{var_label} — {region_desc}", fontsize=13)
        xlabel = var_label + (f" [{unit}]" if unit else "")
        ax.set_xlabel(xlabel, fontsize=12)
        ax.set_ylabel("Normalized weight", fontsize=12)
        ax.legend(fontsize=11)
        fig.tight_layout()
        fname = f"region_{var_name}_L.png"
        fig.savefig(os.path.join(args.output_dir, fname), dpi=200)
        plt.close(fig)
        print(f"  {fname}")

    # ── Region R: N=3 FSI vs 3N SRC ──
    print("\n=== Region R ===")
    for var_name, var_label, unit, range_ in VARIABLES:
        fig, ax = plt.subplots(figsize=(7, 5))
        plotted = False

        if fsi is not None and np.sum(m_n3_R) > 0:
            vals = get_fsi_n3_val(var_name, fsi)
            if vals is not None:
                plot_1d_norm(ax, vals[m_n3_R], fsi["weight"][m_n3_R], 45, range_,
                             label="2N+FSI N=3", linewidth=1.5, color='tab:blue')
                plotted = True

        if src3n is not None and np.sum(m_3n_R) > 0:
            vals = get_3n_val(var_name, src3n)
            if vals is not None:
                plot_1d_norm(ax, vals[m_3n_R], src3n["weight"][m_3n_R], 45, range_,
                             label="3N SRC", linewidth=1.5, color='tab:red')
                plotted = True

        if not plotted:
            plt.close(fig)
            continue

        region_desc = r"Region R (triangular, $A=135°$, $K=4$)"
        ax.set_title(f"{var_label} — {region_desc}", fontsize=13)
        xlabel = var_label + (f" [{unit}]" if unit else "")
        ax.set_xlabel(xlabel, fontsize=12)
        ax.set_ylabel("Normalized weight", fontsize=12)
        ax.legend(fontsize=11)
        fig.tight_layout()
        fname = f"region_{var_name}_R.png"
        fig.savefig(os.path.join(args.output_dir, fname), dpi=200)
        plt.close(fig)
        print(f"  {fname}")

    # ── Region fractions table ──
    print("\n--- Region fractions with kinematic cuts ---")
    print(f"{'Source':<20s} | {'Region L frac':>14s} | {'Region R frac':>14s}")
    print("-" * 55)
    if fsi is not None:
        w_n2 = np.sum(fsi["weight"][m_n2_cuts])
        if w_n2 > 0:
            print(f"{'2N+FSI  N=2':<20s} | {100*np.sum(fsi['weight'][m_n2_L])/w_n2:>13.4f}% | {100*np.sum(fsi['weight'][m_n2_R])/w_n2:>13.4f}%")
        w_n3 = np.sum(fsi["weight"][m_n3_cuts])
        if w_n3 > 0:
            print(f"{'2N+FSI  N=3':<20s} | {100*np.sum(fsi['weight'][m_n3_L])/w_n3:>13.4f}% | {100*np.sum(fsi['weight'][m_n3_R])/w_n3:>13.4f}%")
    if src3n is not None:
        w_3n = np.sum(src3n["weight"][m_3n_cuts])
        if w_3n > 0:
            print(f"{'3N generator':<20s} | {100*np.sum(src3n['weight'][m_3n_L])/w_3n:>13.4f}% | {100*np.sum(src3n['weight'][m_3n_R])/w_3n:>13.4f}%")

    print("\nDone.")


if __name__ == '__main__':
    main()
