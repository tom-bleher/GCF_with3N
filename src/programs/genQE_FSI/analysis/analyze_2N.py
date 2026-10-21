#!/usr/bin/env python3
"""
Comprehensive analysis of 2N SRC events with FSI.
Reads events_2N.root via uproot, applies cuts, produces all plots and tables.

Usage:
    python analyze_2N.py events_2N.root [-o output_dir]
"""
import os
import sys
import argparse
import numpy as np
import uproot
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Polygon

# ── Publication-quality style ──
mpl.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['CMU Serif', 'Computer Modern Roman', 'Times New Roman', 'DejaVu Serif'],
    'mathtext.fontset': 'cm',
    'axes.linewidth': 1.4,
    'xtick.major.width': 1.2,
    'xtick.minor.width': 0.8,
    'ytick.major.width': 1.2,
    'ytick.minor.width': 0.8,
    'xtick.major.size': 6,
    'xtick.minor.size': 3.5,
    'ytick.major.size': 6,
    'ytick.minor.size': 3.5,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.top': True,
    'ytick.right': True,
    'xtick.minor.visible': True,
    'ytick.minor.visible': True,
    'axes.labelsize': 14,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 11,
    'legend.framealpha': 0.9,
    'lines.linewidth': 1.8,
})

# ──────────── Physical constants ────────────
mN = 0.93892  # nucleon mass GeV
kF = 0.25     # Fermi momentum GeV/c

# ──────────── Region definitions ────────────
# One shared definition (A=135 deg, K=4, as in the paper's Fig. 5): the
# triangles drawn below and the regions counted must be the same.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "..", "genQE_3N", "plotting"))
from plot_fig5_theta_heatmap_pair import (region_params,
                                          in_region_R_array as in_region_R,
                                          in_region_L_array as in_region_L)

def in_region_BR(theta12, theta23, A=135.0, K=4.0):
    return in_region_R(theta12, 360.0 - theta12 - theta23, A, K)

def _triangle_area(p1, p2, p3):
    """Area of triangle given 3 vertices."""
    return 0.5 * abs((p2[0]-p1[0])*(p3[1]-p1[1]) - (p3[0]-p1[0])*(p2[1]-p1[1]))

def in_region_center(theta12, theta23, radius):
    """Circle centered at (120, 120)."""
    return (theta12 - 120.0)**2 + (theta23 - 120.0)**2 < radius**2

# def in_region_R(theta12, theta23):
#     return (theta12 > 135) & (theta23 > 135)

# def in_region_L(theta12, theta23):
#     return (theta12 < 45) & (theta23 > 135)


# ──────────── Data loading ────────────
def load_data(filepath):
    """Load TTree into dict of numpy arrays."""
    f = uproot.open(filepath)
    tree = f["events"]
    # Read all scalar branches as numpy, arrays as numpy
    data = tree.arrays(library="np")
    return data


def load_data_3N(filepath):
    """Load 3N TTree into dict of numpy arrays."""
    f = uproot.open(filepath)
    tree = f["events"]
    return tree.arrays(library="np")


def add_derived_3N(d3):
    """Add derived kinematic quantities from raw 3N 4-vectors."""
    e_px, e_py, e_pz = d3["electron"][:, 0], d3["electron"][:, 1], d3["electron"][:, 2]
    e_E = d3["electron"][:, 3]

    l_px, l_py, l_pz = d3["lead"][:, 0], d3["lead"][:, 1], d3["lead"][:, 2]
    r2_px, r2_py, r2_pz = d3["recoil2"][:, 0], d3["recoil2"][:, 1], d3["recoil2"][:, 2]
    r3_px, r3_py, r3_pz = d3["recoil3"][:, 0], d3["recoil3"][:, 1], d3["recoil3"][:, 2]

    q_px, q_py, q_pz = d3["q"][:, 0], d3["q"][:, 1], d3["q"][:, 2]
    q_mag = np.sqrt(q_px**2 + q_py**2 + q_pz**2)

    # p1 = pmiss = lead - q (initial lead momentum)
    p1_x = l_px - q_px
    p1_y = l_py - q_py
    p1_z = l_pz - q_pz
    d3["p1_mag"] = np.sqrt(p1_x**2 + p1_y**2 + p1_z**2)
    d3["p1_after_mag"] = np.sqrt(l_px**2 + l_py**2 + l_pz**2)
    d3["p2_mag"] = np.sqrt(r2_px**2 + r2_py**2 + r2_pz**2)
    d3["p3_mag"] = np.sqrt(r3_px**2 + r3_py**2 + r3_pz**2)

    def angle_deg(ax, ay, az, bx, by, bz):
        dot = ax*bx + ay*by + az*bz
        a_mag = np.sqrt(ax**2 + ay**2 + az**2)
        b_mag = np.sqrt(bx**2 + by**2 + bz**2)
        return np.degrees(np.arccos(np.clip(dot / (a_mag * b_mag + 1e-30), -1, 1)))

    d3["theta12"] = angle_deg(p1_x, p1_y, p1_z, r2_px, r2_py, r2_pz)
    d3["theta13"] = angle_deg(p1_x, p1_y, p1_z, r3_px, r3_py, r3_pz)
    d3["theta23"] = angle_deg(r2_px, r2_py, r2_pz, r3_px, r3_py, r3_pz)

    d3["p1_angle_q"]       = angle_deg(p1_x, p1_y, p1_z, q_px, q_py, q_pz)
    d3["p1_after_angle_q"] = angle_deg(l_px, l_py, l_pz, q_px, q_py, q_pz)
    d3["p2_angle_q"]       = angle_deg(r2_px, r2_py, r2_pz, q_px, q_py, q_pz)
    d3["p3_angle_q"]       = angle_deg(r3_px, r3_py, r3_pz, q_px, q_py, q_pz)
    d3["e_angle_q"]        = angle_deg(e_px, e_py, e_pz, q_px, q_py, q_pz)

    d3["e_mom"]       = e_E
    d3["lead_over_q"] = d3["p1_after_mag"] / (q_mag + 1e-30)

    # Interplane angle: angle between planes {pmiss, p2} and {pmiss, p3}
    pmiss_3 = np.column_stack([p1_x, p1_y, p1_z])
    p2_3 = np.column_stack([r2_px, r2_py, r2_pz])
    p3_3 = np.column_stack([r3_px, r3_py, r3_pz])
    nn1 = np.cross(pmiss_3, p2_3)
    nn2 = np.cross(pmiss_3, p3_3)
    nn1_mag = np.linalg.norm(nn1, axis=1, keepdims=True)
    nn2_mag = np.linalg.norm(nn2, axis=1, keepdims=True)
    cos_ip = np.clip(np.abs(np.sum(nn1 / (nn1_mag + 1e-30) * nn2 / (nn2_mag + 1e-30), axis=1)), 0, 1)
    d3["interplane_angle"] = np.degrees(np.arccos(cos_ip))

    return d3


def add_derived(d):
    """Add derived kinematic quantities from raw 4-vectors."""
    # Shorthand for 4-vector components: [px, py, pz, E]
    e_px, e_py, e_pz = d["electron"][:, 0], d["electron"][:, 1], d["electron"][:, 2]
    e_E = d["electron"][:, 3]

    lp_px, lp_py, lp_pz = d["lead_post"][:, 0], d["lead_post"][:, 1], d["lead_post"][:, 2]
    rp_px, rp_py, rp_pz = d["recoil_post"][:, 0], d["recoil_post"][:, 1], d["recoil_post"][:, 2]

    q_px, q_py, q_pz = d["q"][:, 0], d["q"][:, 1], d["q"][:, 2]
    q_mag = np.sqrt(q_px**2 + q_py**2 + q_pz**2)

    lpre_px, lpre_py, lpre_pz = d["lead_pre"][:, 0], d["lead_pre"][:, 1], d["lead_pre"][:, 2]
    rpre_px, rpre_py, rpre_pz = d["recoil_pre"][:, 0], d["recoil_pre"][:, 1], d["recoil_pre"][:, 2]

    # p1 = pmiss = lead_post - q (reconstructed initial lead)
    p1_x = lp_px - q_px
    p1_y = lp_py - q_py
    p1_z = lp_pz - q_pz
    d["p1_mag"] = np.sqrt(p1_x**2 + p1_y**2 + p1_z**2)

    # p1_after = lead_post (final lead)
    d["p1_after_mag"] = np.sqrt(lp_px**2 + lp_py**2 + lp_pz**2)

    # p2 = recoil_post
    d["p2_mag"] = np.sqrt(rp_px**2 + rp_py**2 + rp_pz**2)

    # p3_hyp = -(p1 + p2) (hypothetical 3rd nucleon, zero-CM assumption)
    p3h_x = -(p1_x + rp_px)
    p3h_y = -(p1_y + rp_py)
    p3h_z = -(p1_z + rp_pz)
    d["p3_hyp_mag"] = np.sqrt(p3h_x**2 + p3h_y**2 + p3h_z**2)

    # Angles between vectors (returns degrees)
    def angle_deg(ax, ay, az, bx, by, bz):
        dot = ax*bx + ay*by + az*bz
        a_mag = np.sqrt(ax**2 + ay**2 + az**2)
        b_mag = np.sqrt(bx**2 + by**2 + bz**2)
        cos_angle = np.clip(dot / (a_mag * b_mag + 1e-30), -1, 1)
        return np.degrees(np.arccos(cos_angle))

    def theta_deg(px, py, pz):
        """Polar angle with z-axis (beam direction)."""
        p_mag = np.sqrt(px**2 + py**2 + pz**2)
        return np.degrees(np.arccos(np.clip(pz / (p_mag + 1e-30), -1, 1)))

    # theta12 = angle(p1, p2)
    d["theta12"] = angle_deg(p1_x, p1_y, p1_z, rp_px, rp_py, rp_pz)
    # theta23 = angle(p2, p3_hyp)
    d["theta23"] = angle_deg(rp_px, rp_py, rp_pz, p3h_x, p3h_y, p3h_z)

    # Angles with q
    d["p1_angle_q"] = angle_deg(p1_x, p1_y, p1_z, q_px, q_py, q_pz)
    d["p1_after_angle_q"] = angle_deg(lp_px, lp_py, lp_pz, q_px, q_py, q_pz)
    d["p2_angle_q"] = angle_deg(rp_px, rp_py, rp_pz, q_px, q_py, q_pz)
    d["p3_hyp_angle_q"] = angle_deg(p3h_x, p3h_y, p3h_z, q_px, q_py, q_pz)
    d["e_angle_q"] = angle_deg(e_px, e_py, e_pz, q_px, q_py, q_pz)

    # Angles with beam (z-axis)
    d["p1_angle_beam"] = theta_deg(p1_x, p1_y, p1_z)
    d["p1_after_angle_beam"] = theta_deg(lp_px, lp_py, lp_pz)
    d["p2_angle_beam"] = theta_deg(rp_px, rp_py, rp_pz)

    # Pre-FSI derived
    d["p1_pre_mag"] = np.sqrt((lpre_px - q_px)**2 + (lpre_py - q_py)**2 + (lpre_pz - q_pz)**2)
    d["p2_pre_mag"] = np.sqrt(rpre_px**2 + rpre_py**2 + rpre_pz**2)
    d["p1_pre_after_mag"] = np.sqrt(lpre_px**2 + lpre_py**2 + lpre_pz**2)

    # Pre-FSI angles with q
    d["p1_pre_angle_q"] = angle_deg(lpre_px - q_px, lpre_py - q_py, lpre_pz - q_pz, q_px, q_py, q_pz)
    d["p2_pre_angle_q"] = angle_deg(rpre_px, rpre_py, rpre_pz, q_px, q_py, q_pz)
    # Pre-FSI hypothetical 3rd: -(p1_pre + p2_pre)
    p3pre_x = -((lpre_px - q_px) + rpre_px)
    p3pre_y = -((lpre_py - q_py) + rpre_py)
    p3pre_z = -((lpre_pz - q_pz) + rpre_pz)
    d["p3_pre_angle_q"] = angle_deg(p3pre_x, p3pre_y, p3pre_z, q_px, q_py, q_pz)

    # Pair observables
    d["pair_cm_mom"] = np.sqrt((p1_x + rp_px)**2 + (p1_y + rp_py)**2 + (p1_z + rp_pz)**2)
    d["pair_rel_mom"] = np.sqrt((p1_x - rp_px)**2 + (p1_y - rp_py)**2 + (p1_z - rp_pz)**2) / 2.

    # p_lead / q ratio
    d["lead_over_q"] = d["p1_after_mag"] / (q_mag + 1e-30)

    # Electron momentum
    d["e_mom"] = e_E

    # Missing mass: mmiss = sqrt((2mN + omega - EN)^2 - |pmiss|^2)
    lp_E = d["lead_post"][:, 3]
    omega = d["nu"]
    mmiss2 = (2*mN + omega - lp_E)**2 - d["p1_mag"]**2
    d["mmiss"] = np.sqrt(np.maximum(mmiss2, 0.0))

    # Find highest-momentum FSI secondary nucleon above kF for each event
    # This requires reading the jagged arrays
    # For now, compute nAboveKF-based theta12/theta23 using the stored nAboveKF
    # The p3_fsi info needs the secondary arrays

    return d


def find_p3_fsi(d):
    """
    Find the highest-momentum FSI secondary nucleon above kF for N>=3 events.
    Returns p3_fsi as (px, py, pz) arrays, and a mask of events where p3_fsi exists.
    """
    n = len(d["weight"])
    p3_px = np.zeros(n)
    p3_py = np.zeros(n)
    p3_pz = np.zeros(n)
    p3_pdg = np.zeros(n, dtype=int)
    has_p3 = np.zeros(n, dtype=bool)

    nSec = d["nSec"]
    sec_pdg = d["sec_pdg"]
    sec_px_arr = d["sec_px"]
    sec_py_arr = d["sec_py"]
    sec_pz_arr = d["sec_pz"]

    for i in range(n):
        if nSec[i] == 0:
            continue
        best_mag = -1.0
        for j in range(nSec[i]):
            pdg = sec_pdg[i][j]
            if pdg == 2212 or pdg == 2112:
                px, py, pz = sec_px_arr[i][j], sec_py_arr[i][j], sec_pz_arr[i][j]
                mag = np.sqrt(px**2 + py**2 + pz**2)
                if mag > kF and mag > best_mag:
                    best_mag = mag
                    p3_px[i], p3_py[i], p3_pz[i] = px, py, pz
                    p3_pdg[i] = pdg
                    has_p3[i] = True

    d["p3_fsi_px"] = p3_px
    d["p3_fsi_py"] = p3_py
    d["p3_fsi_pz"] = p3_pz
    d["p3_fsi_mag"] = np.sqrt(p3_px**2 + p3_py**2 + p3_pz**2)
    d["p3_fsi_pdg"] = p3_pdg
    d["has_p3_fsi"] = has_p3

    # Compute N=3 theta12, theta23 using p3_fsi
    q_px, q_py, q_pz = d["q"][:, 0], d["q"][:, 1], d["q"][:, 2]
    lp_px, lp_py, lp_pz = d["lead_post"][:, 0], d["lead_post"][:, 1], d["lead_post"][:, 2]
    rp_px, rp_py, rp_pz = d["recoil_post"][:, 0], d["recoil_post"][:, 1], d["recoil_post"][:, 2]
    p1_x, p1_y, p1_z = lp_px - q_px, lp_py - q_py, lp_pz - q_pz

    def angle_deg(ax, ay, az, bx, by, bz):
        dot = ax*bx + ay*by + az*bz
        a_mag = np.sqrt(ax**2 + ay**2 + az**2)
        b_mag = np.sqrt(bx**2 + by**2 + bz**2)
        return np.degrees(np.arccos(np.clip(dot / (a_mag * b_mag + 1e-30), -1, 1)))

    d["theta12_n3"] = angle_deg(p1_x, p1_y, p1_z, rp_px, rp_py, rp_pz)
    d["theta23_n3"] = angle_deg(rp_px, rp_py, rp_pz, p3_px, p3_py, p3_pz)
    d["p3_fsi_angle_q"] = angle_deg(p3_px, p3_py, p3_pz, q_px, q_py, q_pz)
    d["pmiss_angle_p3_fsi"] = angle_deg(p1_x, p1_y, p1_z, p3_px, p3_py, p3_pz)

    # Interplane angle for N=3: angle between planes {pmiss, recoil} and {pmiss, p3_fsi}
    pmiss_3 = np.column_stack([p1_x, p1_y, p1_z])
    p2_3 = np.column_stack([rp_px, rp_py, rp_pz])
    p3_3 = np.column_stack([p3_px, p3_py, p3_pz])
    nn1 = np.cross(pmiss_3, p2_3)
    nn2 = np.cross(pmiss_3, p3_3)
    nn1_mag = np.linalg.norm(nn1, axis=1, keepdims=True)
    nn2_mag = np.linalg.norm(nn2, axis=1, keepdims=True)
    cos_ip = np.clip(np.abs(np.sum(nn1 / (nn1_mag + 1e-30) * nn2 / (nn2_mag + 1e-30), axis=1)), 0, 1)
    d["interplane_angle_n3"] = np.degrees(np.arccos(cos_ip))

    return d


# ──────────── Plotting helpers ────────────
def plot_1d(ax, values, weights, bins, range_, **kwargs):
    """Weighted 1D histogram as step plot."""
    counts, edges = np.histogram(values, bins=bins, range=range_, weights=weights)
    centers = 0.5 * (edges[:-1] + edges[1:])
    # Normalize to unit area
    total = np.sum(counts)
    if total > 0:
        counts = counts / total
    ax.step(centers, counts, where='mid', **kwargs)
    return centers, counts

Q2min = 1.0


def build_cuts(cut_list):
    """Combine a list of (label, mask) pairs into (combined_mask, description_string)."""
    mask = cut_list[0][1]
    for _, m in cut_list[1:]:
        mask = mask & m
    desc = ",  ".join(lbl for lbl, _ in cut_list)
    return mask, desc


# ──────────── Main analysis functions ────────────
def plot_theta_heatmaps(d, out_dir):
    """2D theta12-theta23 heatmaps: all events, N=2 only, N=3 only."""
    bins = 400

    # Define cut parameters (change values here → labels auto-update)
    angle_e_max = 45
    angle_pq_max = 8
    lead_q_min = 0.75
    pmiss_lo, pmiss_hi = 0.25, 0.9
    xB_max = 1.2
    interplane_max = 20

    shared_cuts = [
        (rf"$\theta_{{pq}}<{angle_pq_max}^\circ$", d["p1_after_angle_q"] < angle_pq_max),
        (rf"$p_N/q>{lead_q_min}$", d["lead_over_q"] > lead_q_min),
        (rf"${pmiss_lo}<p_{{\rm miss}}<{pmiss_hi}$", (pmiss_lo < d["pmiss"]) & (d["pmiss"] < pmiss_hi)),
        (rf"$x_B<{xB_max}$", d["xB"] < xB_max),
    ]

    cuts_N2_list = [(r"$N_{kF}=2$", d["nAboveKF"] == 2)] + shared_cuts
    cuts_N2, desc_N2 = build_cuts(cuts_N2_list)

    cuts_N3_list = [(r"$N_{kF}\geq3$", (d["nAboveKF"] >= 3) & d["has_p3_fsi"])] + shared_cuts + [
        (rf"$\phi_{{\rm interplane}}<{interplane_max}^\circ$", d["interplane_angle_n3"] < interplane_max),
    ]
    cuts_N3, desc_N3 = build_cuts(cuts_N3_list)

    # Triangle geometry for region boundaries
    A, B = region_params()

    def _line_intersection(m1, b1, m2, b2):
        denom = (m1 - m2)
        if abs(denom) < 1e-15:
            return None
        x = (b2 - b1) / denom
        y = m1 * x + b1
        return float(x), float(y)

    def _add_triangle(ax, points, *, edgecolor, label):
        pts = [p for p in points if p is not None and np.isfinite(p[0]) and np.isfinite(p[1])]
        if len(pts) != 3:
            return
        cx = sum(p[0] for p in pts) / 3.0
        cy = sum(p[1] for p in pts) / 3.0
        pts = sorted(pts, key=lambda p: np.arctan2(p[1] - cy, p[0] - cx))
        poly = Polygon(pts, closed=True, fill=False, edgecolor=edgecolor, linewidth=2.0, label=label)
        ax.add_patch(poly)

    m1_R, b1_R = B, 180.0 * (1.0 - B)
    m2_R, b2_R = (1.0 / B), 180.0 * (1.0 - (1.0 / B))
    m3_R, b3_R = -1.0, (A + 180.0 - (180.0 - A) * B)
    p12_R = _line_intersection(m1_R, b1_R, m2_R, b2_R)
    p13_R = _line_intersection(m1_R, b1_R, m3_R, b3_R)
    p23_R = _line_intersection(m2_R, b2_R, m3_R, b3_R)

    s1_L = -(m1_R / (1.0 + m1_R))
    s2_L = -(m2_R / (1.0 + m2_R))
    b_L = 180.0
    x_vert_L = 360.0 - b3_R
    p12_L = (0.0, 180.0)
    p13_L = (x_vert_L, s1_L * x_vert_L + b_L)
    p23_L = (x_vert_L, s2_L * x_vert_L + b_L)

    # Region BR vertices: transform R vertices (x, y) -> (x, 360-x-y)
    p12_BR = (p12_R[0], 360.0 - p12_R[0] - p12_R[1]) if p12_R else None
    p13_BR = (p13_R[0], 360.0 - p13_R[0] - p13_R[1]) if p13_R else None
    p23_BR = (p23_R[0], 360.0 - p23_R[0] - p23_R[1]) if p23_R else None

    # Center circle: radius chosen so area = triangle area
    tri_area = _triangle_area(p12_R, p13_R, p23_R) if (p12_R and p13_R and p23_R) else 0.0
    center_radius = np.sqrt(tri_area / np.pi) if tri_area > 0 else 10.0

    base_desc = rf"$\theta_e<{angle_e_max}^\circ$,  $Q^2\geq{Q2min}$"

    for label, mask_fn, fname_base, cut_desc in [
        ("All events", lambda: np.ones(len(d["weight"]), dtype=bool), "theta_heatmap_all", base_desc),
        ("N=2 events", lambda: d["nAboveKF"] == 2, "theta_heatmap_N2", rf"$N_{{kF}}=2$,  {base_desc}"),
        ("N=3 events", lambda: (d["nAboveKF"] >= 3) & d["has_p3_fsi"], "theta_heatmap_N3", rf"$N_{{kF}}\geq3$,  {base_desc}"),
        ("N=2 events with cuts", lambda: cuts_N2, "theta_heatmap_N2_cuts", desc_N2 + rf",  {base_desc}"),
        ("N=3 events with cuts", lambda: cuts_N3, "theta_heatmap_N3_cuts", desc_N3 + rf",  {base_desc}"),
    ]:
        mask = mask_fn() & (d["scattering_angle"] < angle_e_max) & (d["Q2"] >= Q2min)
        if "N=3" in label:
            t12_orig = d["theta12_n3"][mask]
            t12_swap = d["pmiss_angle_p3_fsi"][mask]
            t23_orig = d["theta23_n3"][mask]
            t12 = np.concatenate([t12_orig, t12_swap])
            t23 = np.concatenate([t23_orig, t23_orig])
            w = np.concatenate([d["weight"][mask], d["weight"][mask]])
        else:
            t12 = d["theta12"][mask]
            t23 = d["theta23"][mask]
            w = d["weight"][mask]

        w_tot = sum(w) + 1e-30
        fL = sum(w[in_region_L(t12, t23)]) / w_tot
        fR = sum(w[in_region_R(t12, t23)]) / w_tot
        fBR = sum(w[in_region_BR(t12, t23)]) / w_tot
        fC = sum(w[in_region_center(t12, t23, center_radius)]) / w_tot
        print("Current selection: " + label)
        print(f"  L: {fL:.6f}  R: {fR:.6f}  BR: {fBR:.6f}  Center(r={center_radius:.2f}): {fC:.6f}")

        h, xe, ye = np.histogram2d(t12, t23, bins=bins, range=[[0, 180], [0, 180]], weights=w)
        h_norm = h / (np.sum(h) + 1e-30)
        h_norm[h_norm == 0] = np.nan

        for show_tri, suffix in [(False, ""), (True, "_triangles")]:
            fig, ax = plt.subplots(figsize=(7, 6))

            x_line = np.linspace(0, 180, 100)
            ax.plot(x_line, 180 - x_line / 2, 'r--')
            ax.plot(180 - x_line / 2, x_line, 'r--')
            ax.plot(x_line, x_line, 'r--')

            im = ax.pcolormesh(xe, ye, h_norm.T, norm=LogNorm() if np.any(h_norm > 0) else None)
            fig.colorbar(im, ax=ax, label="Normalized weight")
            ax.set_xlabel(r"$\theta_{12}$ (deg)")
            ax.set_ylabel(r"$\theta_{23}$ (deg)")
            ax.set_title(f"2N+FSI: {label}")

            # Center region circle (always drawn)
            circle = plt.Circle((120, 120), center_radius, fill=False,
                                edgecolor='green', linewidth=2.0, label='Region C')
            ax.add_patch(circle)

            if show_tri:
                _add_triangle(ax, [p12_R, p13_R, p23_R], edgecolor='b', label='Region R')
                _add_triangle(ax, [p12_L, p13_L, p23_L], edgecolor='purple', label='Region L')
                _add_triangle(ax, [p12_BR, p13_BR, p23_BR], edgecolor='orange', label='Region BR')
                ax.legend()
            else:
                ax.legend()

            ax.text(0.02, 0.98,
                    f"L: {fL:.4f}\nR: {fR:.4f}\nBR: {fBR:.4f}\nC: {fC:.4f}",
                    transform=ax.transAxes, va='top', ha='left', fontsize=8,
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

            fig.text(0.5, -0.02, f"Cuts: {cut_desc}", ha='center', va='top', fontsize=7,
                     wrap=True, transform=fig.transFigure)
            fig.tight_layout(rect=[0, 0.04, 1, 1])
            fname = f"{fname_base}{suffix}.png"
            fig.savefig(os.path.join(out_dir, fname), dpi=200, bbox_inches='tight')
            plt.close(fig)
            print(f"  {fname}")


def plot_theta_ratio_heatmap(d, d3, out_dir):
    """Ratio heatmap: 3N / (2N+FSI N=3) in theta12-theta23 space.
    Both are symmetrized under 2<->3 swap."""
    bins = 50

    cuts_N3 = ((d["nAboveKF"] >= 3) & d["has_p3_fsi"]
               & (d["p1_after_angle_q"] < 8) & (d["lead_over_q"] > 0.75)
               & (0.25 < d["pmiss"]) & (d["pmiss"] < 0.9)
               & (d["xB"] < 1.2) & (d["interplane_angle_n3"] < 20))

    for pass_label, mask_fsi, mask_3N, fname in [
        ("angle only",
         (d["nAboveKF"] >= 3) & d["has_p3_fsi"] & (d["scattering_angle"] < 45) & (d["Q2"] >= Q2min),
         (d3["p1_mag"] > kF) & (d3["p2_mag"] > kF) & (d3["p3_mag"] > kF) & (d3["scattering_angle"] < 45) & (d3["Q2"] >= Q2min),
         "theta_heatmap_ratio_3N_over_N3_angle_only.png"),
        ("full cuts",
         cuts_N3 & (d["scattering_angle"] < 45) & (d["Q2"] >= Q2min),
         ((d3["p1_mag"] > kF) & (d3["p2_mag"] > kF) & (d3["p3_mag"] > kF)
          & (d3["scattering_angle"] < 45) & (d3["Q2"] >= Q2min)
          & (d3["p1_after_angle_q"] < 8) & (d3["lead_over_q"] > 0.75)
          & (0.25 < d3["p1_mag"]) & (d3["p1_mag"] < 0.9)
          & (d3["xB"] < 1.2) & (d3["interplane_angle"] < 20)),
         "theta_heatmap_ratio_3N_over_N3_full_cuts.png"),
    ]:
        # 2N+FSI N=3 (symmetrized)
        t12_fsi = np.concatenate([d["theta12_n3"][mask_fsi], d["pmiss_angle_p3_fsi"][mask_fsi]])
        t23_fsi = np.concatenate([d["theta23_n3"][mask_fsi], d["theta23_n3"][mask_fsi]])
        w_fsi   = np.concatenate([d["weight"][mask_fsi], d["weight"][mask_fsi]])

        # 3N SRC (symmetrized)
        t12_3N = np.concatenate([d3["theta12"][mask_3N], d3["theta13"][mask_3N]])
        t23_3N = np.concatenate([d3["theta23"][mask_3N], d3["theta23"][mask_3N]])
        w_3N   = np.concatenate([d3["weight"][mask_3N], d3["weight"][mask_3N]])

        h_fsi, xe, ye = np.histogram2d(t12_fsi, t23_fsi, bins=bins,
                                        range=[[0, 180], [0, 180]], weights=w_fsi)
        h_3N, _, _    = np.histogram2d(t12_3N, t23_3N, bins=bins,
                                        range=[[0, 180], [0, 180]], weights=w_3N)

        # Normalize each to unit area
        h_fsi_norm = h_fsi / (np.sum(h_fsi) + 1e-30)
        h_3N_norm  = h_3N  / (np.sum(h_3N) + 1e-30)

        # Ratio: 3N / N=3
        ratio = np.full_like(h_3N_norm, np.nan)
        valid = (h_fsi_norm > 0) & (h_3N_norm > 0)
        ratio[valid] = h_3N_norm[valid] / h_fsi_norm[valid]

        fig, ax = plt.subplots(figsize=(7, 6))
        x_line = np.linspace(0, 180, 100)
        ax.plot(x_line, 180 - x_line / 2, 'r--')
        ax.plot(180 - x_line / 2, x_line, 'r--')
        ax.plot(x_line, x_line, 'r--')

        im = ax.pcolormesh(xe, ye, ratio.T, norm=LogNorm() if np.any(ratio > 0) else None)
        fig.colorbar(im, ax=ax, label="3N SRC / 2N+FSI N=3")
        ax.set_xlabel(r"$\theta_{12}$ (deg)")
        ax.set_ylabel(r"$\theta_{23}$ (deg)")
        ax.set_title(f"3N / N=3 ratio ({pass_label}, symmetrized)")
        fig.tight_layout()
        fig.savefig(os.path.join(out_dir, fname), dpi=200)
        plt.close(fig)
        print(f"  {fname}")


def plot_1d_distributions(d, out_dir, Q2_min=2.0, Q2_max=10.0):
    """Standard 1D kinematic distributions."""
    mask = (d["Q2"] >= Q2_min) & (d["Q2"] < Q2_max)
    w = d["weight"][mask]

    variables = [
        ("xB", d["xB"][mask], (0, 4), r"$x_B$"),
        ("Q2", d["Q2"][mask], (2, 10), r"$Q^2$ [GeV$^2$]"),
        ("scattering_angle", d["scattering_angle"][mask], (0, 50), "Electron angle with beam [deg]"),
        ("e_angle_q", d["e_angle_q"][mask], (30, 90), r"Electron angle with $\vec{q}$ [deg]"),
        ("e_mom", d["e_mom"][mask], (0, 8), "Electron momentum [GeV/c]"),
        ("pmiss", d["pmiss"][mask], (0, 1.5), r"$p_{\mathrm{miss}}$ [GeV/c]"),
        ("p1_angle_beam", d["p1_angle_beam"][mask], (0, 180), "Incoming lead angle with beam [deg]"),
        ("p1_angle_q", d["p1_angle_q"][mask], (0, 180), r"Incoming lead angle with $\vec{q}$ [deg]"),
        ("p2_mag", d["p2_mag"][mask], (0, 1.5), "Recoil momentum [GeV/c]"),
        ("p2_angle_beam", d["p2_angle_beam"][mask], (0, 180), "Recoil angle with beam [deg]"),
        ("p2_angle_q", d["p2_angle_q"][mask], (0, 180), r"Recoil angle with $\vec{q}$ [deg]"),
        ("p1_after_angle_beam", d["p1_after_angle_beam"][mask], (0, 180), "Outgoing lead angle with beam [deg]"),
        ("p1_after_angle_q", d["p1_after_angle_q"][mask], (0, 50), r"Outgoing lead angle with $\vec{q}$ [deg]"),
        ("p1_after_mag", d["p1_after_mag"][mask], (0, 8), "Outgoing lead momentum [GeV/c]"),
        ("theta12", d["theta12"][mask], (0, 180), r"$\theta_{12}$ [deg]"),
        ("pair_cm_mom", d["pair_cm_mom"][mask], (0, 1.0), r"$|p_1+p_2|$ [GeV/c]"),
        ("pair_rel_mom", d["pair_rel_mom"][mask], (0, 1.0), r"$|p_1-p_2|/2$ [GeV/c]"),
        ("p3_hyp_mag", d["p3_hyp_mag"][mask], (0, 1.5), r"$|p_3|$ (hyp.) [GeV/c]"),
        ("p3_hyp_angle_q", d["p3_hyp_angle_q"][mask], (0, 180), r"$p_3$ (hyp.) angle with $\vec{q}$ [deg]"),
    ]

    for name, vals, range_, xlabel in variables:
        fig, ax = plt.subplots(figsize=(7, 5))
        plot_1d(ax, vals, w, 45, range_, label="2N+FSI", linewidth=1.5, color='tab:blue')
        ax.set_xlabel(xlabel, fontsize=12)
        ax.set_ylabel("Normalized weight", fontsize=12)
        ax.set_title(f"2N+FSI: {xlabel}")
        ax.legend()
        fig.tight_layout()
        fig.savefig(os.path.join(out_dir, f"hist_{name}_1D.png"), dpi=200)
        plt.close(fig)

    print(f"  {len(variables)} 1D distribution plots saved")


def plot_xB_Q2_heatmap(d, out_dir, Q2_min=2.0, Q2_max=10.0):
    """2D xB-Q2 heatmap."""
    mask = (d["Q2"] >= Q2_min) & (d["Q2"] < Q2_max) & (d["xB"] >= 0) & (d["xB"] < 4)
    fig, ax = plt.subplots(figsize=(7, 6))
    h, xe, ye = np.histogram2d(d["xB"][mask], d["Q2"][mask], bins=[100, 100],
                                range=[[0, 4], [Q2_min, Q2_max]], weights=d["weight"][mask])
    h[h == 0] = np.nan
    im = ax.pcolormesh(xe, ye, h.T, norm=LogNorm())
    fig.colorbar(im, ax=ax, label="Weight")
    ax.set_xlabel(r"$x_B$")
    ax.set_ylabel(r"$Q^2$ [GeV$^2$]")
    ax.set_title("2N+FSI: xB vs Q2")
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "xB_Q2_heatmap.png"), dpi=200)
    plt.close(fig)
    print("  xB_Q2_heatmap.png")


def plot_region_overlays(d, d3, out_dir):
    """Compare N=2 vs 3N and N=3 vs 3N kinematic distributions.

    Produces plots for:
      - 5 groups: (angle_only × {all, region_L, region_R}) + (full_cuts × {region_L, region_R})
      - 2 comparisons per group: N=2 vs 3N, N=3 vs 3N
      - 12 kinematic variables each
    Total: 120 plots, organized in subdirectories.
    """
    # ── 2N+FSI masks ──
    mask_N2_base = (d["nAboveKF"] == 2) & (d["scattering_angle"] < 45) & (d["Q2"] >= Q2min)
    mask_N2_full = mask_N2_base & (d["p1_after_angle_q"] < 8) & (d["lead_over_q"] > 0.75) & (0.25 < d["pmiss"]) & (d["pmiss"] < 0.9) & (d["xB"] < 1.2)

    mask_N3_base = (d["nAboveKF"] >= 3) & d["has_p3_fsi"] & (d["scattering_angle"] < 45) & (d["Q2"] >= Q2min)
    mask_N3_full = mask_N3_base & (d["p1_after_angle_q"] < 8) & (d["lead_over_q"] > 0.75) & (0.25 < d["pmiss"]) & (d["pmiss"] < 0.9) & (d["xB"] < 1.2) & (d["interplane_angle_n3"] < 20)

    # ── 3N masks ──
    base_3N = (d3["p1_mag"] > kF) & (d3["p2_mag"] > kF) & (d3["p3_mag"] > kF) & (d3["scattering_angle"] < 45) & (d3["Q2"] >= Q2min)
    mask_3N_N2full = base_3N & (d3["p1_after_angle_q"] < 8) & (d3["lead_over_q"] > 0.75) & (0.25 < d3["p1_mag"]) & (d3["p1_mag"] < 0.9) & (d3["xB"] < 1.2)
    mask_3N_N3full = mask_3N_N2full & (d3["interplane_angle"] < 20)

    # ── Region boolean arrays ──
    rl_N2   = in_region_L(d["theta12"],    d["theta23"])
    rr_N2   = in_region_R(d["theta12"],    d["theta23"])
    rbr_N2  = in_region_BR(d["theta12"],   d["theta23"])
    rl_N3   = in_region_L(d["theta12_n3"], d["theta23_n3"])
    rr_N3   = in_region_R(d["theta12_n3"], d["theta23_n3"])
    rbr_N3  = in_region_BR(d["theta12_n3"], d["theta23_n3"])
    rl_3N   = in_region_L(d3["theta12"],   d3["theta23"])
    rr_3N   = in_region_R(d3["theta12"],   d3["theta23"])
    rbr_3N  = in_region_BR(d3["theta12"],  d3["theta23"])

    # ── Variable list ──
    variables = [
        ("init_lead_angle_q",    r"Initial lead angle with $\vec{q}$",    "deg",        (0, 180)),
        ("init_rec1_angle_q",    r"Initial recoil 1 angle with $\vec{q}$","deg",        (0, 180)),
        ("pmiss_LR",             r"$p_{\mathrm{miss}}$",                   "GeV/c",      (0, 1.5)),
        ("p1_after_angle_q_LR",  r"Final lead angle with $\vec{q}$",      "deg",        (0, 50)),
        ("p2_angle_q_LR",        r"Final recoil 1 angle with $\vec{q}$",  "deg",        (0, 180)),
        ("p2_final_mom",         r"$|p_2|$ (final)",                       "GeV/c",      (0, 1.5)),
        ("e_angle_q_LR",         r"Electron angle with $\vec{q}$",         "deg",        (0, 90)),
        ("e_mom_LR",             r"Electron momentum",                     "GeV/c",      (0, 8)),
        ("xB_LR",                r"$x_B$",                                 "",           (0, 4)),
        ("Q2_LR",                r"$Q^2$",                                 r"GeV$^2$",   (2, 10)),
        ("lead_mom_over_q",      r"$|p_{\mathrm{lead}}|/|\vec{q}|$",      "",           (0, 2)),
        ("p1_final_mom",         r"$|p_1|$ (final)",                       "GeV/c",      (0, 8)),
        ("p3_final_mom",         r"$|p_3|$ (final)",                       "GeV/c",      (0, 1.5)),
        ("p3_angle_q",           r"$p_3$ angle with $\vec{q}$",            "deg",        (0, 180)),
        ("pmiss_angle_p2",       r"$\angle(p_{\mathrm{miss}}, p_2)$",      "deg",        (0, 180)),
        ("pmiss_angle_p3",       r"$\angle(p_{\mathrm{miss}}, p_3)$",      "deg",        (0, 180)),
    ]

    # ── Value accessors ──
    N2_map = {
        "init_lead_angle_q":   d["p1_pre_angle_q"],
        "init_rec1_angle_q":   d["p2_pre_angle_q"],
        "pmiss_LR":            d["pmiss"],
        "p1_after_angle_q_LR": d["p1_after_angle_q"],
        "p2_angle_q_LR":       d["p2_angle_q"],
        "p2_final_mom":        d["p2_mag"],
        "e_angle_q_LR":        d["e_angle_q"],
        "e_mom_LR":            d["e_mom"],
        "xB_LR":               d["xB"],
        "Q2_LR":               d["Q2"],
        "lead_mom_over_q":     d["lead_over_q"],
        "p1_final_mom":        d["p1_after_mag"],
        "pmiss_angle_p2":      d["theta12"],
    }
    N3_map = {
        "init_lead_angle_q":   d["p1_pre_angle_q"],
        "init_rec1_angle_q":   d["p2_pre_angle_q"],
        "pmiss_LR":            d["pmiss"],
        "p1_after_angle_q_LR": d["p1_after_angle_q"],
        "p2_angle_q_LR":       d["p2_angle_q"],
        "p2_final_mom":        d["p2_mag"],
        "e_angle_q_LR":        d["e_angle_q"],
        "e_mom_LR":            d["e_mom"],
        "xB_LR":               d["xB"],
        "Q2_LR":               d["Q2"],
        "lead_mom_over_q":     d["lead_over_q"],
        "p1_final_mom":        d["p1_after_mag"],
        "p3_final_mom":        d["p3_fsi_mag"],
        "p3_angle_q":          d["p3_fsi_angle_q"],
        "pmiss_angle_p2":      d["theta12"],
        "pmiss_angle_p3":      d["pmiss_angle_p3_fsi"],
    }
    SRC3_map = {
        "init_lead_angle_q":   d3["p1_angle_q"],
        "init_rec1_angle_q":   d3["p2_angle_q"],
        "pmiss_LR":            d3["p1_mag"],
        "p1_after_angle_q_LR": d3["p1_after_angle_q"],
        "p2_angle_q_LR":       d3["p2_angle_q"],
        "p2_final_mom":        d3["p2_mag"],
        "e_angle_q_LR":        d3["e_angle_q"],
        "e_mom_LR":            d3["e_mom"],
        "xB_LR":               d3["xB"],
        "Q2_LR":               d3["Q2"],
        "lead_mom_over_q":     d3["lead_over_q"],
        "p1_final_mom":        d3["p1_after_mag"],
        "p3_final_mom":        d3["p3_mag"],
        "p3_angle_q":          d3["p3_angle_q"],
        "pmiss_angle_p2":      d3["theta12"],
        "pmiss_angle_p3":      d3["theta13"],
    }

    # ── Groups: (cuts_label, region_label, mask_N2, mask_N3, mask_3N_for_N2, mask_3N_for_N3,
    #             reg_mask_N2, reg_mask_N3, reg_mask_3N) ──
    # For "all" no region mask is applied; for L/R the boolean arrays above are AND-ed in.
    groups = [
        # pass,         region_dir,  mN2_base,     mN3_base,     m3N_base, m3N_base,
        #   reg_N2,  reg_N3,  reg_3N
        ("angle_only", "all",
         mask_N2_base, mask_N3_base, base_3N,        base_3N,
         None, None, None),
        ("angle_only", "region_L",
         mask_N2_base, mask_N3_base, base_3N,        base_3N,
         rl_N2, rl_N3, rl_3N),
        ("angle_only", "region_R",
         mask_N2_base, mask_N3_base, base_3N,        base_3N,
         rr_N2, rr_N3, rr_3N),
        ("angle_only", "region_BR",
         mask_N2_base, mask_N3_base, base_3N,        base_3N,
         rbr_N2, rbr_N3, rbr_3N),
        ("full_cuts", "region_L",
         mask_N2_full, mask_N3_full, mask_3N_N2full, mask_3N_N3full,
         rl_N2, rl_N3, rl_3N),
        ("full_cuts", "region_R",
         mask_N2_full, mask_N3_full, mask_3N_N2full, mask_3N_N3full,
         rr_N2, rr_N3, rr_3N),
        ("full_cuts", "region_BR",
         mask_N2_full, mask_N3_full, mask_3N_N2full, mask_3N_N3full,
         rbr_N2, rbr_N3, rbr_3N),
    ]

    total = 0
    for (pass_name, region_dir,
         m_N2, m_N3, m_3N_for_N2, m_3N_for_N3,
         reg_N2, reg_N3, reg_3N) in groups:

        # Apply region mask where applicable
        mN2  = m_N2  & reg_N2  if reg_N2  is not None else m_N2
        mN3  = m_N3  & reg_N3  if reg_N3  is not None else m_N3
        m3_N2 = m_3N_for_N2 & reg_3N if reg_3N is not None else m_3N_for_N2
        m3_N3 = m_3N_for_N3 & reg_3N if reg_3N is not None else m_3N_for_N3

        for cmp_label, cmp_dir, mask_a, w_a, map_a, color_a, lbl_a, mask_b, w_b, map_b, color_b, lbl_b in [
            ("N2_vs_3N", "N2_vs_3N",
             mN2,  d["weight"],  N2_map,    'tab:blue',  "2N+FSI N=2",
             m3_N2, d3["weight"], SRC3_map,  'tab:red',   "3N SRC"),
            ("N3_vs_3N", "N3_vs_3N",
             mN3,  d["weight"],  N3_map,    'tab:green', "2N+FSI N=3",
             m3_N3, d3["weight"], SRC3_map,  'tab:red',   "3N SRC"),
        ]:
            folder = os.path.join(out_dir, "region_overlays", pass_name, region_dir, cmp_dir)
            os.makedirs(folder, exist_ok=True)

            for name, label, unit, range_ in variables:
                # Skip variables not available on both sides of the comparison
                if name not in map_a or name not in map_b:
                    continue
                fig, ax = plt.subplots(figsize=(7, 5))
                if name in map_a and np.sum(mask_a) > 0:
                    plot_1d(ax, map_a[name][mask_a], w_a[mask_a], 45, range_,
                            label=lbl_a, linewidth=1.5, color=color_a)
                if name in map_b and np.sum(mask_b) > 0:
                    plot_1d(ax, map_b[name][mask_b], w_b[mask_b], 45, range_,
                            label=lbl_b, linewidth=1.5, color=color_b)
                xlabel = f"{label}" + (f" [{unit}]" if unit else "")
                region_title = region_dir.replace("_", " ").title()
                pass_title   = pass_name.replace("_", " ")
                ax.set_xlabel(xlabel, fontsize=12)
                ax.set_ylabel("Normalized weight", fontsize=12)
                ax.set_title(f"{label} — {region_title} ({pass_title})")
                ax.legend()
                fig.tight_layout()
                fname = f"{name}_{cmp_label}_{region_dir}_{pass_name}.png"
                fig.savefig(os.path.join(folder, fname), dpi=200)
                plt.close(fig)
                total += 1

    print(f"  {total} region overlay plots saved")


def compute_region_fractions(d):
    """Print region weight fraction table."""
    # N=2
    mask_N2 = d["nAboveKF"] == 2
    mask_N2_cuts = mask_N2 & (d["scattering_angle"] < 45) & (d["p1_after_angle_q"] < 10) & (d["lead_over_q"] > 0.7)
    w_N2_total = np.sum(d["weight"][mask_N2_cuts])
    t12_n2 = d["theta12"][mask_N2_cuts]
    t23_n2 = d["theta23"][mask_N2_cuts]
    w_N2_L  = np.sum(d["weight"][mask_N2_cuts][in_region_L(t12_n2, t23_n2)])
    w_N2_R  = np.sum(d["weight"][mask_N2_cuts][in_region_R(t12_n2, t23_n2)])
    w_N2_BR = np.sum(d["weight"][mask_N2_cuts][in_region_BR(t12_n2, t23_n2)])

    # N=3
    mask_N3 = (d["nAboveKF"] >= 3) & d["has_p3_fsi"]
    mask_N3_cuts = mask_N3 & (d["scattering_angle"] < 45) & (d["p1_after_angle_q"] < 10) & (d["lead_over_q"] > 0.7) & (d["p2_mag"] > 0.5)
    w_N3_total = np.sum(d["weight"][mask_N3_cuts])
    t12_n3 = d["theta12_n3"][mask_N3_cuts]
    t23_n3 = d["theta23_n3"][mask_N3_cuts]
    w_N3_L  = np.sum(d["weight"][mask_N3_cuts][in_region_L(t12_n3, t23_n3)])
    w_N3_R  = np.sum(d["weight"][mask_N3_cuts][in_region_R(t12_n3, t23_n3)])
    w_N3_BR = np.sum(d["weight"][mask_N3_cuts][in_region_BR(t12_n3, t23_n3)])

    print("\n--- Region fractions with kinematic cuts ---")
    print("  (e angle < 45 deg, lead angle with q < 10 deg, pLead/q > 0.7, |p2|>0.5 for N=3)")
    print(f"{'Source':<20s} | {'Region L frac':>14s} | {'Region R frac':>14s} | {'Region BR frac':>15s}")
    print("-" * 73)
    if w_N2_total > 0:
        print(f"{'2N+FSI  N=2':<20s} | {100*w_N2_L/w_N2_total:>13.4f}% | {100*w_N2_R/w_N2_total:>13.4f}% | {100*w_N2_BR/w_N2_total:>14.4f}%")
    if w_N3_total > 0:
        print(f"{'2N+FSI  N=3':<20s} | {100*w_N3_L/w_N3_total:>13.4f}% | {100*w_N3_R/w_N3_total:>13.4f}% | {100*w_N3_BR/w_N3_total:>14.4f}%")


def plot_pair_cm_momentum(d, out_dir):
    """Pre-FSI pair CM momentum |p1_initial + p2_initial| distribution."""
    lp = d["lead_pre"][:, :3]
    rp = d["recoil_pre"][:, :3]
    q = d["q"][:, :3]
    pmiss_pre = lp - q
    cm = np.sqrt(np.sum((pmiss_pre + rp)**2, axis=1))
    w = d["weight"]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    bns, rng = 80, (0, 0.5)

    # Unweighted
    c1, e1 = np.histogram(cm, bins=bns, range=rng)
    centers = 0.5 * (e1[:-1] + e1[1:])
    c1 = c1 / (np.sum(c1) * (e1[1] - e1[0]))
    ax1.step(centers, c1, where='mid', linewidth=1.5, color='tab:green')
    ax1.set_title('Unweighted (raw event counts)')
    ax1.set_xlabel(r"$|\vec{p}_{1} + \vec{p}_{2}|$ (pre-FSI) [GeV/c]")
    ax1.set_ylabel('Normalized density')

    # Weighted
    c2, _ = np.histogram(cm, bins=bns, range=rng, weights=w)
    c2 = c2 / (np.sum(c2) * (e1[1] - e1[0]))
    ax2.step(centers, c2, where='mid', linewidth=1.5, color='tab:blue')
    ax2.set_title('Weighted (by event weight)')
    ax2.set_xlabel(r"$|\vec{p}_{1} + \vec{p}_{2}|$ (pre-FSI) [GeV/c]")
    ax2.set_ylabel('Normalized density')

    # Per-component std for annotation
    vcm = pmiss_pre + rp
    sig_x = np.std(vcm[:, 0])
    fig.suptitle(rf"Pre-FSI pair CM momentum ($\sigma_{{CM}} \approx {sig_x:.3f}$ GeV/c per component)", fontsize=13)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "pair_cm_pre_fsi.png"), dpi=200)
    plt.close(fig)
    print(f"  pair_cm_pre_fsi.png  (sigma_x={sig_x:.4f})")


def plot_wavefunction_vs_prel(d, out_dir):
    """Mean wavefunction S(p_rel) vs true initial relative momentum (pre-FSI, pre-photon)."""
    if "rho" not in d:
        print("  rho branch not found, skipping wavefunction vs p_rel plot")
        return
    rho = d["rho"]
    lp_pre = d["lead_pre"][:, :3]
    rp_pre = d["recoil_pre"][:, :3]
    q = d["q"][:, :3]
    p_lead_init = lp_pre - q
    p_rel = 0.5 * np.linalg.norm(p_lead_init - rp_pre, axis=1)

    mask = (rho > 0) & np.isfinite(rho)
    p_vals = p_rel[mask]
    r_vals = rho[mask]

    # GCF table is extrapolated nonphysically below ~0.25 GeV/c (S diverges as 1/p^4),
    # so we mask out events below the SRC validity threshold.
    p_rel_min_valid = 0.25
    valid = p_vals >= p_rel_min_valid
    p_vals = p_vals[valid]
    r_vals = r_vals[valid]
    nbins = 75
    bin_edges = np.linspace(0.0, 1.5, nbins + 1)
    centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    sum_rho, _ = np.histogram(p_vals, bins=bin_edges, weights=r_vals)
    counts, _ = np.histogram(p_vals, bins=bin_edges)
    mean_rho = np.divide(sum_rho, counts, out=np.zeros_like(sum_rho), where=counts > 0)

    # Radial probability density: n(p) = 4π p² S(p) / (2π)³
    n_p = 4.0 * np.pi * centers**2 * mean_rho / (2.0 * np.pi)**3
    # Normalize to 1 (∫ n(p) dp = 1)
    bin_width = bin_edges[1] - bin_edges[0]
    integral = np.sum(n_p) * bin_width
    if integral > 0:
        n_p = n_p / integral

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(centers[counts > 0], n_p[counts > 0], 'o-', markersize=3, color='tab:blue')
    ax.set_xlabel(r"$|p_{\rm rel}|$ (true initial relative momentum) [GeV/c]")
    ax.set_ylabel(r"$n(p_{\rm rel}) = \frac{4\pi\, p^2\, S(p)}{(2\pi)^3}$  [1/(GeV/c)]  (normalized)")
    ax.set_title(r"2N+FSI: Initial relative momentum probability density (normalized)")
    ax.grid(True, which='both', alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "wavefunction_vs_prel.png"), dpi=200)
    plt.close(fig)
    print(f"  wavefunction_vs_prel.png  (normalized so ∫n(p)dp = 1)")


def plot_momenta_N3(d, out_dir):
    """Momentum distributions for N=3 events."""
    mask = (d["nAboveKF"] >= 3) & d["has_p3_fsi"]
    w = d["weight"][mask]
    if np.sum(mask) == 0:
        print("  No N=3 events found, skipping momentum plots")
        return

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    for ax, vals, label in [
        (axes[0], d["p1_mag"][mask], r"$|p_1|$ (pmiss)"),
        (axes[1], d["p2_mag"][mask], r"$|p_2|$ (recoil)"),
        (axes[2], d["p3_fsi_mag"][mask], r"$|p_3|$ (FSI sec.)"),
    ]:
        plot_1d(ax, vals, w, 45, (0, 6), label="Post-FSI", linewidth=1.5, color='tab:blue')
        ax.set_xlabel(f"{label} [GeV/c]")
        ax.set_ylabel("Normalized weight")
        ax.legend()
    fig.suptitle("N=3 Momentum Distributions (2N+FSI)", fontsize=14)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "momenta_N3.png"), dpi=200)
    plt.close(fig)
    print("  momenta_N3.png")


def _signed_out_of_plane_angle(p, a, b):
    """Signed angle (deg) between vector p and the plane spanned by a and b."""
    normal = np.cross(a, b)
    n_mag = np.linalg.norm(normal, axis=1)
    p_mag = np.linalg.norm(p, axis=1)
    dot = np.sum(p * normal, axis=1)
    sin_angle = dot / (p_mag * n_mag + 1e-30)
    return np.degrees(np.arcsin(np.clip(sin_angle, -1, 1)))


def plot_coplanarity_N3(d, out_dir):
    """Out-of-plane angle for each nucleon in N=3 events (coplanarity measure)."""
    mask = (d["nAboveKF"] >= 3) & d["has_p3_fsi"]
    w = d["weight"][mask]
    if np.sum(mask) == 0:
        print("  No N=3 events found, skipping coplanarity plots")
        return

    # Nucleon momenta: pmiss for lead, post-FSI for recoil and FSI secondary
    q = np.column_stack([d["q"][mask, 0], d["q"][mask, 1], d["q"][mask, 2]])
    lp_post = np.column_stack([d["lead_post"][mask, 0], d["lead_post"][mask, 1], d["lead_post"][mask, 2]])
    pmiss = lp_post - q  # reconstructed initial lead momentum
    rp = np.column_stack([d["recoil_post"][mask, 0], d["recoil_post"][mask, 1], d["recoil_post"][mask, 2]])
    p3 = np.column_stack([d["p3_fsi_px"][mask], d["p3_fsi_py"][mask], d["p3_fsi_pz"][mask]])

    angle_lead = _signed_out_of_plane_angle(pmiss, rp, p3)
    angle_recoil = _signed_out_of_plane_angle(rp, pmiss, p3)
    angle_p3 = _signed_out_of_plane_angle(p3, pmiss, rp)

    # Fraction of events with all 3 angles below threshold
    w_total = np.sum(w)
    for thresh in (5, 10):
        cop_mask = (np.abs(angle_lead) < thresh) & (np.abs(angle_recoil) < thresh) & (np.abs(angle_p3) < thresh)
        pct = 100 * np.sum(w[cop_mask]) / w_total if w_total > 0 else 0
        print(f"  All coplanarity angles < {thresh}°: {np.sum(cop_mask)}/{len(w)}"
              f"  (weighted: {pct:.2f}%)")

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    for ax, vals, title in [
        (axes[0], angle_lead, r"Lead ($p_{\rm miss}$)"),
        (axes[1], angle_recoil, "Recoil nucleon"),
        (axes[2], angle_p3, "FSI secondary"),
    ]:
        plot_1d(ax, vals, w, 45, (-90, 90), linewidth=1.5, color='tab:blue')
        ax.set_xlabel("Out-of-plane angle [deg]", fontsize=12)
        ax.set_ylabel("Normalized weight", fontsize=12)
        ax.set_title(title)
    fig.suptitle("2N+FSI N=3: signed out-of-plane angle (0° = coplanar)", fontsize=12)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "coplanarity_N3.png"), dpi=200)
    plt.close(fig)
    print("  coplanarity_N3.png")


def _angle_between_planes(shared, a, b):
    """Angle (deg) between planes {shared, a} and {shared, b}, in [0, 90]."""
    n1 = np.cross(shared, a)
    n2 = np.cross(shared, b)
    n1_mag = np.linalg.norm(n1, axis=1, keepdims=True)
    n2_mag = np.linalg.norm(n2, axis=1, keepdims=True)
    n1_hat = n1 / (n1_mag + 1e-30)
    n2_hat = n2 / (n2_mag + 1e-30)
    cos_angle = np.clip(np.abs(np.sum(n1_hat * n2_hat, axis=1)), 0, 1)
    return np.degrees(np.arccos(cos_angle))


def plot_interplane_angle_N3(d, out_dir):
    """Angle between planes {pmiss,p2} and {pmiss,p3} for N=3 events."""
    mask = (d["nAboveKF"] >= 3) & d["has_p3_fsi"] & (d["Q2"] >= Q2min) & (d["scattering_angle"] < 45)
    w = d["weight"][mask]
    if np.sum(mask) == 0:
        print("  No N=3 events found, skipping inter-plane angle plot")
        return

    q = np.column_stack([d["q"][mask, 0], d["q"][mask, 1], d["q"][mask, 2]])
    lp_post = np.column_stack([d["lead_post"][mask, 0], d["lead_post"][mask, 1], d["lead_post"][mask, 2]])
    pmiss = lp_post - q
    rp = np.column_stack([d["recoil_post"][mask, 0], d["recoil_post"][mask, 1], d["recoil_post"][mask, 2]])
    p3 = np.column_stack([d["p3_fsi_px"][mask], d["p3_fsi_py"][mask], d["p3_fsi_pz"][mask]])

    angle = _angle_between_planes(pmiss, rp, p3)

    fig, ax = plt.subplots(figsize=(7, 5))
    plot_1d(ax, angle, w, 45, (0, 90), linewidth=1.5, color='tab:blue')
    ax.set_xlabel(r"Angle between planes $\{p_{\rm miss}, p_2\}$ and $\{p_{\rm miss}, p_3\}$ [deg]")
    ax.set_ylabel("Normalized weight")
    ax.set_title(r"2N+FSI N=3: inter-plane angle ($Q^2 \geq 1$ GeV$^2$, $\theta_e < 45^\circ$)")
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "interplane_angle_N3.png"), dpi=200)
    plt.close(fig)
    print("  interplane_angle_N3.png")


def plot_scattering_plane_angles_N3(d, out_dir):
    """Angles between scattering plane {pmiss,q} and nucleon planes {pmiss,p2}, {pmiss,p3} for N=3."""
    base_mask = (d["nAboveKF"] >= 3) & d["has_p3_fsi"] & (d["Q2"] >= Q2min) & (d["scattering_angle"] < 45)
    if np.sum(base_mask) == 0:
        print("  No N=3 events found, skipping scattering-plane angle plots")
        return

    # Full N=3 cuts (same as theta heatmap)
    angle_pq_max = 8
    lead_q_min = 0.75
    pmiss_lo, pmiss_hi = 0.25, 0.9
    xB_max = 1.2
    interplane_cut = 20
    full_cuts = ((d["p1_after_angle_q"] < angle_pq_max)
                 & (d["lead_over_q"] > lead_q_min)
                 & (pmiss_lo < d["pmiss"]) & (d["pmiss"] < pmiss_hi)
                 & (d["xB"] < xB_max)
                 & (d["interplane_angle_n3"] < interplane_cut))

    for cut_label, extra_mask, suffix in [
        ("no extra cuts", np.ones(len(d["weight"]), dtype=bool), ""),
        (rf"$\theta_{{pq}}<{angle_pq_max}^\circ$, $p_N/q>{lead_q_min}$, ${pmiss_lo}<p_{{\rm miss}}<{pmiss_hi}$, $x_B<{xB_max}$, $\phi_{{\rm ip}}<{interplane_cut}^\circ$",
         full_cuts, "_cuts"),
    ]:
        mask = base_mask & extra_mask
        w = d["weight"][mask]
        print(f"    [{suffix or 'base'}] events: {np.sum(mask)}")

        q = np.column_stack([d["q"][mask, 0], d["q"][mask, 1], d["q"][mask, 2]])
        lp_post = np.column_stack([d["lead_post"][mask, 0], d["lead_post"][mask, 1], d["lead_post"][mask, 2]])
        pmiss = lp_post - q
        rp = np.column_stack([d["recoil_post"][mask, 0], d["recoil_post"][mask, 1], d["recoil_post"][mask, 2]])
        p3 = np.column_stack([d["p3_fsi_px"][mask], d["p3_fsi_py"][mask], d["p3_fsi_pz"][mask]])

        angle_q_p2 = _angle_between_planes(pmiss, q, rp)
        angle_q_p3 = _angle_between_planes(pmiss, q, p3)

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

        plot_1d(ax1, angle_q_p2, w, 50, (0, 90), linewidth=1.5, color='tab:blue')
        ax1.set_xlabel(r"Angle between planes $\{p_{\rm miss}, \vec{q}\}$ and $\{p_{\rm miss}, p_2\}$ [deg]")
        ax1.set_ylabel("Normalized weight")
        ax1.set_title(r"2N+FSI N=3: $\{p_{\rm miss}, \vec{q}\}$ vs $\{p_{\rm miss}, p_2\}$")

        plot_1d(ax2, angle_q_p3, w, 50, (0, 90), linewidth=1.5, color='tab:red')
        ax2.set_xlabel(r"Angle between planes $\{p_{\rm miss}, \vec{q}\}$ and $\{p_{\rm miss}, p_3\}$ [deg]")
        ax2.set_ylabel("Normalized weight")
        ax2.set_title(r"2N+FSI N=3: $\{p_{\rm miss}, \vec{q}\}$ vs $\{p_{\rm miss}, p_3\}$")

        fig.suptitle(rf"2N+FSI N=3 — scattering-plane angles ({cut_label})", fontsize=13)
        fig.tight_layout()
        fname = f"scattering_plane_angles_N3{suffix}.png"
        fig.savefig(os.path.join(out_dir, fname), dpi=200)
        plt.close(fig)
        print(f"  {fname}")

        # 2D heatmap only for the full-cuts version
        if suffix == "_cuts":
            from matplotlib.colors import LogNorm
            fig2, ax2d = plt.subplots(figsize=(7, 6))
            h, xe, ye = np.histogram2d(angle_q_p2, angle_q_p3, bins=50,
                                        range=[[0, 90], [0, 90]], weights=w)
            h[h == 0] = np.nan
            im = ax2d.pcolormesh(xe, ye, h.T)
            fig2.colorbar(im, ax=ax2d, label="Weight")
            ax2d.set_xlabel(r"$\angle(\{p_{\rm miss},\vec{q}\}, \{p_{\rm miss},p_2\})$ [deg]")
            ax2d.set_ylabel(r"$\angle(\{p_{\rm miss},\vec{q}\}, \{p_{\rm miss},p_3\})$ [deg]")
            ax2d.set_title("2N+FSI N=3 — scattering-plane angle heatmap (with cuts)")
            ax2d.set_aspect('equal')
            fig2.tight_layout()
            fig2.savefig(os.path.join(out_dir, "scattering_plane_heatmap_N3.png"), dpi=200)
            plt.close(fig2)
            print("  scattering_plane_heatmap_N3.png")


def plot_interplane_angle_3N_generator(d3, out_dir):
    """Angle between planes {pmiss,p2} and {pmiss,p3} for 3N generator."""
    w = d3["weight"]
    q = np.column_stack([d3["q"][:, 0], d3["q"][:, 1], d3["q"][:, 2]])
    lead = np.column_stack([d3["lead"][:, 0], d3["lead"][:, 1], d3["lead"][:, 2]])
    pmiss = lead - q
    p2 = np.column_stack([d3["recoil2"][:, 0], d3["recoil2"][:, 1], d3["recoil2"][:, 2]])
    p3 = np.column_stack([d3["recoil3"][:, 0], d3["recoil3"][:, 1], d3["recoil3"][:, 2]])

    angle = _angle_between_planes(pmiss, p2, p3)

    fig, ax = plt.subplots(figsize=(7, 5))
    plot_1d(ax, angle, w, 45, (0, 90), linewidth=1.5, color='tab:orange')
    ax.set_xlabel(r"Angle between planes $\{p_{\rm miss}, p_2\}$ and $\{p_{\rm miss}, p_3\}$ [deg]")
    ax.set_ylabel("Normalized weight")
    ax.set_title(r"3N Generator: inter-plane angle (sanity check, expect $\delta$(0°))")
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "interplane_angle_3N_generator.png"), dpi=200)
    plt.close(fig)
    print("  interplane_angle_3N_generator.png")


def plot_interplane_ratio_3N_over_FSI(d, d3, out_dir):
    """Plot cumulative weight ratio 3N/FSI as a function of interplane angle cut."""
    # 3N generator: interplane angle with Q2 and theta_e cuts
    mask_3N = (d3["Q2"] >= Q2min) & (d3["scattering_angle"] < 45)
    q3 = np.column_stack([d3["q"][mask_3N, 0], d3["q"][mask_3N, 1], d3["q"][mask_3N, 2]])
    lead3 = np.column_stack([d3["lead"][mask_3N, 0], d3["lead"][mask_3N, 1], d3["lead"][mask_3N, 2]])
    pmiss3 = lead3 - q3
    p2_3 = np.column_stack([d3["recoil2"][mask_3N, 0], d3["recoil2"][mask_3N, 1], d3["recoil2"][mask_3N, 2]])
    p3_3 = np.column_stack([d3["recoil3"][mask_3N, 0], d3["recoil3"][mask_3N, 1], d3["recoil3"][mask_3N, 2]])
    angle_3N = _angle_between_planes(pmiss3, p2_3, p3_3)
    w_3N = d3["weight"][mask_3N]

    # 2N+FSI N=3: interplane angle with Q2, theta_e, nAboveKF>=3 cuts
    mask_FSI = (d["nAboveKF"] >= 3) & d["has_p3_fsi"] & (d["Q2"] >= Q2min) & (d["scattering_angle"] < 45)
    q_f = np.column_stack([d["q"][mask_FSI, 0], d["q"][mask_FSI, 1], d["q"][mask_FSI, 2]])
    lp_f = np.column_stack([d["lead_post"][mask_FSI, 0], d["lead_post"][mask_FSI, 1], d["lead_post"][mask_FSI, 2]])
    pmiss_f = lp_f - q_f
    rp_f = np.column_stack([d["recoil_post"][mask_FSI, 0], d["recoil_post"][mask_FSI, 1], d["recoil_post"][mask_FSI, 2]])
    p3_f = np.column_stack([d["p3_fsi_px"][mask_FSI], d["p3_fsi_py"][mask_FSI], d["p3_fsi_pz"][mask_FSI]])
    angle_FSI = _angle_between_planes(pmiss_f, rp_f, p3_f)
    w_FSI = d["weight"][mask_FSI]

    # Scan phi_cut and compute cumulative weight ratio
    phi_cuts = np.linspace(1, 90, 200)
    total_3N = np.sum(w_3N)
    total_FSI = np.sum(w_FSI)
    ratios = []
    for phi_c in phi_cuts:
        cum_3N = np.sum(w_3N[angle_3N < phi_c]) / total_3N
        cum_FSI = np.sum(w_FSI[angle_FSI < phi_c]) / total_FSI
        ratios.append(cum_3N / cum_FSI if cum_FSI > 0 else 0)
    ratios = np.array(ratios)

    best_idx = np.argmax(ratios)
    best_phi = phi_cuts[best_idx]
    best_ratio = ratios[best_idx]

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(phi_cuts, ratios, 'k-', linewidth=1.5)
    ax.axvline(best_phi, color='red', linestyle='--', linewidth=1,
               label=rf"max at $\phi_{{\rm cut}}={best_phi:.1f}^\circ$ (ratio={best_ratio:.2f})")
    ax.set_xlabel(r"$\phi_{\rm cut}$ [deg]")
    ax.set_ylabel(r"$\sum w_{3N}(\phi<\phi_{\rm cut}) \;/\; \sum w_{\rm FSI}(\phi<\phi_{\rm cut})$  (both normalized)")
    ax.set_title(r"Cumulative weight ratio: 3N / (2N+FSI $N\geq3$) vs interplane angle cut")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "interplane_ratio_3N_over_FSI.png"), dpi=200)
    plt.close(fig)
    print(f"  interplane_ratio_3N_over_FSI.png  (best phi_cut={best_phi:.1f} deg, ratio={best_ratio:.2f})")


def plot_scatplane_ratio_3N_over_FSI(d, d3, out_dir):
    """Cumulative weight ratio 3N/FSI as function of scattering-plane angle cut (with full SRC cuts)."""
    # Cut parameters (same as theta heatmap N=3)
    angle_pq_max = 8
    lead_q_min = 0.75
    pmiss_lo, pmiss_hi = 0.25, 0.9
    xB_max = 1.2
    interplane_cut = 20

    # --- 3N generator ---
    mask_3N = ((d3["Q2"] >= Q2min) & (d3["scattering_angle"] < 45)
               & (d3["p1_after_angle_q"] < angle_pq_max)
               & (d3["lead_over_q"] > lead_q_min)
               & (pmiss_lo < d3["p1_mag"]) & (d3["p1_mag"] < pmiss_hi)
               & (d3["xB"] < xB_max)
               & (d3["interplane_angle"] < interplane_cut))
    q3 = np.column_stack([d3["q"][mask_3N, 0], d3["q"][mask_3N, 1], d3["q"][mask_3N, 2]])
    lead3 = np.column_stack([d3["lead"][mask_3N, 0], d3["lead"][mask_3N, 1], d3["lead"][mask_3N, 2]])
    pmiss3 = lead3 - q3
    p2_3N = np.column_stack([d3["recoil2"][mask_3N, 0], d3["recoil2"][mask_3N, 1], d3["recoil2"][mask_3N, 2]])
    p3_3N = np.column_stack([d3["recoil3"][mask_3N, 0], d3["recoil3"][mask_3N, 1], d3["recoil3"][mask_3N, 2]])
    phi1_3N = _angle_between_planes(pmiss3, q3, p2_3N)
    phi2_3N = _angle_between_planes(pmiss3, q3, p3_3N)
    w_3N = d3["weight"][mask_3N]

    # --- 2N+FSI N=3 ---
    mask_FSI = ((d["nAboveKF"] >= 3) & d["has_p3_fsi"]
                & (d["Q2"] >= Q2min) & (d["scattering_angle"] < 45)
                & (d["p1_after_angle_q"] < angle_pq_max)
                & (d["lead_over_q"] > lead_q_min)
                & (pmiss_lo < d["pmiss"]) & (d["pmiss"] < pmiss_hi)
                & (d["xB"] < xB_max)
                & (d["interplane_angle_n3"] < interplane_cut))
    q_f = np.column_stack([d["q"][mask_FSI, 0], d["q"][mask_FSI, 1], d["q"][mask_FSI, 2]])
    lp_f = np.column_stack([d["lead_post"][mask_FSI, 0], d["lead_post"][mask_FSI, 1], d["lead_post"][mask_FSI, 2]])
    pmiss_f = lp_f - q_f
    rp_f = np.column_stack([d["recoil_post"][mask_FSI, 0], d["recoil_post"][mask_FSI, 1], d["recoil_post"][mask_FSI, 2]])
    p3_f = np.column_stack([d["p3_fsi_px"][mask_FSI], d["p3_fsi_py"][mask_FSI], d["p3_fsi_pz"][mask_FSI]])
    phi1_FSI = _angle_between_planes(pmiss_f, q_f, rp_f)
    phi2_FSI = _angle_between_planes(pmiss_f, q_f, p3_f)
    w_FSI = d["weight"][mask_FSI]

    # Scan phi_cut: require BOTH phi1 < phi_cut AND phi2 < phi_cut
    phi_cuts = np.linspace(1, 90, 200)
    total_3N = np.sum(w_3N)
    total_FSI = np.sum(w_FSI)
    ratios = []
    for phi_c in phi_cuts:
        sel_3N = (phi1_3N < phi_c) & (phi2_3N < phi_c)
        sel_FSI = (phi1_FSI < phi_c) & (phi2_FSI < phi_c)
        cum_3N = np.sum(w_3N[sel_3N]) / total_3N if total_3N > 0 else 0
        cum_FSI = np.sum(w_FSI[sel_FSI]) / total_FSI if total_FSI > 0 else 0
        ratios.append(cum_3N / cum_FSI if cum_FSI > 0 else 0)
    ratios = np.array(ratios)

    best_idx = np.argmax(ratios)
    best_phi = phi_cuts[best_idx]
    best_ratio = ratios[best_idx]

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(phi_cuts, ratios, 'k-', linewidth=1.5)
    ax.axvline(best_phi, color='red', linestyle='--', linewidth=1,
               label=rf"max at $\phi_{{\rm cut}}={best_phi:.1f}^\circ$ (ratio={best_ratio:.2f})")
    ax.set_xlabel(r"$\phi_{\rm cut}$ [deg]")
    ax.set_ylabel(r"$\frac{\sum w_{3N}(\phi_1,\phi_2<\phi_{\rm cut})}{\sum w_{\rm FSI}(\phi_1,\phi_2<\phi_{\rm cut})}$  (both normalized)",
                   fontsize=11)
    ax.set_title(r"Cumulative weight ratio: 3N / FSI vs scattering-plane angle cut (with SRC cuts)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "scatplane_ratio_3N_over_FSI.png"), dpi=200)
    plt.close(fig)
    print(f"  scatplane_ratio_3N_over_FSI.png  (best phi_cut={best_phi:.1f} deg, ratio={best_ratio:.2f})")

    # S^2/B version
    ratios_s2b = []
    for phi_c in phi_cuts:
        sel_3N = (phi1_3N < phi_c) & (phi2_3N < phi_c)
        sel_FSI = (phi1_FSI < phi_c) & (phi2_FSI < phi_c)
        cum_3N = np.sum(w_3N[sel_3N]) / total_3N if total_3N > 0 else 0
        cum_FSI = np.sum(w_FSI[sel_FSI]) / total_FSI if total_FSI > 0 else 0
        ratios_s2b.append(cum_3N**2 / cum_FSI if cum_FSI > 0 else 0)
    ratios_s2b = np.array(ratios_s2b)

    best_idx2 = np.argmax(ratios_s2b)
    best_phi2 = phi_cuts[best_idx2]
    best_ratio2 = ratios_s2b[best_idx2]

    fig2, ax2 = plt.subplots(figsize=(8, 5))
    ax2.plot(phi_cuts, ratios_s2b, 'k-', linewidth=1.5)
    ax2.axvline(best_phi2, color='red', linestyle='--', linewidth=1,
                label=rf"max at $\phi_{{\rm cut}}={best_phi2:.1f}^\circ$ (ratio={best_ratio2:.3f})")
    ax2.set_xlabel(r"$\phi_{\rm cut}$ [deg]")
    ax2.set_ylabel(r"$\frac{(\sum w_{3N})^2}{\sum w_{\rm FSI}}$  (both normalized)", fontsize=11)
    ax2.set_title(r"$S^2/B$: 3N$^2$ / FSI vs scattering-plane angle cut (with SRC cuts)")
    ax2.legend()
    fig2.tight_layout()
    fig2.savefig(os.path.join(out_dir, "scatplane_s2b_3N_over_FSI.png"), dpi=200)
    plt.close(fig2)
    print(f"  scatplane_s2b_3N_over_FSI.png  (best phi_cut={best_phi2:.1f} deg, ratio={best_ratio2:.3f})")


def _compute_lightcone_pD(d_in, is_fsi=False):
    """Compute lightcone variables for pD topology. Returns (mask, alpha_q, alpha_p_final, alpha_d, alpha_p_initial, alpha_sum, weights)."""
    w = d_in["weight"]

    if is_fsi:
        # FSI: p2 = recoil_post, p3 = p3_fsi
        r2 = d_in["recoil_post"][:, :3]
        r3 = np.column_stack([d_in["p3_fsi_px"], d_in["p3_fsi_py"], d_in["p3_fsi_pz"]])
        lead = d_in["lead_post"][:, :3]
        q_3 = d_in["q"][:, :3]
        pmiss = lead - q_3
        # pn check: rec_type and p3_fsi_pdg
        is_pn = (((d_in["rec_type"] == 2212) & (d_in["p3_fsi_pdg"] == 2112)) |
                 ((d_in["rec_type"] == 2112) & (d_in["p3_fsi_pdg"] == 2212)))
        lead_is_proton = d_in["lead_type"] == 2212
        n3_mask = (d_in["nAboveKF"] >= 3) & d_in["has_p3_fsi"]
    else:
        # 3N generator
        r2 = d_in["recoil2"][:, :3]
        r3 = d_in["recoil3"][:, :3]
        lead = d_in["lead"][:, :3]
        q_3 = d_in["q"][:, :3]
        pmiss = lead - q_3
        pCode, nCode = 2212, 2112
        is_pn = (((d_in["N2_type"] == pCode) & (d_in["N3_type"] == nCode)) |
                 ((d_in["N2_type"] == nCode) & (d_in["N3_type"] == pCode)))
        lead_is_proton = d_in["N1_type"] == 2212
        n3_mask = np.ones(len(w), dtype=bool)

    pd_vec = r2 + r3
    pd_mag = np.linalg.norm(pd_vec, axis=1)
    pmiss_mag = np.linalg.norm(pmiss, axis=1)

    # Relative momentum
    p_rel = np.linalg.norm(r2 - r3, axis=1) / 2.0
    p_rel_max = 0.02
    is_deuteron = is_pn & (p_rel < p_rel_max)

    # lead kinematics
    lead_mag = np.linalg.norm(lead, axis=1)
    q_mag_all = np.linalg.norm(q_3, axis=1)
    lead_over_q = lead_mag / (q_mag_all + 1e-30)
    dot_lead_q = np.sum(lead * q_3, axis=1)
    angle_lead_q = np.degrees(np.arccos(np.clip(dot_lead_q / (lead_mag * q_mag_all + 1e-30), -1, 1)))

    pd_theta = np.degrees(np.arccos(np.clip(pd_vec[:, 2] / (pd_mag + 1e-30), -1, 1)))

    e_angle_cut = (d_in["scattering_angle"] > 7) & (d_in["scattering_angle"] < 45)

    mask = (n3_mask & is_deuteron & lead_is_proton
            & (0.55 < pd_mag) & (pd_mag < 0.9)
            & (d_in["xB"] > 1.3) & (d_in["Q2"] > 1.0)
            & (0.65 < lead_over_q) & (lead_over_q < 0.95)
            & (angle_lead_q < 30.0)
            & (50.0 < pd_theta) & (pd_theta < 110.0)
            & e_angle_cut)

    if np.sum(mask) == 0:
        return mask, None, None, None, None, None, None

    q_3m = q_3[mask]
    q_mag = np.linalg.norm(q_3m, axis=1)
    q_hat = q_3m / (q_mag[:, None] + 1e-30)
    nu = d_in["nu"][mask]

    alpha_q = (nu - q_mag) / mN

    p_f = lead[mask]
    p_f_mag = np.linalg.norm(p_f, axis=1)
    E_p = np.sqrt(mN**2 + p_f_mag**2)
    p_f_proj = np.sum(p_f * q_hat, axis=1)
    alpha_p_final = (E_p - p_f_proj) / mN

    m_d = 1.875613
    pd_m = pd_vec[mask]
    pd_mag_m = pd_mag[mask]
    pd_proj = np.sum(pd_m * q_hat, axis=1)
    E_d = np.sqrt(m_d**2 + pd_mag_m**2)
    alpha_deuteron = (E_d - pd_proj) / mN

    alpha_p_initial = alpha_p_final - alpha_q
    alpha_sum = alpha_p_initial + alpha_deuteron

    return mask, alpha_q, alpha_p_final, alpha_deuteron, alpha_p_initial, alpha_sum, w[mask]


def plot_lightcone_pD_N3(d, out_dir):
    """Lightcone variables for 2N+FSI N=3 pD topology."""
    mask, alpha_q, alpha_p_final, alpha_d, alpha_p_initial, alpha_sum, w_m = _compute_lightcone_pD(d, is_fsi=True)
    if alpha_q is None:
        print("  No pD topology events in FSI, skipping")
        return

    variables = [
        (alpha_q, r"$\alpha_q = (\nu - |\vec{q}|)/m_N$", (-1.6, -0.2)),
        (alpha_p_final, r"$\alpha_{p}^{\rm final}$", (0.0, 1.0)),
        (alpha_d, r"$\alpha_{d}$", (1.0, 2.0)),
        (alpha_p_initial, r"$\alpha_{p}^{\rm initial}$", (1.0, 2.0)),
        (alpha_sum, r"$\alpha_{p}^{\rm initial} + \alpha_{d}$", (1.5, 4.0)),
    ]

    fig, axes = plt.subplots(len(variables), 1, figsize=(8, 3 * len(variables)))
    for ax, (vals, label, rng) in zip(axes, variables):
        counts, edges = np.histogram(vals, bins=80, range=rng, weights=w_m)
        centers = 0.5 * (edges[:-1] + edges[1:])
        ax.step(centers, counts, where='mid', linewidth=1.2, color='tab:blue', label='Total')
        wmean = np.average(vals, weights=w_m)
        ax.text(0.98, 0.95, rf"$\langle {label[1:-1]} \rangle = {wmean:.3f}$",
                transform=ax.transAxes, ha='right', va='top', fontsize=10,
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        ax.set_xlabel(label)
        ax.set_ylabel("Total Weight")
        ax.legend(fontsize=8)
    fig.suptitle(rf"2N+FSI N=3: Lightcone — pD topology ({np.sum(mask)} events)", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "lightcone_pD_N3.png"), dpi=200)
    plt.close(fig)
    print(f"  lightcone_pD_N3.png  ({np.sum(mask)} events)")


def plot_lightcone_pD_combined(d, d3, out_dir, lam=0.5):
    """Combined 3N + FSI lightcone distributions with tunable lambda."""
    mask_3N, aq_3N, apf_3N, ad_3N, api_3N, asum_3N, w_3N = _compute_lightcone_pD(d3, is_fsi=False)
    mask_FSI, aq_FSI, apf_FSI, ad_FSI, api_FSI, asum_FSI, w_FSI = _compute_lightcone_pD(d, is_fsi=True)

    if aq_3N is None or aq_FSI is None:
        print("  Insufficient events for combined lightcone plot, skipping")
        return

    variables = [
        (aq_3N, aq_FSI, r"$\alpha_q$", (-1.6, -0.2)),
        (apf_3N, apf_FSI, r"$\alpha_{p}^{\rm final}$", (0.0, 1.0)),
        (ad_3N, ad_FSI, r"$\alpha_{d}$", (1.0, 2.0)),
        (api_3N, api_FSI, r"$\alpha_{p}^{\rm initial}$", (1.0, 2.0)),
        (asum_3N, asum_FSI, r"$\alpha_{p}^{\rm initial} + \alpha_{d}$", (1.5, 4.0)),
    ]

    fig, axes = plt.subplots(len(variables), 1, figsize=(8, 3 * len(variables)))
    for ax, (v3N, vFSI, label, rng) in zip(axes, variables):
        nbins = 80
        h3N, edges = np.histogram(v3N, bins=nbins, range=rng, weights=w_3N)
        hFSI, _ = np.histogram(vFSI, bins=nbins, range=rng, weights=w_FSI)
        centers = 0.5 * (edges[:-1] + edges[1:])

        # Normalize each to 1
        h3N_n = h3N / (np.sum(h3N) + 1e-30)
        hFSI_n = hFSI / (np.sum(hFSI) + 1e-30)
        h_comb = lam * h3N_n + (1 - lam) * hFSI_n

        ax.step(centers, h3N_n, where='mid', linewidth=1, color='tab:orange', alpha=0.4, label='3N')
        ax.step(centers, hFSI_n, where='mid', linewidth=1, color='tab:blue', alpha=0.4, label='2N+FSI')
        ax.step(centers, h_comb, where='mid', linewidth=2, color='black',
                label=rf'$\lambda={lam:.1f}$')
        ax.set_xlabel(label)
        ax.set_ylabel("Normalized")
        ax.legend(fontsize=8)

    fig.suptitle(rf"Lightcone pD: combined $\lambda \cdot$3N + $(1-\lambda) \cdot$FSI ($\lambda={lam}$)", fontsize=12)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "lightcone_pD_combined.png"), dpi=200)
    plt.close(fig)
    print(f"  lightcone_pD_combined.png  (lambda={lam})")


def plot_coplanarity_3N_generator(d3, out_dir):
    """Out-of-plane angle for 3N generator events (sanity check: should be delta at 0)."""
    w = d3["weight"]

    # pmiss = lead - q, recoil2, recoil3
    q = np.column_stack([d3["q"][:, 0], d3["q"][:, 1], d3["q"][:, 2]])
    lead = np.column_stack([d3["lead"][:, 0], d3["lead"][:, 1], d3["lead"][:, 2]])
    pmiss = lead - q
    p2 = np.column_stack([d3["recoil2"][:, 0], d3["recoil2"][:, 1], d3["recoil2"][:, 2]])
    p3 = np.column_stack([d3["recoil3"][:, 0], d3["recoil3"][:, 1], d3["recoil3"][:, 2]])

    angle_lead = _signed_out_of_plane_angle(pmiss, p2, p3)
    angle_recoil2 = _signed_out_of_plane_angle(p2, pmiss, p3)
    angle_recoil3 = _signed_out_of_plane_angle(p3, pmiss, p2)

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    for ax, vals, title in [
        (axes[0], angle_lead, r"Lead ($p_{\rm miss}$)"),
        (axes[1], angle_recoil2, "Recoil 2"),
        (axes[2], angle_recoil3, "Recoil 3"),
    ]:
        plot_1d(ax, vals, w, 45, (-90, 90), linewidth=1.5, color='tab:orange')
        ax.set_xlabel("Out-of-plane angle [deg]", fontsize=12)
        ax.set_ylabel("Normalized weight", fontsize=12)
        ax.set_title(title)
    fig.suptitle("3N Generator: signed out-of-plane angle (sanity check, expect δ(0))", fontsize=12)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "coplanarity_3N_generator.png"), dpi=200)
    plt.close(fig)
    print("  coplanarity_3N_generator.png")


def plot_SRC_momenta_with_cuts(d, out_dir):
    """TEMPORARY: Lead and recoil momentum distributions under SRC-style cuts."""
    q_mag = np.sqrt(d["q"][:, 0]**2 + d["q"][:, 1]**2 + d["q"][:, 2]**2)

    mask = (
        (d["xB"] > 1.2)
        & (d["Q2"] >= 1.7)
        & (d["Q2"] < 4.0)
        & (d["p1_after_angle_q"] < 25)
        & (d["lead_over_q"] > 0.62) & (d["lead_over_q"] < 0.92)
        & (d["p1_mag"] > 0.4) & (d["p1_mag"] < 1.0)
        & (d["mmiss"] <= 1.1)
        & (d["p2_mag"] > 0.35)
        & (d["scattering_angle"] < 45)
    )
    w = d["weight"][mask]
    n_pass = np.sum(mask)
    print(f"  Events passing SRC cuts: {n_pass}  (weighted: {np.sum(w):.2f})")

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Lead nucleon momentum
    ax = axes[0]
    plot_1d(ax, d["p1_after_mag"][mask], w, 50, (0, 2.5),
            label="Lead (FSI)", linewidth=1.5, color='tab:blue')
    plot_1d(ax, d["p1_pre_after_mag"][mask], w, 50, (0, 2.5),
            label="Lead (PWIA)", linewidth=1.5, color='tab:blue', linestyle='--')
    ax.set_xlabel(r"$|p_N|$ [GeV/c]", fontsize=12)
    ax.set_ylabel("Normalized weight", fontsize=12)
    ax.set_title("Lead nucleon momentum")
    ax.legend()

    # Recoil nucleon momentum
    ax = axes[1]
    plot_1d(ax, d["p2_mag"][mask], w, 50, (0.35, 1),
            label="Recoil (FSI)", linewidth=1.5, color='tab:red')
    plot_1d(ax, d["p2_pre_mag"][mask], w, 50, (0.35, 1),
            label="Recoil (PWIA)", linewidth=1.5, color='tab:red', linestyle='--')
    ax.set_xlabel(r"$|p_{\mathrm{recoil}}|$ [GeV/c]", fontsize=12)
    ax.set_ylabel("Normalized weight", fontsize=12)
    ax.set_title("Recoil nucleon momentum")
    ax.legend()

    fig.suptitle(r"2N+FSI: SRC cuts ($x_B>1.2$, $\theta_{pq}<25°$, $0.62<|p_N|/|q|<0.92$,"
                 "\n" r"$0.4<|p_{\rm miss}|<1.0$, $m_{\rm miss}\leq1.1$, $|p_{\rm rec}|>0.35$)",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "SRC_cuts_momenta.png"), dpi=200)
    plt.close(fig)
    print("  SRC_cuts_momenta.png")


# ──────────── Main ────────────
def main():
    parser = argparse.ArgumentParser(description="Analyze 2N SRC events from ROOT TTree")
    parser.add_argument("input", help="Path to events ROOT file (e.g. events/hN/events_2N_fsi_501.root)",
                        default="events/hN/events_2N_fsi_501.root")
    parser.add_argument("--input-3n", default="../genQE_3N/events_3N.root",
                        help="Path to events_3N.root for cross-generator comparisons")
    parser.add_argument("-o", "--output-dir", default="analysis_output_2N/png_files",
                        help="Output directory for plots")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    print(f"Loading {args.input}...")
    d = load_data(args.input)
    print(f"  {len(d['weight'])} events loaded")

    print("Computing derived quantities...")
    d = add_derived(d)
    d = find_p3_fsi(d)

    print(f"\nPlots will be saved to: {args.output_dir}\n")

    print("=== SRC-cut momentum distributions (TEMPORARY) ===")
    plot_SRC_momenta_with_cuts(d, args.output_dir)

    print("=== Theta heatmaps ===")
    plot_theta_heatmaps(d, args.output_dir)

    print("\n=== 1D distributions ===")
    plot_1d_distributions(d, args.output_dir)

    print("\n=== xB-Q2 heatmap ===")
    plot_xB_Q2_heatmap(d, args.output_dir)

    print("\n=== Region overlays (N=2 vs 3N, N=3 vs 3N) ===")
    try:
        print(f"Loading 3N data from {args.input_3n}...")
        d3 = load_data_3N(args.input_3n)
        d3 = add_derived_3N(d3)
        print(f"  {len(d3['weight'])} 3N events loaded")
        print("\n=== 3N / N=3 ratio heatmaps ===")
        plot_theta_ratio_heatmap(d, d3, args.output_dir)
        plot_region_overlays(d, d3, args.output_dir)
        print("\n=== 3N generator coplanarity (sanity check) ===")
        plot_coplanarity_3N_generator(d3, args.output_dir)
        print("\n=== 3N generator inter-plane angle (sanity check) ===")
        plot_interplane_angle_3N_generator(d3, args.output_dir)
        print("\n=== Interplane angle ratio: 3N / FSI ===")
        plot_interplane_ratio_3N_over_FSI(d, d3, args.output_dir)
        print("\n=== Scattering-plane angle ratio: 3N / FSI ===")
        plot_scatplane_ratio_3N_over_FSI(d, d3, args.output_dir)
        print("\n=== Lightcone pD combined (3N + FSI) ===")
        plot_lightcone_pD_combined(d, d3, args.output_dir, lam=0.1)
    except FileNotFoundError:
        print(f"  WARNING: 3N file not found ({args.input_3n}), skipping region overlays")

    print("\n=== Pair CM momentum (pre-FSI) ===")
    plot_pair_cm_momentum(d, args.output_dir)

    print("\n=== Wavefunction vs p_rel ===")
    plot_wavefunction_vs_prel(d, args.output_dir)

    print("\n=== N=3 momentum distributions ===")
    plot_momenta_N3(d, args.output_dir)

    print("\n=== N=3 coplanarity ===")
    plot_coplanarity_N3(d, args.output_dir)

    print("\n=== N=3 inter-plane angle ===")
    plot_interplane_angle_N3(d, args.output_dir)

    print("\n=== N=3 scattering-plane angles ===")
    plot_scattering_plane_angles_N3(d, args.output_dir)

    print("\n=== N=3 lightcone pD topology ===")
    plot_lightcone_pD_N3(d, args.output_dir)

    compute_region_fractions(d)

    print("\nDone.")


if __name__ == "__main__":
    main()
