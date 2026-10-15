#!/usr/bin/env python3
"""Side-by-side theta12-theta23 heatmap: 3N+FSI vs 2N+FSI.

Both panels use the same 2N-style angle conventions (p_miss = p_lead - q,
p3_recon = -(p_miss + p_recoil_measured)) and the same cut chain. The only
difference is the input sample.

Cuts (default; flags can override):
    weight > 0
    detector tag:
        - 3N: |p_lead| > kF AND at least one recoil above kF;
              if both above kF, pick the higher-|p| one as the measured recoil.
        - 2N: nAboveKF == 2 (both nucleons above kF; the only 2N-quality tag).
    theta_e < 45 deg
    Q^2 > 1.0 GeV^2
    xB < 1.2

The tighter SRC-trigger cuts (theta_pq, lead/q, |p_miss| range) are NOT
applied here -- this matches the user's "drop those three cuts" request.

Both panels share a single colour scale so visual comparison is direct.
"""
import argparse
import os
import sys

import numpy as np
import uproot
import awkward as ak
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paper_style import apply_style, figure_size, CMAP, OVERLAY_COLOR


kF = 0.25
mN = 0.93892

# Standard SRC-trigger cut values (for the pass-fraction report only)
SRC_THETA_PQ_MAX = 8.0
SRC_LEAD_OVER_Q_MIN = 0.75
SRC_PMISS_LO = 0.25
SRC_PMISS_HI = 0.90


def angle_deg(ax, ay, az, bx, by, bz):
    dot = ax * bx + ay * by + az * bz
    a = np.sqrt(ax * ax + ay * ay + az * az)
    b = np.sqrt(bx * bx + by * by + bz * bz)
    return np.degrees(np.arccos(np.clip(dot / (a * b + 1e-30), -1.0, 1.0)))


# ---------- 3N sample (raw genQE_3N_FSI tree) -------------------------------
def load_3N(path, ebeam, theta_e_max, Q2_min, xB_max, kF_thr,
            theta_pq_max=180.0, lead_over_q_min=0.0,
            pmiss_lo=0.0, pmiss_hi=1e9, mode='ge2',
            theta_e_min=0.0, stats=None,
            lead_over_q_max=1e9, xB_min=0.0, src_from_args=False):
    print(f'Reading 3N+FSI: {path}')
    with uproot.open(path) as f:
        tree = f['genT']
        arr = tree.arrays(['weight', 'pe', 'pLead', 'p2', 'p3'], library='np')
        # Secondaries added 2026-07-20 so 3N uses the identical detection rule
        # as 2N/MF. Older samples lack them -> fall back to primaries only.
        has_sec_branches = 'sec_px' in tree.keys()
        sec = (tree.arrays(['sec_px', 'sec_py', 'sec_pz', 'sec_pdg'], library='ak')
               if has_sec_branches else None)
    w  = arr['weight']
    pe = arr['pe']
    pL = arr['pLead']; p2 = arr['p2']; p3 = arr['p3']

    qx = -pe[:, 0]; qy = -pe[:, 1]; qz = ebeam - pe[:, 2]
    q_mag = np.sqrt(qx**2 + qy**2 + qz**2)
    pe_mag = np.sqrt((pe**2).sum(axis=1))
    omega = ebeam - pe_mag
    Q2 = q_mag**2 - omega**2
    xB = Q2 / (2.0 * mN * omega + 1e-30)
    theta_e = np.degrees(np.arctan2(np.sqrt(pe[:, 0]**2 + pe[:, 1]**2),
                                    pe[:, 2]))

    pL_mag = np.sqrt((pL**2).sum(axis=1))
    p2_mag = np.sqrt((p2**2).sum(axis=1))
    p3_mag = np.sqrt((p3**2).sum(axis=1))

    above_L  = pL_mag > kF_thr
    above_p2 = p2_mag > kF_thr
    above_p3 = p3_mag > kF_thr
    only_p2 = above_p2 & ~above_p3
    only_p3 = above_p3 & ~above_p2
    both    = above_p2 & above_p3

    # Leading above-kF secondary nucleon and how many exist (for the detection
    # rule + the nAboveKF count), identical to the 2N/MF treatment.
    n = len(w)
    if sec is not None:
        smag = np.sqrt(sec['sec_px']**2 + sec['sec_py']**2 + sec['sec_pz']**2)
        sec_above = ((sec['sec_pdg'] == 2212) | (sec['sec_pdg'] == 2112)) & (smag > kF_thr)
        n_sec_above = ak.to_numpy(ak.sum(sec_above, axis=1))
        skey = ak.where(sec_above, smag, -1.0)
        has_sec = ak.to_numpy(ak.fill_none(ak.max(skey, axis=1) > 0, False))
        sidx = ak.argmax(skey, axis=1, keepdims=True)
        sec_x = ak.to_numpy(ak.fill_none(ak.flatten(sec['sec_px'][sidx]), 0.0))
        sec_y = ak.to_numpy(ak.fill_none(ak.flatten(sec['sec_py'][sidx]), 0.0))
        sec_z = ak.to_numpy(ak.fill_none(ak.flatten(sec['sec_pz'][sidx]), 0.0))
    else:
        n_sec_above = np.zeros(n, dtype=int)
        has_sec = np.zeros(n, dtype=bool)
        sec_x = sec_y = sec_z = np.zeros(n)

    nAbove = (above_L.astype(int) + above_p2.astype(int)
              + above_p3.astype(int) + n_sec_above)

    if mode == 'eq2':
        # Exactly 2 nucleons above kF among {lead, p2, p3, secondaries}; the
        # lead is one, and the detected recoil is the UNIQUE other above-kF
        # nucleon -- a surviving primary if it stayed above kF, else the leading
        # above-kF ejected secondary. Same rule as load_2N.
        detector = above_L & (nAbove == 2)
        use_p2  = above_p2
        use_p3  = above_p3 & ~above_p2
        use_sec = (~above_p2) & (~above_p3) & has_sec
    elif mode == 'ge2':
        # lead above kF AND at least one recoil above kF; pick higher-|p|
        detector = above_L & (above_p2 | above_p3)
        use_p2 = detector & (only_p2 | (both & (p2_mag >= p3_mag)))
        use_p3 = detector & (only_p3 | (both & (p3_mag >  p2_mag)))
        use_sec = np.zeros(n, dtype=bool)
    else:
        raise ValueError(f"mode must be 'eq2' or 'ge2', got {mode!r}")

    rx = np.where(use_p2, p2[:, 0], np.where(use_p3, p3[:, 0], sec_x))
    ry = np.where(use_p2, p2[:, 1], np.where(use_p3, p3[:, 1], sec_y))
    rz = np.where(use_p2, p2[:, 2], np.where(use_p3, p3[:, 2], sec_z))

    p_miss_x = pL[:, 0] - qx
    p_miss_y = pL[:, 1] - qy
    p_miss_z = pL[:, 2] - qz
    p_miss_mag = np.sqrt(p_miss_x**2 + p_miss_y**2 + p_miss_z**2)

    p3rec_x = -(p_miss_x + rx)
    p3rec_y = -(p_miss_y + ry)
    p3rec_z = -(p_miss_z + rz)

    theta12 = angle_deg(p_miss_x, p_miss_y, p_miss_z, rx, ry, rz)
    theta23 = angle_deg(rx, ry, rz, p3rec_x, p3rec_y, p3rec_z)

    theta_pq    = angle_deg(pL[:, 0], pL[:, 1], pL[:, 2], qx, qy, qz)
    lead_over_q = pL_mag / (q_mag + 1e-30)

    mask = (np.isfinite(w) & (w > 0) & detector
            & (theta_e > theta_e_min) & (theta_e < theta_e_max)
            & (Q2 > Q2_min)
            & (xB > xB_min) & (xB < xB_max)
            & (theta_pq < theta_pq_max)
            & (lead_over_q > lead_over_q_min)
            & (lead_over_q < lead_over_q_max)
            & (p_miss_mag > pmiss_lo) & (p_miss_mag < pmiss_hi))
    print(f'  total entries: {len(w)}; passing 3N cuts ({mode}): {int(mask.sum())}')

    if stats is not None:
        # Base cuts contain NO xB requirement: theta_e window, Q2, and the
        # detector tag only. Any xB window belongs to the chosen SRC cut set.
        ok = np.isfinite(w) & (w > 0)
        base = (ok & detector
                & (theta_e > theta_e_min) & (theta_e < theta_e_max)
                & (Q2 > Q2_min))
        if src_from_args:
            src = (base & (theta_pq < theta_pq_max)
                   & (lead_over_q > lead_over_q_min)
                   & (lead_over_q < lead_over_q_max)
                   & (p_miss_mag > pmiss_lo) & (p_miss_mag < pmiss_hi)
                   & (xB > xB_min) & (xB < xB_max))
        else:
            src = (base & (theta_pq < SRC_THETA_PQ_MAX)
                   & (lead_over_q > SRC_LEAD_OVER_Q_MIN)
                   & (p_miss_mag > SRC_PMISS_LO) & (p_miss_mag < SRC_PMISS_HI))
        stats.update(w_all=float(w[ok].sum()), n_all=int(ok.sum()),
                     w_base=float(w[base].sum()), n_base=int(base.sum()),
                     w_src=float(w[src].sum()), n_src=int(src.sum()))

    return theta12[mask], theta23[mask], w[mask]


# ---------- 2N sample (processed events_2N.root) ----------------------------
def load_2N(path, theta_e_max, Q2_min, xB_max, kF_thr=kF,
            theta_pq_max=180.0, lead_over_q_min=0.0,
            pmiss_lo=0.0, pmiss_hi=1e9, mode='ge2',
            theta_e_min=0.0, stats=None,
            lead_over_q_max=1e9, xB_min=0.0, src_from_args=False):
    print(f'Reading 2N+FSI: {path}')
    with uproot.open(path) as f:
        ev = f['events']
        arr = ev.arrays(
            ['weight', 'lead_post', 'recoil_post', 'q',
             'Q2', 'xB', 'scattering_angle', 'pmiss', 'nAboveKF'],
            library='np',
        )
        sec = ev.arrays(['sec_px', 'sec_py', 'sec_pz', 'sec_pdg'], library='ak')
    w        = arr['weight']
    Q2       = arr['Q2']
    xB       = arr['xB']
    theta_e  = arr['scattering_angle']
    pmiss    = arr['pmiss']
    nAbove   = arr['nAboveKF']
    lp = arr['lead_post']; rp = arr['recoil_post']; q = arr['q']

    lp_mag = np.sqrt(lp[:, 0]**2 + lp[:, 1]**2 + lp[:, 2]**2)
    rp_mag = np.sqrt(rp[:, 0]**2 + rp[:, 1]**2 + rp[:, 2]**2)
    q_mag  = np.sqrt(q[:, 0]**2 + q[:, 1]**2 + q[:, 2]**2)

    # Detection-based recoil: the two nucleons whose opening angles we form are
    # the lead (fixed by the detector cuts, consistent with q) and the UNIQUE
    # nucleon above kF. That is recoil_post when it survived FSI above kF;
    # otherwise (recoil knocked below kF, an ejected secondary took its place in
    # the nAboveKF==2 count) it is the leading above-kF secondary nucleon. Using
    # the below-kF recoil_post for the angle would be using an undetected nucleon.
    spx, spy, spz, spdg = sec['sec_px'], sec['sec_py'], sec['sec_pz'], sec['sec_pdg']
    sp_mag = np.sqrt(spx**2 + spy**2 + spz**2)
    sec_above = ((spdg == 2212) | (spdg == 2112)) & (sp_mag > kF_thr)
    sp_key = ak.where(sec_above, sp_mag, -1.0)
    has_sec = ak.to_numpy(ak.fill_none(ak.max(sp_key, axis=1) > 0, False))
    sidx = ak.argmax(sp_key, axis=1, keepdims=True)
    sec_x = ak.to_numpy(ak.fill_none(ak.flatten(spx[sidx]), 0.0))
    sec_y = ak.to_numpy(ak.fill_none(ak.flatten(spy[sidx]), 0.0))
    sec_z = ak.to_numpy(ak.fill_none(ak.flatten(spz[sidx]), 0.0))
    use_recoil = rp_mag > kF_thr
    rdx = np.where(use_recoil, rp[:, 0], sec_x)
    rdy = np.where(use_recoil, rp[:, 1], sec_y)
    rdz = np.where(use_recoil, rp[:, 2], sec_z)
    have_recoil = use_recoil | has_sec   # a detected above-kF recoil exists

    if mode == 'eq2':
        # Strict tag: tree's nAboveKF == 2, with the detected above-kF recoil
        # well-defined. Mirrors the paper's 2N+FSI cut.
        detector = (nAbove == 2) & have_recoil
    elif mode == 'ge2':
        # Loose tag: lead and recoil both above kF (FSI may have kicked
        # extra nucleons above kF; we don't have their 4-vectors anyway).
        detector = (lp_mag > kF_thr) & (rp_mag > kF_thr)
    else:
        raise ValueError(f"mode must be 'eq2' or 'ge2', got {mode!r}")

    p_miss_x = lp[:, 0] - q[:, 0]
    p_miss_y = lp[:, 1] - q[:, 1]
    p_miss_z = lp[:, 2] - q[:, 2]

    p3rec_x = -(p_miss_x + rdx)
    p3rec_y = -(p_miss_y + rdy)
    p3rec_z = -(p_miss_z + rdz)

    theta12 = angle_deg(p_miss_x, p_miss_y, p_miss_z,
                        rdx, rdy, rdz)
    theta23 = angle_deg(rdx, rdy, rdz,
                        p3rec_x, p3rec_y, p3rec_z)

    theta_pq    = angle_deg(lp[:, 0], lp[:, 1], lp[:, 2],
                            q[:, 0], q[:, 1], q[:, 2])
    lead_over_q = lp_mag / (q_mag + 1e-30)

    mask = (np.isfinite(w) & (w > 0)
            & detector
            & (theta_e > theta_e_min) & (theta_e < theta_e_max)
            & (Q2 > Q2_min)
            & (xB > xB_min) & (xB < xB_max)
            & (theta_pq < theta_pq_max)
            & (lead_over_q > lead_over_q_min)
            & (lead_over_q < lead_over_q_max)
            & (pmiss > pmiss_lo) & (pmiss < pmiss_hi))
    print(f'  total entries: {len(w)}; passing 2N cuts ({mode}): {int(mask.sum())}')

    if stats is not None:
        # Base cuts contain NO xB requirement (see load_3N).
        ok = np.isfinite(w) & (w > 0)
        base = (ok & detector
                & (theta_e > theta_e_min) & (theta_e < theta_e_max)
                & (Q2 > Q2_min))
        if src_from_args:
            src = (base & (theta_pq < theta_pq_max)
                   & (lead_over_q > lead_over_q_min)
                   & (lead_over_q < lead_over_q_max)
                   & (pmiss > pmiss_lo) & (pmiss < pmiss_hi)
                   & (xB > xB_min) & (xB < xB_max))
        else:
            src = (base & (theta_pq < SRC_THETA_PQ_MAX)
                   & (lead_over_q > SRC_LEAD_OVER_Q_MIN)
                   & (pmiss > SRC_PMISS_LO) & (pmiss < SRC_PMISS_HI))
        stats.update(w_all=float(w[ok].sum()), n_all=int(ok.sum()),
                     w_base=float(w[base].sum()), n_base=int(base.sum()),
                     w_src=float(w[src].sum()), n_src=int(src.sum()))

    return theta12[mask], theta23[mask], w[mask]


def plot_panel(ax, t12, t23, w, bins, vmin, vmax):
    h, xe, ye = np.histogram2d(
        t12, t23, bins=bins, range=[[0, 180], [0, 180]], weights=w,
    )
    h_norm = h / (h.sum() + 1e-30)
    h_plot = np.where(h_norm > 0, h_norm, np.nan)
    im = ax.pcolormesh(xe, ye, h_plot.T, cmap=CMAP,
                       norm=LogNorm(vmin=vmin, vmax=vmax),
                       shading='auto', rasterized=True)

    x = np.linspace(0, 180, 200)
    line_kw = dict(color=OVERLAY_COLOR, linestyle='--', lw=0.7, alpha=0.85)
    ax.plot(x, 180.0 - x / 2.0, **line_kw)
    ax.plot(180.0 - x / 2.0, x, **line_kw)
    ax.plot(x, x, **line_kw)

    ax.set_xlim(0, 180); ax.set_ylim(0, 180)
    ax.set_aspect('equal')
    ax.set_xticks(np.arange(0, 181, 30))
    ax.set_yticks(np.arange(0, 181, 30))
    ax.set_xlabel(r'$\theta_{12}$ (deg)')
    return im


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--input-3N', default='events/3N_FSI_5M_12C.root')
    p.add_argument('--input-2N',
                   default='../genQE_FSI/events/misc/events_2N.root')
    p.add_argument('--output',
                   default='analysis/Plots/theta_heatmap_3N_vs_2N.pdf')
    p.add_argument('--ebeam', type=float, default=6.0)
    p.add_argument('--bins',  type=int,   default=36)
    p.add_argument('--dpi',   type=int,   default=400)
    p.add_argument('--theta-e-min', type=float, default=8.0)
    p.add_argument('--theta-e-max', type=float, default=45.0)
    p.add_argument('--Q2-min',      type=float, default=1.0)
    p.add_argument('--xB-max',      type=float, default=1.2)
    p.add_argument('--kF',          type=float, default=kF)
    p.add_argument('--vmin', type=float, default=None,
                   help='default: vmax * vmin-frac')
    p.add_argument('--vmin-frac', type=float, default=1e-4,
                   help='vmin as a fraction of each panel vmax '
                        '(e.g. 1e-8 for a deep high-res color floor)')
    p.add_argument('--vmax', type=float, default=None,
                   help='default: largest normalized bin over both panels')
    p.add_argument('--mode', choices=('eq2', 'ge2'), default='ge2',
                   help='detector tag: eq2 = exactly 2 above kF, ge2 = at least 2')
    p.add_argument('--theta-pq-max',    type=float, default=180.0)
    p.add_argument('--lead-over-q-min', type=float, default=0.0)
    p.add_argument('--lead-over-q-max', type=float, default=1e9)
    p.add_argument('--pmiss-lo',        type=float, default=0.0)
    p.add_argument('--pmiss-hi',        type=float, default=1e9)
    p.add_argument('--xB-min',          type=float, default=0.0)
    p.add_argument('--src-stats-from-args', action='store_true',
                   help='report the SRC pass fraction using the cut values '
                        'given on the command line instead of the built-in '
                        'SRC_* constants')
    p.add_argument('--triangles', action='store_true',
                   help='outline the L (top-left) and R (top-right) regions')
    p.add_argument('--label-3N', default=r'3N+FSI on $^{12}$C')
    p.add_argument('--label-2N', default=r'2N+FSI')
    args = p.parse_args()

    apply_style()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    root_dir   = os.path.abspath(os.path.join(script_dir, '..'))
    in3 = args.input_3N if os.path.isabs(args.input_3N) else os.path.join(root_dir, args.input_3N)
    in2 = args.input_2N if os.path.isabs(args.input_2N) else os.path.join(root_dir, args.input_2N)
    out_path = args.output if os.path.isabs(args.output) else os.path.join(root_dir, args.output)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    stats_3, stats_2 = {}, {}
    t12_3, t23_3, w_3 = load_3N(in3, args.ebeam,
                                args.theta_e_max, args.Q2_min, args.xB_max,
                                args.kF,
                                theta_pq_max=args.theta_pq_max,
                                lead_over_q_min=args.lead_over_q_min,
                                lead_over_q_max=args.lead_over_q_max,
                                pmiss_lo=args.pmiss_lo,
                                pmiss_hi=args.pmiss_hi,
                                xB_min=args.xB_min,
                                mode=args.mode,
                                theta_e_min=args.theta_e_min,
                                src_from_args=args.src_stats_from_args,
                                stats=stats_3)
    t12_2, t23_2, w_2 = load_2N(in2,
                                args.theta_e_max, args.Q2_min, args.xB_max,
                                args.kF,
                                theta_pq_max=args.theta_pq_max,
                                lead_over_q_min=args.lead_over_q_min,
                                lead_over_q_max=args.lead_over_q_max,
                                pmiss_lo=args.pmiss_lo,
                                pmiss_hi=args.pmiss_hi,
                                xB_min=args.xB_min,
                                mode=args.mode,
                                theta_e_min=args.theta_e_min,
                                src_from_args=args.src_stats_from_args,
                                stats=stats_2)

    n3 = len(w_3); n2 = len(w_2)

    # ---- L/R region fractions + base/SRC pass-fraction report ----
    from plot_fig5_theta_heatmap_pair import (in_region_L_array,
                                              in_region_R_array,
                                              triangle_vertices_L,
                                              triangle_vertices_R,
                                              _add_triangle)
    summary_lines = [f'# cuts: mode={args.mode}, kF={args.kF}, '
                     f'{args.theta_e_min}<theta_e<{args.theta_e_max}, '
                     f'Q2>{args.Q2_min}, {args.xB_min}<xB<{args.xB_max}, '
                     f'theta_pq<{args.theta_pq_max}, '
                     f'{args.lead_over_q_min}<lead/q<{args.lead_over_q_max}, '
                     f'{args.pmiss_lo}<pmiss<{args.pmiss_hi}']
    for name, t12, t23, w, st in [('3N', t12_3, t23_3, w_3, stats_3),
                                  ('2N', t12_2, t23_2, w_2, stats_2)]:
        wt = w.sum() + 1e-30
        fL = w[in_region_L_array(t12, t23)].sum() / wt
        fR = w[in_region_R_array(t12, t23)].sum() / wt
        ratio = (fL / fR) if fR > 0 else float('inf')
        f_base = st['w_base'] / (st['w_all'] + 1e-30)
        f_src = st['w_src'] / (st['w_all'] + 1e-30)
        for line in (
            f'[{name}] f_L={fL:.4f}  f_R={fR:.4f}  L/R={ratio:.3f}',
            f'[{name}] base-cut pass fraction (weighted): {f_base:.4f} '
            f'({st["n_base"]}/{st["n_all"]} events)',
            f'[{name}] SRC-cut  pass fraction (weighted): {f_src:.5f} '
            f'({st["n_src"]}/{st["n_all"]} events)',
        ):
            print('  ' + line)
            summary_lines.append(line)
    summary_path = os.path.splitext(out_path)[0] + '_summary.txt'

    def _hist_max(t12, t23, w):
        h, _, _ = np.histogram2d(t12, t23, bins=args.bins,
                                 range=[[0, 180], [0, 180]], weights=w)
        h = h / (h.sum() + 1e-30)
        pos = h[h > 0]
        return pos.max() if pos.size else 1e-2

    # Independent color scale per panel (the two samples are separately
    # normalized, so a shared scale would be misleading).
    vmax_3 = args.vmax if args.vmax is not None else _hist_max(t12_3, t23_3, w_3)
    vmax_2 = args.vmax if args.vmax is not None else _hist_max(t12_2, t23_2, w_2)
    vmin_3 = args.vmin if args.vmin is not None else vmax_3 * args.vmin_frac
    vmin_2 = args.vmin if args.vmin is not None else vmax_2 * args.vmin_frac

    fig, axes = plt.subplots(1, 2, sharey=True,
                             figsize=figure_size(cols=2, ratio=0.50))
    im_left  = plot_panel(axes[0], t12_3, t23_3, w_3,
                          args.bins, vmin_3, vmax_3)
    im_right = plot_panel(axes[1], t12_2, t23_2, w_2,
                          args.bins, vmin_2, vmax_2)

    axes[0].set_ylabel(r'$\theta_{23}$ (deg)')
    axes[1].tick_params(labelleft=False)

    tag_lbl = r'$n_{>k_F}=2$' if args.mode == 'eq2' else r'$n_{>k_F}\geq 2$'
    axes[0].set_title(r'%s (%s, $N=%d$)' % (args.label_3N, tag_lbl, n3),
                      fontsize=8)
    axes[1].set_title(r'%s (%s, $N=%d$)' % (args.label_2N, tag_lbl, n2),
                      fontsize=8)

    if args.triangles:
        tri_L = triangle_vertices_L()
        tri_R = triangle_vertices_R()
        for ax in axes:
            _add_triangle(ax, tri_L, edgecolor='k', lw=1.0, label_text='L')
            _add_triangle(ax, tri_R, edgecolor='k', lw=1.0, label_text='R')

    cb1 = fig.colorbar(im_left, ax=axes[0], pad=0.02, fraction=0.046)
    cb2 = fig.colorbar(im_right, ax=axes[1], label='Normalized weight',
                       pad=0.02, fraction=0.046)
    for cb in (cb1, cb2):
        cb.ax.tick_params(labelsize=8)
    cb2.ax.yaxis.label.set_size(9)

    fig.savefig(out_path, dpi=args.dpi, bbox_inches='tight')
    print(f'Saved {out_path}')

    with open(summary_path, 'w') as fh:
        fh.write('\n'.join(summary_lines) + '\n')
    print(f'Saved {summary_path}')


if __name__ == '__main__':
    main()
