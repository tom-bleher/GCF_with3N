#!/usr/bin/env python3
"""Single-panel TRUE opening-angle heatmap from a raw genQE_3N(_FSI) tree
(intended for the no-FSI / PWIA sample), with four signal regions integrated
(numbers go to the _summary.txt only; no overlays are drawn on the plot):

    L      top-left triangle   (A=135 deg, K=4; fig5 geometry)
    R      top-right triangle
    BR     bottom-right triangle (R reflected about theta23 -> 360-t12-t23)
    Star   circle at (120, 120) deg, equal in area to the triangles

TRUE angles of the three generated nucleons — no detector tag, no
reconstruction (for a -K sample every initial nucleon already has k > kF):

    theta12 = angle(p1, p2)  with p1 = pLead - q (initial struck momentum)
    theta23 = angle(p2, p3)

Cuts: weight > 0, theta_e in (theta-e-min, theta-e-max), Q2 > Q2-min,
xB < xB-max.  The summary also reports the weighted base-cut and SRC-cut
(theta_pq < 8, |p_lead|/|q| > 0.75, 0.25 < pmiss < 0.9) pass fractions.
"""
import argparse
import os
import sys

import numpy as np
import uproot
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paper_style import apply_style, figure_size, CMAP, OVERLAY_COLOR
from plot_fig5_theta_heatmap_pair import (in_region_L_array, in_region_R_array,
                                          triangle_vertices_R)

mN = 0.93892

SRC_THETA_PQ_MAX = 8.0
SRC_LEAD_OVER_Q_MIN = 0.75
SRC_PMISS_LO = 0.25
SRC_PMISS_HI = 0.90


def angle_deg(ax, ay, az, bx, by, bz):
    dot = ax * bx + ay * by + az * bz
    a = np.sqrt(ax * ax + ay * ay + az * az)
    b = np.sqrt(bx * bx + by * by + bz * bz)
    return np.degrees(np.arccos(np.clip(dot / (a * b + 1e-30), -1.0, 1.0)))


def in_region_BR_array(t12, t23):
    return in_region_R_array(t12, 360.0 - t12 - t23)


def triangle_area(pts):
    (x1, y1), (x2, y2), (x3, y3) = pts
    return 0.5 * abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1))


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--input', default='events/3N_PWIA_20M_12C_kF220_6GeV.root')
    p.add_argument('--output',
                   default='analysis/Plots/theta_heatmap_3N_noFSI_regions.pdf')
    p.add_argument('--title', default=r'3N (no FSI) on $^{12}$C')
    p.add_argument('--ebeam', type=float, default=6.0)
    p.add_argument('--bins',  type=int,   default=72)
    p.add_argument('--dpi',   type=int,   default=400)
    p.add_argument('--theta-e-min', type=float, default=8.0)
    p.add_argument('--theta-e-max', type=float, default=45.0)
    p.add_argument('--Q2-min',      type=float, default=1.0)
    p.add_argument('--xB-max',      type=float, default=1.2)
    p.add_argument('--vmin', type=float, default=None)
    p.add_argument('--vmin-frac', type=float, default=1e-4)
    p.add_argument('--vmax', type=float, default=None)
    args = p.parse_args()

    apply_style()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    root_dir = os.path.abspath(os.path.join(script_dir, '..'))
    in_path = args.input if os.path.isabs(args.input) else os.path.join(root_dir, args.input)
    out_path = args.output if os.path.isabs(args.output) else os.path.join(root_dir, args.output)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    print(f'Reading 3N: {in_path}')
    with uproot.open(in_path) as f:
        tree = f['genT']
        # Branch names: the current generator writes p1_gen/p2_gen/p3_gen,
        # older trees wrote pLead/p2/p3 (and pLead_pre/p2_pre/p3_pre pre-FSI).
        have = set(tree.keys())
        if 'p1_gen' in have:
            gen = ('p1_gen', 'p2_gen', 'p3_gen')
        elif 'pLead' in have:
            gen = ('pLead', 'p2', 'p3')
        else:
            gen = ('pLead_pre', 'p2_pre', 'p3_pre')
        arr = tree.arrays(['weight', 'pe', *gen], library='np')
    w = arr['weight']
    pe = arr['pe']; pL = arr[gen[0]]; p2 = arr[gen[1]]; p3 = arr[gen[2]]

    qx = -pe[:, 0]; qy = -pe[:, 1]; qz = args.ebeam - pe[:, 2]
    q_mag = np.sqrt(qx**2 + qy**2 + qz**2)
    pe_mag = np.sqrt((pe**2).sum(axis=1))
    omega = args.ebeam - pe_mag
    Q2 = q_mag**2 - omega**2
    xB = Q2 / (2.0 * mN * omega + 1e-30)
    theta_e = np.degrees(np.arctan2(np.sqrt(pe[:, 0]**2 + pe[:, 1]**2),
                                    pe[:, 2]))

    # TRUE initial-state angles: p1 = pLead - q; recoils as generated.
    p1x = pL[:, 0] - qx; p1y = pL[:, 1] - qy; p1z = pL[:, 2] - qz
    theta12_all = angle_deg(p1x, p1y, p1z, p2[:, 0], p2[:, 1], p2[:, 2])
    theta23_all = angle_deg(p2[:, 0], p2[:, 1], p2[:, 2],
                            p3[:, 0], p3[:, 1], p3[:, 2])

    ok = np.isfinite(w) & (w > 0)
    base = (ok & (theta_e > args.theta_e_min) & (theta_e < args.theta_e_max)
            & (Q2 > args.Q2_min) & (xB < args.xB_max))

    pmiss = np.sqrt(p1x**2 + p1y**2 + p1z**2)
    pL_mag = np.sqrt((pL**2).sum(axis=1))
    theta_pq = angle_deg(pL[:, 0], pL[:, 1], pL[:, 2], qx, qy, qz)
    src = (base & (theta_pq < SRC_THETA_PQ_MAX)
           & (pL_mag / (q_mag + 1e-30) > SRC_LEAD_OVER_Q_MIN)
           & (pmiss > SRC_PMISS_LO) & (pmiss < SRC_PMISS_HI))

    t12, t23, ww = theta12_all[base], theta23_all[base], w[base]
    print(f'  total entries: {len(w)}; passing base cuts: {int(base.sum())}')

    # ---- region integrals ----
    star_radius = np.sqrt(triangle_area(triangle_vertices_R()) / np.pi)
    in_star = ((t12 - 120.0) ** 2 + (t23 - 120.0) ** 2) < star_radius ** 2

    wt = ww.sum() + 1e-30
    regions = [
        ('L (top-left)',      in_region_L_array(t12, t23)),
        ('R (top-right)',     in_region_R_array(t12, t23)),
        ('BR (bottom-right)', in_region_BR_array(t12, t23)),
        (f'Star (120,120), r={star_radius:.1f} deg', in_star),
    ]
    fracs = {name: float(ww[m].sum()) / wt for name, m in regions}
    fL = fracs['L (top-left)']
    fR = fracs['R (top-right)']
    ratio = fL / fR if fR > 0 else float('inf')
    f_base = float(w[base].sum()) / (float(w[ok].sum()) + 1e-30)
    f_src = float(w[src].sum()) / (float(w[ok].sum()) + 1e-30)

    summary_lines = [
        f'# input: {in_path}',
        f'# TRUE 3-nucleon angles (p1 = pLead - q); no detector tag, no '
        f'reconstruction',
        f'# cuts: {args.theta_e_min}<theta_e<{args.theta_e_max}, '
        f'Q2>{args.Q2_min}, xB<{args.xB_max}',
    ]
    for name, _ in regions:
        summary_lines.append(f'{name}: {100.0 * fracs[name]:.2f}% of plotted weight')
    summary_lines += [
        f'L/R = {ratio:.3f}',
        f'base-cut pass fraction (weighted): {f_base:.4f} '
        f'({int(base.sum())}/{int(ok.sum())} events)',
        f'SRC-cut  pass fraction (weighted): {f_src:.5f} '
        f'({int(src.sum())}/{int(ok.sum())} events)',
    ]
    print('\n'.join('  ' + s for s in summary_lines))

    # ---- plot ----
    h, xe, ye = np.histogram2d(t12, t23, bins=args.bins,
                               range=[[0, 180], [0, 180]], weights=ww)
    h_norm = h / (h.sum() + 1e-30)
    h_plot = np.where(h_norm > 0, h_norm, np.nan)
    pos = h_norm[h_norm > 0]
    vmax = args.vmax if args.vmax is not None else (pos.max() if pos.size else 1e-2)
    vmin = args.vmin if args.vmin is not None else vmax * args.vmin_frac

    fig, ax = plt.subplots(figsize=figure_size(cols=1, ratio=0.92))
    im = ax.pcolormesh(xe, ye, h_plot.T, cmap=CMAP,
                       norm=LogNorm(vmin=vmin, vmax=vmax),
                       shading='auto', rasterized=True)
    cbar = fig.colorbar(im, ax=ax, label='Normalized weight',
                        pad=0.02, fraction=0.046)
    cbar.ax.tick_params(labelsize=8)
    cbar.ax.yaxis.label.set_size(9)

    x = np.linspace(0, 180, 200)
    line_kw = dict(color=OVERLAY_COLOR, linestyle='--', lw=0.7, alpha=0.85)
    ax.plot(x, 180.0 - x / 2.0, **line_kw)
    ax.plot(180.0 - x / 2.0, x, **line_kw)
    ax.plot(x, x, **line_kw)

    ax.set_xlim(0, 180)
    ax.set_ylim(0, 180)
    ax.set_aspect('equal')
    ax.set_xticks(np.arange(0, 181, 30))
    ax.set_yticks(np.arange(0, 181, 30))
    ax.set_xlabel(r'$\theta_{12}$ (deg)')
    ax.set_ylabel(r'$\theta_{23}$ (deg)')
    ax.set_title(r'%s ($N=%d$)' % (args.title, len(ww)), fontsize=8)

    fig.tight_layout()
    fig.savefig(out_path, dpi=args.dpi)
    print(f'Saved {out_path}')

    summary_path = os.path.splitext(out_path)[0] + '_summary.txt'
    with open(summary_path, 'w') as fh:
        fh.write('\n'.join(summary_lines) + '\n')
    print(f'Saved {summary_path}')


if __name__ == '__main__':
    main()
