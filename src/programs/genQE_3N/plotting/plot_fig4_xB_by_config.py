#!/usr/bin/env python3
"""Fig. 4 of the 3N-SRC paper.

xB distribution for the three angular configurations of the 3N-SRC, defined by
the line-bounded regions in (theta12, theta23) (same as Fig. 2):

    Lead Rocket     region 0  (head, near (180, 0))
    Recoil Rocket   region 2  (tail, near (180, 180))
    Star            region 6  (star, around (120, 120))

Cuts (paper convention): theta_e<45 deg, Q^2>1 GeV^2, no CM momentum.
Same 99.9-percentile weight clipping as Fig. 2 to tame heavy-tailed weights.
"""
import argparse
import os
import sys
import numpy as np
import uproot
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paper_style import apply_style, figure_size, LINE_COLORS
from regions import in_region


# (label, region_idx, color) — line-bounded regions from
# analysis/SRC_analysis_3N.cpp (desmos link in regions.py).
REGIONS = [
    ('Lead Rocket',   0, LINE_COLORS[0]),   # head rocket (near (180, 0))
    ('Recoil Rocket', 2, LINE_COLORS[1]),   # tail rocket (near (180, 180))
    ('Star',          6, LINE_COLORS[2]),   # star (around (120, 120))
]


def angle_deg(ax, ay, az, bx, by, bz):
    dot = ax * bx + ay * by + az * bz
    a = np.sqrt(ax * ax + ay * ay + az * az)
    b = np.sqrt(bx * bx + by * by + bz * bz)
    return np.degrees(np.arccos(np.clip(dot / (a * b + 1e-30), -1.0, 1.0)))


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--input', default='events_3N.root')
    p.add_argument('--output', default='analysis_output/png_files/Paper plots/fig4_xB_by_config.pdf')
    p.add_argument('--bins', type=int, default=50)
    p.add_argument('--xB-min', type=float, default=0.0)
    p.add_argument('--xB-max', type=float, default=3.0)
    p.add_argument('--theta-e-max', type=float, default=45.0)
    p.add_argument('--Q2-min', type=float, default=1.0)
    p.add_argument('--weight-clip-pct', type=float, default=99.9)
    p.add_argument('--dpi', type=int, default=400)
    p.add_argument('--logy', action='store_true',
                   help='log y axis (the high-xB tail is ~1e-3 of the peak and '
                        'is invisible on a linear scale)')
    args = p.parse_args()

    apply_style()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    root_dir = os.path.abspath(os.path.join(script_dir, '..'))
    in_path = args.input if os.path.isabs(args.input) else os.path.join(root_dir, args.input)
    out_path = args.output if os.path.isabs(args.output) else os.path.join(root_dir, args.output)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    print(f'Reading {in_path}')
    tree = uproot.open(in_path)['events']
    arr = tree.arrays(
        ['weight', 'lead', 'recoil2', 'recoil3', 'q', 'Q2', 'xB',
         'scattering_angle'],
        library='np',
    )
    print(f'  {len(arr["weight"])} events')

    w = arr['weight']
    Q2 = arr['Q2']
    xB = arr['xB']
    theta_e = arr['scattering_angle']
    lead, r2, r3, q = arr['lead'], arr['recoil2'], arr['recoil3'], arr['q']

    p1x = lead[:, 0] - q[:, 0]
    p1y = lead[:, 1] - q[:, 1]
    p1z = lead[:, 2] - q[:, 2]
    theta12 = angle_deg(p1x, p1y, p1z, r2[:, 0], r2[:, 1], r2[:, 2])
    theta23 = angle_deg(r2[:, 0], r2[:, 1], r2[:, 2],
                        r3[:, 0], r3[:, 1], r3[:, 2])

    base = (np.isfinite(w) & (w > 0)
            & (theta_e < args.theta_e_max)
            & (Q2 > args.Q2_min))
    print(f'  base events after cuts: {int(base.sum())}')

    if args.weight_clip_pct < 100.0:
        clip = float(np.percentile(w[base], args.weight_clip_pct))
        w = np.minimum(w, clip)
        print(f'  weight clip at {args.weight_clip_pct}% = {clip:.3e}')

    edges = np.linspace(args.xB_min, args.xB_max, args.bins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    bin_width = edges[1] - edges[0]

    fig, ax = plt.subplots(figsize=figure_size(cols=2, ratio=0.38))

    w_base_sum = float(w[base].sum())
    print(f'  base weight sum (post-clip) = {w_base_sum:.6e}')

    for label, ireg, color in REGIONS:
        m = base & in_region(theta12, theta23, ireg)
        n_in = int(m.sum())
        w_in = float(w[m].sum())
        frac = w_in / w_base_sum if w_base_sum > 0 else 0.0
        print(f'  [{label}] region={ireg} : {n_in} events, '
              f'sum(w)={w_in:.6e}  ({100*frac:.2f}% of base)')

        counts, _ = np.histogram(xB[m], bins=edges, weights=w[m])
        ax.plot(centers, counts, color=color, lw=1.4, label=label)

    ax.set_xlim(args.xB_min, args.xB_max)
    if args.logy:
        ax.set_yscale('log')
    else:
        ax.set_ylim(bottom=0)
    ax.set_xlabel(r'$x_B$')
    ax.set_ylabel(r'Total weight')
    ax.legend(loc='upper right', frameon=False, handlelength=1.6)

    fig.tight_layout()
    fig.savefig(out_path, dpi=args.dpi)
    png_path = os.path.splitext(out_path)[0] + '.png'
    fig.savefig(png_path, dpi=args.dpi)
    print(f'Saved {out_path}')
    print(f'Saved {png_path} (preview)')


if __name__ == '__main__':
    main()
