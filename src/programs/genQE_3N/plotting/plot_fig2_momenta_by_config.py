#!/usr/bin/env python3
"""Fig. 2 of the 3N-SRC paper.

Three subplots showing the ground-state momentum distributions of the three
nucleons in three angular configurations of the 3N-SRC, defined as the
line-bounded regions of regions.py in the (theta12, theta23) plane:

    Lead Rocket     region 0, near (180, 0)    -> lead is the rocket nucleon
    Recoil Rocket   region 2, near (180, 180)  -> recoil2 is the rocket nucleon
    Star            region 6, near (120, 120)  -> all three at ~120 deg apart

Cuts (paper convention): theta_e<45 deg, Q^2>1 GeV^2, no CM momentum.
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


# (label, region_idx) — using the line-bounded regions from
# analysis/SRC_analysis_3N.cpp (desmos link in regions.py).
REGIONS = [
    ('Lead Rocket',   0),   # head rocket (near (180, 0))
    ('Recoil Rocket', 2),   # tail rocket (near (180, 180))
    ('Star',          6),   # star (around (120, 120))
]


def angle_deg(ax, ay, az, bx, by, bz):
    dot = ax * bx + ay * by + az * bz
    a = np.sqrt(ax * ax + ay * ay + az * az)
    b = np.sqrt(bx * bx + by * by + bz * bz)
    return np.degrees(np.arccos(np.clip(dot / (a * b + 1e-30), -1.0, 1.0)))


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--input', default='events_3N.root')
    p.add_argument('--output', default='analysis_output/png_files/Paper plots/fig2_momenta_by_config.pdf')
    p.add_argument('--bins', type=int, default=20)
    p.add_argument('--p-min', type=float, default=0.3, help='Momentum range minimum [GeV/c]')
    p.add_argument('--p-max', type=float, default=1.0, help='Momentum range maximum [GeV/c]')
    p.add_argument('--weight-clip-pct', type=float, default=99.9,
                   help='Clip weights at this percentile (computed on the post-cut sample)')
    p.add_argument('--theta-e-max', type=float, default=45.0)
    p.add_argument('--Q2-min', type=float, default=1.0)
    p.add_argument('--dpi', type=int, default=400)
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
        ['weight', 'lead', 'recoil2', 'recoil3', 'q', 'Q2', 'scattering_angle'],
        library='np',
    )
    n_total = len(arr['weight'])
    print(f'  {n_total} events')

    w = arr['weight']
    Q2 = arr['Q2']
    theta_e = arr['scattering_angle']
    lead, r2, r3, q = arr['lead'], arr['recoil2'], arr['recoil3'], arr['q']

    # Initial lead momentum = lead (final) - q
    p1x = lead[:, 0] - q[:, 0]
    p1y = lead[:, 1] - q[:, 1]
    p1z = lead[:, 2] - q[:, 2]
    p1_mag = np.sqrt(p1x**2 + p1y**2 + p1z**2)
    r2_mag = np.sqrt(r2[:, 0]**2 + r2[:, 1]**2 + r2[:, 2]**2)
    r3_mag = np.sqrt(r3[:, 0]**2 + r3[:, 1]**2 + r3[:, 2]**2)

    theta12 = angle_deg(p1x, p1y, p1z, r2[:, 0], r2[:, 1], r2[:, 2])
    theta23 = angle_deg(r2[:, 0], r2[:, 1], r2[:, 2],
                        r3[:, 0], r3[:, 1], r3[:, 2])

    base = (np.isfinite(w) & (w > 0)
            & (theta_e < args.theta_e_max)
            & (Q2 > args.Q2_min))
    print(f'  base events after cuts: {int(base.sum())}')

    # Clip weights at a high percentile (computed on the post-cut sample) to
    # tame extreme outliers that dominate weighted histograms despite millions
    # of events. This trades a small bias in the heavy tail for stability.
    if args.weight_clip_pct < 100.0:
        clip = float(np.percentile(w[base], args.weight_clip_pct))
        n_clipped = int((w[base] > clip).sum())
        w = np.minimum(w, clip)
        print(f'  weight clip at {args.weight_clip_pct}% = {clip:.3e} '
              f'({n_clipped} events affected)')

    # The recoil-rocket panel uses region 2, near (180, 180): theta12 ~ 180
    # and theta23 ~ 180, so p1 || r3 and r2 is opposite both, i.e. r2 is the
    # rocket. (Region 3, near (0, 180), is the mirror image with r3 as the
    # rocket.) Map for plotting so the "Recoil 1" curve is the rocket recoil:
    #   Lead     -> p1
    #   Recoil 1 -> r2 (the high-momentum recoil at (180, 180))
    #   Recoil 2 -> r3
    nucleons = [
        ('Lead',     p1_mag, LINE_COLORS[0], '-'),
        ('Recoil 1', r2_mag, LINE_COLORS[1], '-'),
        ('Recoil 2', r3_mag, LINE_COLORS[2], '-'),
    ]

    fig, axes = plt.subplots(1, 3, figsize=figure_size(cols=2, ratio=0.42),
                             sharey=False)

    edges = np.linspace(args.p_min, args.p_max, args.bins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    bin_width = edges[1] - edges[0]

    w_base_sum = float(w[base].sum())
    print(f'  base weight sum (post-clip) = {w_base_sum:.6e}')

    # All seven line-bounded regions (0-6) — main: 0,2,3,6 ; transitions: 1,4,5
    all_region_names = {
        0: 'Lead Rocket (head, 180/0)',
        1: 'between head & star',
        2: 'Recoil Rocket (tail, 180/180)',
        3: 'Recoil Rocket (tail, 0/180)',
        4: 'between (180,180) & star',
        5: 'between (0,180) & star',
        6: 'Star (120/120)',
    }
    print('  --- all regions ---')
    total_check = 0.0
    for ireg, rname in all_region_names.items():
        m_all = base & in_region(theta12, theta23, ireg)
        n_all = int(m_all.sum())
        w_all = float(w[m_all].sum())
        frac_all = w_all / w_base_sum if w_base_sum > 0 else 0.0
        total_check += frac_all
        print(f'  region={ireg} {rname:<32s}: {n_all:>9d} events, '
              f'sum(w)={w_all:.6e}  ({100*frac_all:.2f}% of base)')
    print(f'  total across regions 0-6: {100*total_check:.2f}% of base')
    print('  --- panel regions ---')

    for ax, (label, ireg) in zip(axes, REGIONS):
        m = base & in_region(theta12, theta23, ireg)
        n_in = int(m.sum())
        w_in = float(w[m].sum())
        frac = w_in / w_base_sum if w_base_sum > 0 else 0.0
        print(f'  [{label}] region={ireg} : {n_in} events, '
              f'sum(w)={w_in:.6e}  ({100*frac:.2f}% of base)')

        for nname, vals, color, ls in nucleons:
            counts, _ = np.histogram(vals[m], bins=edges, weights=w[m])
            ax.plot(centers, counts, color=color, ls=ls,
                    lw=1.4, label=nname)

        ax.set_xlim(args.p_min, args.p_max)
        ax.set_xlabel(r'Momentum (GeV/$c$)')
        ax.text(0.5, 1.02, label, transform=ax.transAxes,
                ha='center', va='bottom', fontsize=9)
        ax.set_ylim(bottom=0)

    axes[0].set_ylabel(r'Total weight')
    axes[-1].legend(loc='upper right', frameon=False, handlelength=1.6,
                    fontsize=8)

    fig.tight_layout()
    fig.savefig(out_path, dpi=args.dpi)
    png_path = os.path.splitext(out_path)[0] + '.png'
    fig.savefig(png_path, dpi=args.dpi)
    print(f'Saved {out_path}')
    print(f'Saved {png_path} (preview)')


if __name__ == '__main__':
    main()
