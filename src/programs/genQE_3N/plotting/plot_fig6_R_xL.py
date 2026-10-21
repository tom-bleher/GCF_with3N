#!/usr/bin/env python3
"""Fig. 6 of the 3N-SRC paper.

Adapted from genQE_FSI/analysis/plot_R_xL.py.

For an event mixture of 2N+FSI background and 3N-SRC signal, with

    x_L = N_{L,2N} / N_{L,3N}

the 2N-to-3N count RATIO in the L region of the (theta12, theta23) plane
(not a fraction: x_L = 0 is pure 3N, x_L -> infinity pure 2N), the observed
L/R event-count ratio is

    R(x_L) = (1 + x_L) / (1/R_3N + x_L/R_2N)

where R_3N = N_{L,3N} / N_{R,3N} and R_2N = N_{L,2N} / N_{R,2N} are the
"pure-signal" and "pure-background" L/R ratios obtained from Fig. 5.

make_paper_figures.sh passes the R values printed by
plot_fig5_theta_heatmap_pair.py. The defaults are the shipped ones
(Paper_Figures_indepFSI/Fig5_summary.txt):
    R_3N = 0.451   (3N+FSI)
    R_2N = 0.107   (2N+FSI, detection-based)
"""
import argparse
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paper_style import apply_style, figure_size


# Default values from the shipped Fig. 5 (Paper_Figures_indepFSI/Fig5_summary.txt,
# paper cuts); make_paper_figures.sh overrides them with the values it measures.
R_3N_DEFAULT = 0.451
R_2N_DEFAULT = 0.107


def R(xL, R_3N, R_2N):
    return (1.0 + xL) / (1.0 / R_3N + xL / R_2N)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--output',
                   default='analysis_output/png_files/Paper plots/fig6_R_xL.pdf')
    p.add_argument('--R-3N', type=float, default=R_3N_DEFAULT)
    p.add_argument('--R-2N', type=float, default=R_2N_DEFAULT)
    p.add_argument('--xL-max', type=float, default=0.5)
    p.add_argument('--dpi', type=int, default=400)
    p.add_argument('--cols', type=int, default=1, choices=(1, 2))
    args = p.parse_args()

    apply_style()

    print(f'  R_3N = {args.R_3N:.4f}  (= N_L,3N / N_R,3N)')
    print(f'  R_2N = {args.R_2N:.4f}  (= N_L,2N / N_R,2N)')
    print(f'  R(xL=0)    = {R(0.0, args.R_3N, args.R_2N):.4f}')
    print(f'  R(xL=0.10) = {R(0.10, args.R_3N, args.R_2N):.4f}')
    print(f'  R(xL=0.50) = {R(args.xL_max, args.R_3N, args.R_2N):.4f}')

    script_dir = os.path.dirname(os.path.abspath(__file__))
    root_dir = os.path.abspath(os.path.join(script_dir, '..'))
    out_path = args.output if os.path.isabs(args.output) else os.path.join(root_dir, args.output)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    xL = np.linspace(0.0, args.xL_max, 500)
    y = R(xL, args.R_3N, args.R_2N)

    fig, ax = plt.subplots(figsize=figure_size(cols=args.cols, ratio=0.78))
    ax.plot(xL, y, 'k-', linewidth=1.6)

    # No-FSI marker (xL = 0)
    ax.plot(0.0, R(0.0, args.R_3N, args.R_2N), 'o',
            color='#2471a3', markersize=6, zorder=5,
            markeredgecolor='k', markeredgewidth=0.6)
    ax.annotate('No FSI', xy=(0.0, R(0.0, args.R_3N, args.R_2N)),
                xytext=(0.03, R(0.0, args.R_3N, args.R_2N) - 0.06),
                fontsize=9, color='#2471a3', va='top')

    # 10%-FSI marker
    xL_fsi = 0.10
    ax.plot(xL_fsi, R(xL_fsi, args.R_3N, args.R_2N), 'o',
            color='#c0392b', markersize=6, zorder=5,
            markeredgecolor='k', markeredgewidth=0.6)
    ax.annotate(r'$10\%$ FSI', xy=(xL_fsi, R(xL_fsi, args.R_3N, args.R_2N)),
                xytext=(xL_fsi + 0.02, R(xL_fsi, args.R_3N, args.R_2N) + 0.06),
                fontsize=9, color='#c0392b', va='bottom')

    ax.set_xlabel(r'$x_L \equiv N_{2N} / N_{3N}$ in region L')
    ax.set_ylabel(r'$R \equiv N_L / N_R$')
    ax.set_xlim(-0.02, args.xL_max + 0.02)
    ax.set_ylim(0.0, max(1.05, 1.05 * args.R_3N))

    # Annotate the input R values in a small box
    txt = (rf'$R_{{3N}} = {args.R_3N:.3f}$' '\n'
           rf'$R_{{2N}} = {args.R_2N:.3f}$')
    ax.text(0.97, 0.97, txt, transform=ax.transAxes,
            ha='right', va='top', fontsize=8,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                      edgecolor='0.7', linewidth=0.5, alpha=0.9))

    fig.tight_layout()
    fig.savefig(out_path, dpi=args.dpi)
    png_path = os.path.splitext(out_path)[0] + '.png'
    fig.savefig(png_path, dpi=args.dpi)
    print(f'Saved {out_path}')
    print(f'Saved {png_path} (preview)')


if __name__ == '__main__':
    main()
