#!/usr/bin/env python3
"""R(xL): observed L/R event-count ratio vs the 2N-vs-3N mixing parameter xL.

The default pure-sample L/R ratios come from an OLD per-region report with a
looser lead cut (eq2 mode, theta_e in [8, 45], Q^2 > 1, theta_pq < 12 deg):

    R_3N  =  f_L^3N / f_R^3N  =  24.6 / 2.8    = 8.786
    R_2N  =  f_L^2N / f_R^2N  =  0.11 / 11.0   = 0.010

They are NOT the paper's values. The paper cuts (theta_pq < 8 deg, see
Paper_Figures_indepFSI/Fig5_summary.txt) give R_3N = 0.451, R_2N = 0.107;
pass those with --R-3N/--R-2N (or use plot_fig6_R_xL.py). The plots label
the defaults "theta_pq < 12 deg".

For an event mixture with relative 2N-to-3N weighting parameterised by xL,

    R(xL) = (1 + xL) / (1/R_3N + xL/R_2N)

so xL=0 -> pure 3N (R -> R_3N), xL->infty -> pure 2N (R -> R_2N).
"""
import argparse
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paper_style import apply_style, figure_size, LINE_COLORS


# L/R fractions from the eq2_no_cuts region report (theta_pq < 12 deg; NOT the
# paper's theta_pq < 8 deg selection -- see the module docstring)
F_L_3N_PCT = 24.6      # per cent of total 3N weight inside Region L
F_R_3N_PCT = 2.8
F_L_2N_PCT = 0.11
F_R_2N_PCT = 11.0


def R_curve(xL, R_3N, R_2N):
    return (1.0 + xL) / (1.0 / R_3N + xL / R_2N)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--output',
                   default='analysis/Plots/R_xL.pdf')
    p.add_argument('--R-3N', type=float, default=F_L_3N_PCT / F_R_3N_PCT)
    p.add_argument('--R-2N', type=float, default=F_L_2N_PCT / F_R_2N_PCT)
    p.add_argument('--xL-min', type=float, default=1e-4)
    p.add_argument('--xL-max', type=float, default=1e3)
    p.add_argument('--npts',   type=int,   default=400)
    p.add_argument('--dpi',    type=int,   default=400)
    args = p.parse_args()

    apply_style()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    root_dir   = os.path.abspath(os.path.join(script_dir, '..'))
    out_path = args.output if os.path.isabs(args.output) else os.path.join(root_dir, args.output)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    R_3N = args.R_3N
    R_2N = args.R_2N
    print(f'R_3N = {R_3N:.4f}    R_2N = {R_2N:.4f}')

    # ------------- Log-scale version -----------------------------------------
    xL_log = np.logspace(np.log10(args.xL_min), np.log10(args.xL_max), args.npts)
    R_log  = R_curve(xL_log, R_3N, R_2N)

    fig, ax = plt.subplots(figsize=figure_size(cols=1, ratio=0.78))
    ax.plot(xL_log, R_log, color=LINE_COLORS[0], lw=1.5)

    # Horizontal asymptotes (no text annotations to avoid clutter)
    ax.axhline(R_3N, color='0.4', lw=0.6, ls=':')
    ax.axhline(R_2N, color='0.4', lw=0.6, ls=':')

    # Compact stats box anchored away from the curve
    txt = (f'$R_{{3N}}={R_3N:.2f}$\n'
           f'$R_{{2N}}={R_2N:.3f}$')
    if (R_3N, R_2N) == (F_L_3N_PCT / F_R_3N_PCT, F_L_2N_PCT / F_R_2N_PCT):
        txt += '\n' r'($\theta_{pq}<12^\circ$ defaults)'
    ax.text(0.97, 0.05, txt, transform=ax.transAxes,
            ha='right', va='bottom', fontsize=8,
            bbox=dict(boxstyle='round,pad=0.3', fc='white',
                      ec='0.6', lw=0.5, alpha=0.9))

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(args.xL_min, args.xL_max)
    ax.set_ylim(min(R_2N, R_3N) * 0.5, max(R_2N, R_3N) * 2.0)
    ax.set_xlabel(r'$x_L$')
    ax.set_ylabel(r'$R(x_L) = N_L / N_R$')
    ax.set_title('Mixing curve in Region L vs Region R (log-log)',
                 fontsize=9)
    ax.grid(True, which='both', linewidth=0.25, alpha=0.35)
    fig.tight_layout()
    fig.savefig(out_path, dpi=args.dpi)
    fig.savefig(os.path.splitext(out_path)[0] + '.png', dpi=args.dpi)
    plt.close(fig)
    print(f'Saved (log)    {out_path}')

    # ------------- Linear-scale version --------------------------------------
    # The transition happens around xL ~ R_2N/R_3N ~ 1e-3, so pick a linear
    # x-range that shows the transition clearly without compressing it.
    xL_max_lin = max(0.5, 50.0 * R_2N / R_3N)  # ~30x the half-transition point
    xL_lin = np.linspace(0.0, xL_max_lin, args.npts)
    R_lin  = R_curve(xL_lin, R_3N, R_2N)

    fig, ax = plt.subplots(figsize=figure_size(cols=1, ratio=0.78))
    ax.plot(xL_lin, R_lin, color=LINE_COLORS[0], lw=1.5)
    ax.axhline(R_3N, color='0.4', lw=0.6, ls=':')
    ax.axhline(R_2N, color='0.4', lw=0.6, ls=':')

    ax.text(0.97, 0.95, txt, transform=ax.transAxes,
            ha='right', va='top', fontsize=8,
            bbox=dict(boxstyle='round,pad=0.3', fc='white',
                      ec='0.6', lw=0.5, alpha=0.9))

    ax.set_xlim(0.0, xL_max_lin)
    ax.set_ylim(0.0, R_3N * 1.1)
    ax.set_xlabel(r'$x_L$')
    ax.set_ylabel(r'$R(x_L) = N_L / N_R$')
    ax.set_title('Mixing curve in Region L vs Region R (linear)',
                 fontsize=9)
    ax.grid(True, which='both', linewidth=0.25, alpha=0.35)
    fig.tight_layout()
    out_lin = os.path.splitext(out_path)[0] + '_linear.pdf'
    fig.savefig(out_lin, dpi=args.dpi)
    fig.savefig(os.path.splitext(out_lin)[0] + '.png', dpi=args.dpi)
    plt.close(fig)
    print(f'Saved (linear) {out_lin}')

    # ------------- Linear-scale R vs 1/xL ------------------------------------
    # 1/xL = (3N count in Region L) / (2N count in Region L).
    # Transition midpoint at 1/xL = R_3N/R_2N; pick a range that goes a few
    # times past it so the R_3N asymptote is visible.
    inv_max = 3.0 * R_3N / R_2N
    inv_xL = np.linspace(1e-6, inv_max, args.npts)  # 1/xL > 0
    xL_for_inv = 1.0 / inv_xL
    R_inv = R_curve(xL_for_inv, R_3N, R_2N)

    fig, ax = plt.subplots(figsize=figure_size(cols=1, ratio=0.78))
    ax.plot(inv_xL, R_inv, color=LINE_COLORS[0], lw=1.5)
    ax.axhline(R_3N, color='0.4', lw=0.6, ls=':')
    ax.axhline(R_2N, color='0.4', lw=0.6, ls=':')

    # Stats box anchored bottom-right (curve rises from left, asymptotes top)
    ax.text(0.97, 0.05, txt, transform=ax.transAxes,
            ha='right', va='bottom', fontsize=8,
            bbox=dict(boxstyle='round,pad=0.3', fc='white',
                      ec='0.6', lw=0.5, alpha=0.9))

    ax.set_xlim(0.0, inv_max)
    ax.set_ylim(0.0, R_3N * 1.1)
    ax.set_xlabel(r'$1/x_L$  =  3N / 2N count ratio in Region L')
    ax.set_ylabel(r'$R(x_L) = N_L / N_R$')
    ax.set_title('Mixing curve in Region L vs Region R (linear, vs 1/$x_L$)',
                 fontsize=9)
    ax.grid(True, which='both', linewidth=0.25, alpha=0.35)
    fig.tight_layout()
    out_inv = os.path.splitext(out_path)[0] + '_linear_invxL.pdf'
    fig.savefig(out_inv, dpi=args.dpi)
    fig.savefig(os.path.splitext(out_inv)[0] + '.png', dpi=args.dpi)
    plt.close(fig)
    print(f'Saved (linear vs 1/xL) {out_inv}')


if __name__ == '__main__':
    main()
