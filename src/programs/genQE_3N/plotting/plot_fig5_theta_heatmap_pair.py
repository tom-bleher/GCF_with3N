#!/usr/bin/env python3
"""Fig. 5 of the 3N-SRC paper.

Side-by-side theta12-theta23 heatmaps:

    Left:  3N+FSI generator
    Right: 2N+FSI generator

Both samples are read from detection-based `events` trees (built by
make_events_from_genT.py for the 3N genT, by SRC_analysis_2N for the 2N
sample) and pass through the SAME (e,e'pp) selection: exactly 2 detected
protons above kF, lead = detected proton closest to q with theta_pq < 8 deg
and |p_lead| > 0.75 |q|; the third (neutron) momentum is RECONSTRUCTED from
the missing momentum, p3 = -(p1 + p_recoil) with p1 = p_lead - q, and must
be above kF. theta12 = angle(p1, p_recoil), theta23 = angle(p_recoil, p3).

Cuts (paper convention, applied identically to BOTH samples):
  theta_e < 45 deg, Q^2 > 1 GeV^2, xB < 1.2,
  theta(p_lead, q) < 8 deg, |p_lead|/|q| > 0.75, 0.25 < pmiss < 0.9 GeV/c.

Two output variants are written: plain heatmaps (default) and a "_triangles"
version with the L (top-left) and R (top-right) regions outlined in solid
black.  Independent log color scales per panel, viridis colormap, paper styling
from paper_style.py.
"""
import argparse
import os
import sys
import numpy as np
import uproot
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Polygon

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paper_style import apply_style, figure_size, CMAP, OVERLAY_COLOR


mN = 0.93892


# ---------- region geometry (mirrors analyze_2N.py:region_params) ----------
def region_params(A=135.0, K=4.0):
    A_rad = np.radians(A)
    B = (np.degrees(np.arctan(np.sin(A_rad) / (K + np.cos(A_rad))))
         / (180.0 - A))
    return A, B


def _line_intersection(m1, b1, m2, b2):
    denom = m1 - m2
    if abs(denom) < 1e-15:
        return None
    x = (b2 - b1) / denom
    return float(x), float(m1 * x + b1)


def triangle_vertices_R(A=135.0, K=4.0):
    """Vertices of the top-right (Region R) triangle in (theta12, theta23)."""
    _, B = region_params(A, K)
    m1, b1 = B, 180.0 * (1.0 - B)
    m2, b2 = (1.0 / B), 180.0 * (1.0 - (1.0 / B))
    m3, b3 = -1.0, A + 180.0 - (180.0 - A) * B
    return [_line_intersection(m1, b1, m2, b2),
            _line_intersection(m1, b1, m3, b3),
            _line_intersection(m2, b2, m3, b3)]


def triangle_vertices_L(A=135.0, K=4.0):
    """Vertices of the top-left (Region L) triangle.

    Mirror of Region R about the x = (180 - y) axis -- we just reflect each
    vertex via (x, y) -> (360 - x - y, y) and confirm with two boundary lines.
    """
    _, B = region_params(A, K)
    m1, b1 = B, 180.0 * (1.0 - B)
    m2, b2 = (1.0 / B), 180.0 * (1.0 - (1.0 / B))
    m3, b3 = -1.0, A + 180.0 - (180.0 - A) * B
    s1_L = -(m1 / (1.0 + m1))
    s2_L = -(m2 / (1.0 + m2))
    b_L = 180.0
    x_vert_L = 360.0 - b3
    return [(0.0, 180.0),
            (x_vert_L, s1_L * x_vert_L + b_L),
            (x_vert_L, s2_L * x_vert_L + b_L)]


def _add_triangle(ax, points, *, edgecolor, lw=1.0, label_text=None):
    pts = [p for p in points
           if p is not None and np.isfinite(p[0]) and np.isfinite(p[1])]
    if len(pts) != 3:
        return
    cx = sum(p[0] for p in pts) / 3.0
    cy = sum(p[1] for p in pts) / 3.0
    pts = sorted(pts, key=lambda p: np.arctan2(p[1] - cy, p[0] - cx))
    ax.add_patch(Polygon(pts, closed=True, fill=False,
                         edgecolor=edgecolor, linewidth=lw))
    if label_text is not None:
        ax.text(cx, cy, label_text, ha='center', va='center',
                color=edgecolor, fontsize=8,
                bbox=dict(boxstyle='round,pad=0.15', facecolor='white',
                          edgecolor='none', alpha=0.7))


# ---------- angle helpers ----------
def angle_deg(ax, ay, az, bx, by, bz):
    dot = ax * bx + ay * by + az * bz
    a = np.sqrt(ax * ax + ay * ay + az * az)
    b = np.sqrt(bx * bx + by * by + bz * bz)
    return np.degrees(np.arccos(np.clip(dot / (a * b + 1e-30), -1.0, 1.0)))


# ---------- unified (e,e'pp) loader (3N and 2N samples alike) ----------
def load_events(path, label, theta_e_max, Q2_min, theta_pq_max,
                lead_over_q_min, pmiss_lo, pmiss_hi, xB_max, kf):
    """Load a detection-based `events` tree (make_events_from_genT.py /
    SRC_analysis_2N schema) and apply the paper (e,e'pp) selection —
    IDENTICAL for the 3N and 2N samples:

      * exactly 2 detected protons above kF (nProtonAboveKF == 2)
      * lead cuts: theta(p_lead, q) < 8 deg, |p_lead| > 0.75 |q|
      * 0.25 < pmiss < 0.9 GeV/c, xB < 1.2, Q2 > 1, theta_e < 45 deg
      * third nucleon RECONSTRUCTED from missing momentum:
        p3 = -(p1 + p_recoil) with p1 = p_lead - q; require |p3| > kF
        (all nucleons above the Fermi level, per the paper text)

    theta12 = angle(p1, p_recoil), theta23 = angle(p_recoil, p3).
    """
    print(f'Reading {label}: {path}')
    arr = uproot.open(path)['events'].arrays(
        ['weight', 'lead', 'recoil', 'q',
         'Q2', 'xB', 'pmiss', 'scattering_angle', 'nProtonAboveKF'],
        library='np',
    )
    w = arr['weight']
    Q2, xB, pmiss, theta_e = arr['Q2'], arr['xB'], arr['pmiss'], arr['scattering_angle']
    nP = arr['nProtonAboveKF']
    lead, rec, q = arr['lead'], arr['recoil'], arr['q']

    p1x, p1y, p1z = lead[:, 0] - q[:, 0], lead[:, 1] - q[:, 1], lead[:, 2] - q[:, 2]
    p3x = -(p1x + rec[:, 0])
    p3y = -(p1y + rec[:, 1])
    p3z = -(p1z + rec[:, 2])
    p3_mag = np.sqrt(p3x**2 + p3y**2 + p3z**2)

    lead_mag = np.sqrt(lead[:, 0]**2 + lead[:, 1]**2 + lead[:, 2]**2)
    q_mag = np.sqrt(q[:, 0]**2 + q[:, 1]**2 + q[:, 2]**2)
    lead_over_q = lead_mag / (q_mag + 1e-30)
    lead_angle_q = angle_deg(lead[:, 0], lead[:, 1], lead[:, 2],
                             q[:, 0], q[:, 1], q[:, 2])

    theta12 = angle_deg(p1x, p1y, p1z, rec[:, 0], rec[:, 1], rec[:, 2])
    theta23 = angle_deg(rec[:, 0], rec[:, 1], rec[:, 2], p3x, p3y, p3z)

    mask = (np.isfinite(w) & (w > 0)
            & (nP == 2)
            & (theta_e < theta_e_max) & (Q2 > Q2_min)
            & (xB < xB_max)
            & (lead_angle_q < theta_pq_max)
            & (lead_over_q > lead_over_q_min)
            & (pmiss > pmiss_lo) & (pmiss < pmiss_hi)
            & (p3_mag > kf))

    print(f'  total events: {len(w)} ; (e,e\'pp) events passing cuts: {int(mask.sum())}')
    return theta12[mask], theta23[mask], w[mask]


# ---------- region fractions ----------
def in_region_R_array(theta12, theta23, A=135.0, K=4.0):
    _, B = region_params(A, K)
    line1 = B * (theta12 - 180.0) + 180.0
    line2 = (1.0 / B) * (theta12 - 180.0) + 180.0
    line3 = -(theta12 - A) + 180.0 - (180.0 - A) * B
    return (theta23 <= line1) & (theta23 >= line2) & (theta23 >= line3)


def in_region_L_array(theta12, theta23, A=135.0, K=4.0):
    return in_region_R_array(360.0 - theta12 - theta23, theta23, A, K)


# ---------- plotting ----------
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

    ax.set_xlim(0, 180)
    ax.set_ylim(0, 180)
    ax.set_aspect('equal')
    ax.set_xticks(np.arange(0, 181, 30))
    ax.set_yticks(np.arange(0, 181, 30))
    ax.set_xlabel(r'$\theta_{12}$ (deg)')
    return im


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--input-3N', default='events_3N.root')
    p.add_argument('--input-2N',
                   default='../genQE_FSI/events/misc/events_2N.root')
    p.add_argument('--output', default='analysis_output/png_files/Paper plots/fig5_theta_heatmap_pair.pdf')
    p.add_argument('--bins-3N', type=int, default=400)
    p.add_argument('--bins-2N', type=int, default=200)
    p.add_argument('--dpi', type=int, default=400)
    # cuts
    p.add_argument('--theta-e-max', type=float, default=45.0)
    p.add_argument('--Q2-min', type=float, default=1.0)
    p.add_argument('--xB-max', type=float, default=1.2)
    p.add_argument('--theta-pq-max', type=float, default=8.0)
    p.add_argument('--lead-over-q-min', type=float, default=0.75)
    p.add_argument('--pmiss-lo', type=float, default=0.25)
    p.add_argument('--pmiss-hi', type=float, default=0.90)
    p.add_argument('--kF', type=float, default=0.22,
                   help='Fermi level for the reconstructed-third cut (GeV/c)')
    # color scales (independent per panel)
    p.add_argument('--vmin-3N', type=float, default=1e-9)
    p.add_argument('--vmax-3N', type=float, default=2e-2)
    p.add_argument('--vmin-2N', type=float, default=1e-10)
    p.add_argument('--vmax-2N', type=float, default=3e-4)
    args = p.parse_args()

    apply_style()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    root_dir = os.path.abspath(os.path.join(script_dir, '..'))
    in3 = args.input_3N if os.path.isabs(args.input_3N) else os.path.join(root_dir, args.input_3N)
    in2 = args.input_2N if os.path.isabs(args.input_2N) else os.path.join(root_dir, args.input_2N)
    out_path = args.output if os.path.isabs(args.output) else os.path.join(root_dir, args.output)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    cut_kwargs = dict(theta_e_max=args.theta_e_max, Q2_min=args.Q2_min,
                      theta_pq_max=args.theta_pq_max,
                      lead_over_q_min=args.lead_over_q_min,
                      pmiss_lo=args.pmiss_lo, pmiss_hi=args.pmiss_hi,
                      xB_max=args.xB_max, kf=args.kF)

    t12_3, t23_3, w_3 = load_events(in3, '3N+FSI', **cut_kwargs)
    t12_2, t23_2, w_2 = load_events(in2, '2N+FSI', **cut_kwargs)

    # L/R weight fractions on the cut sample (for the caption / paper Fig. 6)
    for name, t12, t23, w in [('3N+FSI', t12_3, t23_3, w_3),
                              ('2N+FSI', t12_2, t23_2, w_2)]:
        wt = w.sum() + 1e-30
        inL = in_region_L_array(t12, t23)
        inR = in_region_R_array(t12, t23)
        fL = w[inL].sum() / wt
        fR = w[inR].sum() / wt
        ratio = (fL / fR) if fR > 0 else float('inf')
        print(f'  [{name}] f_L={fL:.4f}  f_R={fR:.4f}  L/R={ratio:.3f}')
        # Weighted-MC statistical uncertainties (L and R are disjoint, so the
        # sum-of-w^2 errors on N_L and N_R are independent). Kept off the
        # [name] line, which make_paper_figures.sh parses for L/R.
        s2L, s2R = (w[inL]**2).sum(), (w[inR]**2).sum()
        s2 = (w**2).sum()
        efL = np.sqrt(s2L * (1 - 2 * fL) + s2 * fL**2) / wt
        efR = np.sqrt(s2R * (1 - 2 * fR) + s2 * fR**2) / wt
        eratio = (ratio * np.sqrt(s2L / (w[inL].sum()**2 + 1e-300)
                                  + s2R / (w[inR].sum()**2 + 1e-300))
                  if fR > 0 else float('nan'))
        print(f'    stat. unc. {name}: f_L={fL:.5f}+-{efL:.5f}  '
              f'f_R={fR:.5f}+-{efR:.5f}  L/R={ratio:.4f}+-{eratio:.4f}  '
              f'(n_L={int(inL.sum())}, n_R={int(inR.sum())}, n={len(w)})')

    # ---- Build figure (two panels, two-column PRR width) ----
    tri_R = triangle_vertices_R()
    tri_L = triangle_vertices_L()

    for variant, draw_triangles in [('', False), ('_triangles', True)]:
        fig, axes = plt.subplots(1, 2, sharey=True,
                                 figsize=figure_size(cols=2, ratio=0.50))
        im1 = plot_panel(axes[0], t12_3, t23_3, w_3, args.bins_3N,
                         args.vmin_3N, args.vmax_3N)
        im2 = plot_panel(axes[1], t12_2, t23_2, w_2, args.bins_2N,
                         args.vmin_2N, args.vmax_2N)

        axes[0].set_ylabel(r'$\theta_{23}$ (deg)')
        axes[1].tick_params(labelleft=False)

        cb1 = fig.colorbar(im1, ax=axes[0], pad=0.02, fraction=0.046)
        cb2 = fig.colorbar(im2, ax=axes[1], label='Normalized weight',
                           pad=0.02, fraction=0.046)
        for cb in (cb1, cb2):
            cb.ax.tick_params(labelsize=8)
        cb2.ax.yaxis.label.set_size(9)

        if draw_triangles:
            for ax in axes:
                _add_triangle(ax, tri_L, edgecolor='k', lw=1.0, label_text='L')
                _add_triangle(ax, tri_R, edgecolor='k', lw=1.0, label_text='R')

        fig.tight_layout()
        out_pdf = out_path.replace('.pdf', f'{variant}.pdf')
        out_png = out_pdf.replace('.pdf', '.png')
        fig.savefig(out_pdf, dpi=args.dpi)
        fig.savefig(out_png, dpi=args.dpi)
        plt.close(fig)
        print(f'Saved {out_pdf}')
        print(f'Saved {out_png} (preview)')


if __name__ == '__main__':
    main()
