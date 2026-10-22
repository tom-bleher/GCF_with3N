#!/usr/bin/env python3
"""Overlay L/R/Star regions on the 3N/2N heatmaps and print fractions incl. 'other'."""
import os, sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_theta_heatmap_3N_vs_2N import load_3N, load_2N
from plot_fig5_theta_heatmap_pair import (in_region_L_array, in_region_R_array,
                                          triangle_vertices_L, triangle_vertices_R)
from report_region_weights import in_star, STAR_CENTER, STAR_RADIUS

root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
in3 = os.path.join(root, 'events/3N_FSI_15M_12C.root')
in2 = os.path.join(root, '../genQE_FSI/events/misc/events_2N.root')

common = dict(theta_e_max=45.0, Q2_min=1.0, xB_max=1.2, kF_thr=0.25,
              mode='eq2', theta_e_min=8.0)
t12_3, t23_3, w3 = load_3N(in3, 6.0, **common)
t12_2, t23_2, w2 = load_2N(in2, **common)

def fracs(t12, t23, w):
    tot = w.sum()
    L = w[in_region_L_array(t12, t23)].sum()/tot
    R = w[in_region_R_array(t12, t23)].sum()/tot
    S = w[in_star(t12, t23)].sum()/tot
    return L, R, S, 1.0 - L - R - S

for name, t12, t23, w in [('3N', t12_3, t23_3, w3), ('2N', t12_2, t23_2, w2)]:
    L, R, S, O = fracs(t12, t23, w)
    # L, at (0, 180): p1 || p2 with the reconstructed third nucleon opposite,
    # so p3 is the rocket. R, at (180, 180): the recoil p2 is the rocket.
    print(f'{name}: L(3rd-rkt)={L:6.2%}  R(recoil-rkt)={R:6.2%}  '
          f'Star={S:6.2%}  Other={O:6.2%}  (sum={L+R+S+O:.2%})')

# ---- overlay figure ----
def panel(ax, t12, t23, w, title):
    h, xe, ye = np.histogram2d(t12, t23, bins=36, range=[[0,180],[0,180]], weights=w)
    h = h/(h.sum()+1e-30)
    ax.pcolormesh(xe, ye, np.where(h>0,h,np.nan).T, cmap='viridis',
                  norm=LogNorm(vmin=1e-5, vmax=5e-2), shading='auto')
    x = np.linspace(0,180,200)
    for yy in (180-x/2, 360-2*x, x):
        ax.plot(x, yy, 'k--', lw=0.6, alpha=0.6)
    # L/R triangles
    for verts, col, lab in [(triangle_vertices_L(), 'red', 'L'),
                            (triangle_vertices_R(), 'magenta', 'R')]:
        vv = np.array(verts + [verts[0]])
        ax.plot(vv[:,0], vv[:,1], color=col, lw=1.6)
        cx, cy = np.mean([v[0] for v in verts]), np.mean([v[1] for v in verts])
        ax.text(cx, cy, lab, color=col, fontsize=11, ha='center', va='center', weight='bold')
    th = np.linspace(0,2*np.pi,100)
    ax.plot(STAR_CENTER[0]+STAR_RADIUS*np.cos(th),
            STAR_CENTER[1]+STAR_RADIUS*np.sin(th), color='cyan', lw=1.6)
    ax.text(STAR_CENTER[0], STAR_CENTER[1], 'Star', color='cyan', fontsize=9,
            ha='center', va='center', weight='bold')
    ax.set_xlim(0,180); ax.set_ylim(0,180); ax.set_aspect('equal')
    ax.set_xlabel(r'$\theta_{12}$ (deg)'); ax.set_title(title, fontsize=9)

fig, axes = plt.subplots(1,2, figsize=(11,5), sharey=True)
panel(axes[0], t12_3, t23_3, w3, '3N+FSI 12C (15M)')
panel(axes[1], t12_2, t23_2, w2, '2N+FSI')
axes[0].set_ylabel(r'$\theta_{23}$ (deg)')
out = os.path.join(root, 'analysis/Plots/region_overlay_check.png')
fig.savefig(out, dpi=150, bbox_inches='tight')
print('Saved', out)
