#!/usr/bin/env python3
"""Opening angle between the odd proton and the (p,n) relative momentum,
computed straight from the 3N density tables -- no generator, no kinematics.

See pd_jacobi.py for the Jacobi conventions and the cell integration.

Row 1  full table support : the whole wavefunction.
Row 2  generator support  : all three internal momenta above
                            fkMinInternal = 0.25 GeV/c, which is what
                            QEGenerator_3N.cc:224 actually keeps, so this is the
                            like-for-like row to compare with generated events.
Column 3 overlays the two choices of which proton is the odd one.  They must
coincide: the two protons are identical, swapping them is theta -> 180-theta,
and every table is symmetric under that to ~1e-6 (storage rounding).

Writes opening_angle_density.png and opening_angle_density_summary.txt .
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pd_jacobi as J

HERE = os.path.dirname(os.path.abspath(__file__))
BINS = np.linspace(0, 180, 91)
XL = r"$\theta$ between the odd proton and $(\vec p_p-\vec p_n)/2$  [deg]"

cache = {}
for name, fn in J.TABLES:
    for pd_only in (False, True):
        cache[(name, pd_only)] = J.histograms(J.load(fn, pd_only), BINS)
        print(f"  done {'pD ' if pd_only else '   '}{name}")

fig, axes = plt.subplots(2, 3, figsize=(16.5, 9))
out = []
for row, support in enumerate([False, True]):
    sfx = "_sup" if support else ""
    for col, pd_only in enumerate([False, True]):
        ax = axes[row][col]
        for name, _ in J.TABLES:
            d = cache[(name, pd_only)]["thA" + sfx]
            ax.stairs(d["h"], BINS, color=J.COL[name], lw=1.9, label=name)
            out.append(f"{'pD  ' if pd_only else 'full'} {name:11s} "
                       f"{'gen-support ' if support else 'full-support'} "
                       f"<theta_p,(pn)> = {d['mean']:6.2f} deg")
        ax.set_title(f"{'pD component only' if pd_only else 'full wavefunction'}"
                     f"   ({'generator support: all $|p_i|>0.25$' if support else 'full table support'})",
                     fontsize=11)
        ax.set_xlabel(XL); ax.set_ylabel("normalised")
        ax.set_xlim(0, 180); ax.set_ylim(bottom=0); ax.legend(frameon=False, fontsize=8)

    ax, d = axes[row][2], cache[("AV18", False)]
    hA, hB, hn = (d[k + sfx]["h"] for k in ("thA", "thB", "th_native"))
    ax.stairs(hA, BINS, color="#2b6cb0", lw=2.8, label=r"choice A: odd proton = $v_2$")
    ax.stairs(hB, BINS, color="#dd8800", lw=1.4, ls="--", label=r"choice B: odd proton = $v_3$")
    # th_native IS the table's theta grid (61 values), so it aliases against the
    # 90 output bins -- rebin by 3 for a readable reference curve
    hn3 = hn.reshape(-1, 3).mean(axis=1)
    ax.stairs(hn3, BINS[::3], color="0.45", lw=1.5, ls=":",
              label=r"table axis: neutron vs $(\vec p_p-\vec p_p)/2$")
    dmax = np.abs(hA - hB).max() / max(hA.max(), 1e-300)
    ax.set_title(f"AV18: A vs B   (max |A-B|/max = {dmax:.1e})", fontsize=11)
    ax.set_xlabel(r"$\theta$ [deg]"); ax.set_ylabel("normalised")
    ax.set_xlim(0, 180); ax.set_ylim(bottom=0); ax.legend(frameon=False, fontsize=8)
    out.append(f"A/B identity check ({'generator support' if support else 'full support'}): "
               f"max|A-B|/max(A) = {dmax:.3e}")

fig.suptitle(r"3N density: opening angle between the odd proton and the relative momentum of the remaining $(p,n)$ pair"
             "\nfrom the density tables alone (no generator, no kinematics). Tables are $^3$He 3N densities, "
             r"used as the universal 3N SRC density; the generator's target is $^{12}$C", fontsize=12.5)
fig.tight_layout(rect=(0, 0, 1, 0.94))
png = os.path.join(HERE, "opening_angle_density.png")
fig.savefig(png, dpi=150)
open(os.path.join(HERE, "opening_angle_density_summary.txt"), "w").write("\n".join(out) + "\n")
print("\n".join(out)); print("saved", png)
