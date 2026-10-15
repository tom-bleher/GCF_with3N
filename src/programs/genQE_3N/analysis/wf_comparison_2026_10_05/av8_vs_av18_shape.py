#!/usr/bin/env python3
"""Is old AV8_1D.dnst proportional to the new Faddeev AV18 table, and where?

Writes av8_vs_av18_shape.png: the pointwise ratio AV8/AV18 divided by the ratio
of their full integrals C, as a function of k_a and of k_b (density-weighted,
i.e. averaged with the AV18 density so the physically populated region dominates).
Proportionality <=> curve flat at 1.
"""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ARR = "/home/ilayws/Ilay-Generators/src/programs/genQE_3N/array"
OUT = "/home/ilayws/Ilay-Generators/src/programs/genQE_3N/analysis/wf_comparison_2026_10_05/av8_vs_av18_shape.png"
NT, NK, HBARC = 61, 101, 0.1973269804
th = np.radians(np.linspace(0, 180, NT)); k = np.linspace(0, 10, NK)
load = lambda fn: np.loadtxt(f"{ARR}/{fn}").reshape(NT, NK, NK)
a18, a8 = load("he3_3Ndensity_AV18.dnst"), load("AV8_1D.dnst")

w = (k**2)[None, :, None]*(k**2)[None, None, :]*np.sin(th)[:, None, None]
tri = lambda x: np.trapezoid(np.trapezoid(np.trapezoid(x, k, axis=2), k, axis=1), th)
C = tri(a8*w)/tri(a18*w)
dens = a18*w
r = np.where((a18 > 0) & (a8 > 0), a8/np.maximum(a18, 1e-300), np.nan)/C

fig, axes = plt.subplots(1, 3, figsize=(15, 4.4))
for ax, axis_sum, lab in ((axes[0], (0, 2), r"$k_a$  (neutron momentum)"),
                          (axes[1], (0, 1), r"$k_b$  ($pp$ half-relative)")):
    num = np.nansum(np.where(np.isfinite(r), r*dens, 0.0), axis=axis_sum)
    den = np.nansum(np.where(np.isfinite(r), dens, 0.0), axis=axis_sum)
    frac = dens.sum(axis=axis_sum)/dens.sum()
    ax.axhline(1.0, color="0.6", lw=1, ls="--")
    ax.plot(k*HBARC, num/np.maximum(den, 1e-300), color="#cc3355", lw=2)
    ax2 = ax.twinx()
    ax2.fill_between(k*HBARC, frac/frac.max(), color="#5588cc", alpha=0.18, lw=0)
    ax2.set_ylim(0, 3); ax2.set_yticks([])
    ax.set_xlabel(lab + "  [GeV/c]"); ax.set_xlim(0, 1.2); ax.set_ylim(0, 1.6)
    ax.set_ylabel(r"$(\rho_{AV8}/\rho_{AV18})\,/\,C$")
    ax.text(0.03, 0.06, "shading: fraction of the\nAV18 norm (arb.)", transform=ax.transAxes, fontsize=8, color="#336699")

ia = np.argmin(abs(k*HBARC - 0.0)) # 2D map over (k_a,k_b), theta-averaged
num = np.nansum(np.where(np.isfinite(r), r*dens, 0.0), axis=0)
den = np.nansum(np.where(np.isfinite(r), dens, 0.0), axis=0)
m = np.where(den > 0, num/den, np.nan)
im = axes[2].pcolormesh(k*HBARC, k*HBARC, m.T, vmin=0.4, vmax=1.6, cmap="RdBu_r")
axes[2].set_xlabel(r"$k_a$ [GeV/c]"); axes[2].set_ylabel(r"$k_b$ [GeV/c]")
axes[2].set_xlim(0, 1.2); axes[2].set_ylim(0, 1.2)
fig.colorbar(im, ax=axes[2], label=r"$(\rho_{AV8}/\rho_{AV18})/C$")
fig.suptitle(f"old AV8_1D.dnst vs Faddeev AV18:  C = integral ratio = {C:.4g}"
             rf"  $\simeq 2\pi(\hbar c)^6$ = {2*np.pi*HBARC**6:.4g}", fontsize=11)
fig.tight_layout()
fig.savefig(OUT, dpi=160)
print("C =", C, " C/2pi(hbarc)^6 =", C/(2*np.pi*HBARC**6))
print("saved", OUT)
