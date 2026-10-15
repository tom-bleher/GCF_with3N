#!/usr/bin/env python3
"""Direct comparison of the 3N density tables (no generator involved).

Grid (QEGenerator_3N.hh): theta 61 x [0,180] deg, k_cm 101 x [0,10] fm^-1,
k_rel 101 x [0,10] fm^-1, one value per line, ordered theta -> k_cm -> k_rel.

New tables (u=1..5 full, u=6..10 pD-only) store the RAW unpolarized density
rho(theta, k_cm, k_rel). The old AV8/N2LO/G3/AV4 tables use an unknown
convention, so their normalisation is NOT comparable -- only shapes are, and
even those only loosely.
"""
import os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ARR = "/home/ilayws/Ilay-Generators/src/programs/genQE_3N/array"
OUT = os.path.dirname(os.path.abspath(__file__))
NT, NC, NR = 61, 101, 101
th = np.linspace(0, 180, NT); kcm = np.linspace(0, 10, NC); krel = np.linspace(0, 10, NR)

TABLES = [("AV18 (u=1)", "he3_3Ndensity_AV18.dnst", "#3366bb", "-"),
          ("AV18+UIX (u=2)", "he3_3Ndensity_AV18UIX.dnst", "#5588cc", "-"),
          ("CD-Bonn (u=4)", "he3_3Ndensity_CDBonn.dnst", "#118844", "-"),
          ("pD AV18 (u=6)", "he3_3Ndensity_pD_AV18.dnst", "#cc3355", "--"),
          ("pD AV18+UIX (u=7)", "he3_3Ndensity_pD_AV18UIX.dnst", "#dd7799", "--"),
          ("pD CD-Bonn (u=9)", "he3_3Ndensity_pD_CDBonn.dnst", "#ee9944", "--"),
          ("AV8 (old)", "AV8_1D.dnst", "#888888", ":")]


def load(fn):
    v = np.loadtxt(os.path.join(ARR, fn), usecols=0)
    assert v.size == NT*NC*NR, (fn, v.size)
    return v.reshape(NT, NC, NR)


def measure_integral(rho):
    """int rho * k_cm^2 k_rel^2 sin(theta) dtheta dk_cm dk_rel (the phase-space
    measure get_rho() multiplies in, up to constant factors)."""
    w = (kcm**2)[None, :, None] * (krel**2)[None, None, :] * np.sin(np.radians(th))[:, None, None]
    return np.trapezoid(np.trapezoid(np.trapezoid(rho*w, krel, axis=2), kcm, axis=1), np.radians(th))


data, lines = {}, []
for lab, fn, _, _ in TABLES:
    p = os.path.join(ARR, fn)
    if not os.path.exists(p):
        print("missing", fn); continue
    data[lab] = load(fn)
    print(f"loaded {lab:20s} {fn}", flush=True)

lines.append("=== table integrals  I = int rho k_cm^2 k_rel^2 sin(th) dth dk_cm dk_rel ===")
ref = measure_integral(data["AV18 (u=1)"])
for lab in data:
    I = measure_integral(data[lab])
    lines.append(f"  {lab:20s} I = {I:12.5g}   I/I(AV18) = {I/ref:8.4f}   max rho = {data[lab].max():.4g}")
lines.append("  (the old AV8 table's convention is unknown -> its ratio is meaningless)")

# pD / full ratios per potential
lines.append("")
lines.append("=== pD-only fraction of the full table (same potential) ===")
for full, pd in (("AV18 (u=1)", "pD AV18 (u=6)"), ("AV18+UIX (u=2)", "pD AV18+UIX (u=7)"),
                 ("CD-Bonn (u=4)", "pD CD-Bonn (u=9)")):
    if full in data and pd in data:
        lines.append(f"  {pd:20s} / {full:20s} = {measure_integral(data[pd])/measure_integral(data[full]):.4f}")

# ---- marginals ----
fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
for lab, fn, col, ls in TABLES:
    if lab not in data:
        continue
    rho = data[lab]
    w_th = (kcm**2)[None, :, None]*(krel**2)[None, None, :]
    m_th = np.trapezoid(np.trapezoid(rho*w_th, krel, axis=2), kcm, axis=1)*np.sin(np.radians(th))
    w_c = np.sin(np.radians(th))[:, None, None]*(krel**2)[None, None, :]
    m_c = np.trapezoid(np.trapezoid(rho*w_c, krel, axis=2), np.radians(th), axis=0)*kcm**2
    w_r = np.sin(np.radians(th))[:, None, None]*(kcm**2)[None, :, None]
    m_r = np.trapezoid(np.trapezoid(rho*w_r, kcm, axis=1), np.radians(th), axis=0)*krel**2
    for a, x, m in ((ax[0], th, m_th), (ax[1], kcm, m_c), (ax[2], krel, m_r)):
        s = m/np.trapezoid(m, x)
        a.plot(x, s, color=col, ls=ls, lw=1.8, label=lab)
ax[0].set_xlabel(r"$\theta_{ab}(p_a,p_b)$ [deg]"); ax[1].set_xlabel(r"$p_a = p_1 - P_{cm}/3$ [fm$^{-1}$]")
ax[2].set_xlabel(r"$p_b = (p_2-p_3)/2$ [fm$^{-1}$]")
ax[1].set_xlim(0, 5); ax[2].set_xlim(0, 5)
for a in ax: a.set_ylabel("normalised marginal"); a.set_ylim(bottom=0)
ax[0].legend(frameon=False, fontsize=8)
fig.suptitle("3N density tables: shape comparison (each marginal area-normalised)", fontsize=12)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "tables_marginals.png"), dpi=170); plt.close(fig)

# log-scale k_rel tail, the SRC-relevant region
fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))
for lab, fn, col, ls in TABLES:
    if lab not in data: continue
    rho = data[lab]
    w_r = np.sin(np.radians(th))[:, None, None]*(kcm**2)[None, :, None]
    m_r = np.trapezoid(np.trapezoid(rho*w_r, kcm, axis=1), np.radians(th), axis=0)*krel**2
    w_c = np.sin(np.radians(th))[:, None, None]*(krel**2)[None, None, :]
    m_c = np.trapezoid(np.trapezoid(rho*w_c, krel, axis=2), np.radians(th), axis=0)*kcm**2
    ax[0].semilogy(krel, m_r/np.trapezoid(m_r, krel), color=col, ls=ls, lw=1.8, label=lab)
    ax[1].semilogy(kcm, m_c/np.trapezoid(m_c, kcm), color=col, ls=ls, lw=1.8)
ax[0].set_xlabel(r"$p_b = (p_2-p_3)/2$ [fm$^{-1}$]"); ax[1].set_xlabel(r"$p_a = p_1 - P_{cm}/3$ [fm$^{-1}$]")
ax[0].set_ylabel("normalised marginal"); ax[0].legend(frameon=False, fontsize=8)
for a in ax: a.set_xlim(0, 8); a.set_ylim(1e-8, None)
fig.suptitle("High-momentum tails (log scale)", fontsize=12)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "tables_tails.png"), dpi=170); plt.close(fig)

# mean values in the SRC window
lines.append("")
lines.append("=== k_rel marginal: m(k_rel) = k_rel^2 * int rho k_cm^2 sin(th) dk_cm dth ===")
lines.append("    (Jacobi convention: table k_cm = p_a = p1 - P_cm/3, k_rel = p_b = (p2-p3)/2,")
lines.append("     theta = angle(p_a,p_b); the measure above is the Jacobi phase space in get_rho())")
lines.append("    thresholds: 1.115 fm^-1 = 0.22 GeV/c (the kF used in the runs), 1.27 fm^-1 = 0.25 GeV/c")
sel = krel > 1.27
sel22 = krel > 1.115
for lab in data:
    rho = data[lab]
    w_r = np.sin(np.radians(th))[:, None, None]*(kcm**2)[None, :, None]
    m_r = np.trapezoid(np.trapezoid(rho*w_r, kcm, axis=1), np.radians(th), axis=0)*krel**2
    frac = np.trapezoid(m_r[sel], krel[sel])/np.trapezoid(m_r, krel)
    frac22 = np.trapezoid(m_r[sel22], krel[sel22])/np.trapezoid(m_r, krel)
    mean = np.trapezoid(m_r*krel, krel)/np.trapezoid(m_r, krel)
    w_th = (kcm**2)[None, :, None]*(krel**2)[None, None, :]
    m_th = np.trapezoid(np.trapezoid(rho*w_th, krel, axis=2), kcm, axis=1)*np.sin(np.radians(th))
    mth = np.trapezoid(m_th*th, th)/np.trapezoid(m_th, th)
    lines.append(f"  {lab:20s} frac(p_b>0.22GeV) = {100*frac22:6.3f} %   frac(p_b>0.25GeV) = {100*frac:6.3f} %"
                 f"   <p_b> = {mean:5.3f} fm^-1   <theta_ab> = {mth:5.1f} deg")

open(os.path.join(OUT, "tables_summary.txt"), "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
print("\nsaved tables_marginals.png, tables_tails.png, tables_summary.txt in", OUT)
