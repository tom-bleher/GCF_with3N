#!/usr/bin/env python3
"""(1) Is the old AV8 table the same as the new Faddeev AV18 up to a constant?
   (2) What exactly is the new tables' normalisation?

Grid: theta 61 x [0,180] deg, k_a 101 x [0,10] fm^-1, k_b 101 x [0,10] fm^-1,
one value per line, ordered theta -> k_a -> k_b (QEGenerator_3N::fill_array).

Full Jacobi measure, for f depending only on |p_a|, |p_b| and the angle between:
  int d3p_a d3p_b f = 4pi * 2pi * int f p_a^2 p_b^2 sin(theta) dp_a dp_b dtheta
                    = 8 pi^2 * A
"""
import numpy as np
ARR = "/home/ilayws/Ilay-Generators/src/programs/genQE_3N/array"
NT, NK = 61, 101
HBARC = 0.1973269804          # GeV*fm
th = np.radians(np.linspace(0, 180, NT))
k = np.linspace(0, 10, NK)     # fm^-1
kG = k * HBARC                 # same grid expressed in GeV/c


def load(fn):
    return np.loadtxt(f"{ARR}/{fn}").reshape(NT, NK, NK)


def A_of(rho, kk):
    """int rho * k_a^2 k_b^2 sin(theta) dk_a dk_b dtheta, k axis in the given units"""
    w = (kk**2)[None, :, None] * (kk**2)[None, None, :] * np.sin(th)[:, None, None]
    return np.trapezoid(np.trapezoid(np.trapezoid(rho*w, kk, axis=2), kk, axis=1), th)


av18 = load("he3_3Ndensity_AV18.dnst")
av8 = load("AV8_1D.dnst")
EIGHT_PI2 = 8*np.pi**2

print("=== (2) normalisation ===")
for nm, rho in (("Faddeev AV18 (new)", av18), ("AV8 (old)", av8)):
    A_fm, A_GeV = A_of(rho, k), A_of(rho, kG)
    print(f"{nm}")
    print(f"    int rho d3pa d3pb, momenta in fm^-1  = 8pi^2 * {A_fm:.6g} = {EIGHT_PI2*A_fm:.6g}")
    print(f"    int rho d3pa d3pb, momenta in GeV/c  = 8pi^2 * {A_GeV:.6g} = {EIGHT_PI2*A_GeV:.6g}")
    print(f"      ... / (1/3) = {EIGHT_PI2*A_GeV*3:.6g} (GeV)   {EIGHT_PI2*A_fm*3:.6g} (fm^-1)")
    print(f"      ... / (2pi/3) = {EIGHT_PI2*A_fm/(2*np.pi/3):.6g} (fm^-1)")

print()
print("=== (1) AV8 vs AV18: same shape up to a constant? ===")
both = (av18 > 0) & (av8 > 0)
print(f"  grid points with both > 0: {both.sum()} / {av18.size}"
      f"   (AV18>0: {(av18>0).sum()}, AV8>0: {(av8>0).sum()})")
r = av8[both] / av18[both]
w = av18[both]                      # density-weighted, so the bulk dominates
print(f"  ratio AV8/AV18:  median {np.median(r):.6g}   mean {r.mean():.6g}")
print(f"    percentiles  1% {np.percentile(r,1):.4g}   25% {np.percentile(r,25):.4g}"
      f"   50% {np.median(r):.4g}   75% {np.percentile(r,75):.4g}   99% {np.percentile(r,99):.4g}")
print(f"    min {r.min():.4g}   max {r.max():.4g}   spread(max/min) = {r.max()/r.min():.4g}")
wm = np.average(r, weights=w)
print(f"  density-weighted mean ratio {wm:.6g}, weighted rms/mean = "
      f"{np.sqrt(np.average((r-wm)**2, weights=w))/wm:.4g}")
# best single constant and the residual after removing it
c = np.sum(w*av8[both])/np.sum(w*av18[both])
res = av8[both] - c*av18[both]
print(f"  best-fit constant c (density-weighted LS): {c:.6g}")
print(f"  residual |AV8 - c*AV18| / AV8, density-weighted mean = "
      f"{np.average(np.abs(res)/np.maximum(av8[both],1e-300), weights=w):.4g}")
print(f"  correlation of log densities: {np.corrcoef(np.log(av18[both]), np.log(av8[both]))[0,1]:.6f}")
# where do they differ most (by region of the grid)?
print()
print("  ratio AV8/AV18 by k_b (pp half-relative momentum) band, density-weighted:")
for lo, hi in ((0, 1), (1, 2), (2, 3), (3, 5), (5, 10)):
    m = both & ((k >= lo) & (k < hi))[None, None, :]
    if m.sum():
        rr = av8[m]/av18[m]; ww = av18[m]
        print(f"    k_b in [{lo:2d},{hi:2d}) fm^-1 : weighted mean ratio {np.average(rr, weights=ww):10.4g}"
              f"   median {np.median(rr):10.4g}")
