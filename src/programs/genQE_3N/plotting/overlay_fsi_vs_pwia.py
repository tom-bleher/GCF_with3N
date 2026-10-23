"""Overlay of the 12C/4He alpha3N ratio with and without FSI on 12C.

Inputs default to the v3 samples (post-Pauli-veto-removal, internal-momentum
sampling cut). Run from anywhere; paths are resolved relative to this file.
"""
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_alpha3N_ratio import load_alpha3N

ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
EV = os.path.join(ROOT_DIR, 'events/new_fsi_v3')
OUT = os.path.join(ROOT_DIR, 'analysis/Plots/alpha3N_ratio_FSI_vs_PWIA_overlay_v3')
EBEAM = 10.6
# w_max trims the extreme importance-sampling weight tail (a few events per
# 100M at raw w ~ 1e-8; one 4He event at 5e-8 dominated its alpha bin).
kw = dict(Q2_min=4.1, Q2_max=4.3, w_max=1e-8)

a_fsi, w_fsi, _ = load_alpha3N(f'{EV}/fsi_hN_106.root', EBEAM, **kw)
a_pw,  w_pw,  _ = load_alpha3N(f'{EV}/pwia_12C_106.root', EBEAM, **kw)
a_he,  w_he,  _ = load_alpha3N(f'{EV}/pwia_4He_106.root', EBEAM, **kw)

# Auto-range the alpha3N axis to the populated region (the kinematic cuts,
# e.g. a narrow Q2 slice at fixed theta_e, can move it far from [0.6, 2.6]).
a_all = np.concatenate([a_fsi, a_pw, a_he])
lo, hi = np.percentile(a_all, [0.02, 99.98])
pad = 0.02 * (hi - lo)
edges = np.linspace(lo - pad, hi + pad, 41)
centers = 0.5 * (edges[:-1] + edges[1:])


def ratio_with_err(a_num, w_num, a_den, w_den, min_n=20):
    H_n, _ = np.histogram(a_num, bins=edges, weights=w_num)
    H_d, _ = np.histogram(a_den, bins=edges, weights=w_den)
    V_n, _ = np.histogram(a_num, bins=edges, weights=w_num ** 2)
    V_d, _ = np.histogram(a_den, bins=edges, weights=w_den ** 2)
    N_n, _ = np.histogram(a_num, bins=edges)
    N_d, _ = np.histogram(a_den, bins=edges)
    good = (N_n >= min_n) & (N_d >= min_n) & (H_d > 0)
    r = np.where(good, H_n / np.where(H_d > 0, H_d, 1.0), np.nan)
    rel = np.where(good, np.sqrt(np.where(H_n > 0, V_n / H_n ** 2, 0.0)
                                 + np.where(H_d > 0, V_d / H_d ** 2, 0.0)), np.nan)
    return r, r * rel


r_fsi, e_fsi = ratio_with_err(a_fsi, w_fsi, a_he, w_he)
r_pw,  e_pw  = ratio_with_err(a_pw,  w_pw,  a_he, w_he)

fig, (ax, axd) = plt.subplots(2, 1, figsize=(7, 6.5), sharex=True,
                              gridspec_kw={'height_ratios': [1.4, 1]})
ax.errorbar(centers, r_fsi, yerr=e_fsi, fmt='o', ms=4, lw=0.8, capsize=2,
            color='C0', label=r'$^{12}$C 3N+FSI hN / $^{4}$He PWIA')
ax.errorbar(centers, r_pw, yerr=e_pw, fmt='s', ms=4, lw=0.8, capsize=2,
            color='C3', label=r'$^{12}$C 3N PWIA / $^{4}$He PWIA')
ax.set_ylabel(r'$\sigma(^{12}\mathrm{C})\,/\,\sigma(^{4}\mathrm{He})$')
ax.legend(frameon=False, fontsize=10)
ax.grid(alpha=0.3, lw=0.4)
ax.set_title(r'$E_{beam}=10.6$ GeV, $\theta_e\in[11°,12°]$, $Q^2\in[4.1,4.3]$ GeV$^2$')

# (12C FSI / 4He) / (12C PWIA / 4He) = 12C FSI / 12C PWIA: the shared 4He
# histogram cancels, and so does its error. Take the double ratio from the two
# (independently generated) 12C histograms directly, on the bins where both
# single ratios are shown.
dbl, dbl_err = ratio_with_err(a_fsi, w_fsi, a_pw, w_pw)
shown = np.isfinite(r_fsi) & np.isfinite(r_pw)
dbl = np.where(shown, dbl, np.nan)
dbl_err = np.where(shown, dbl_err, np.nan)
axd.errorbar(centers, dbl, yerr=dbl_err, fmt='o', ms=4, lw=0.8, capsize=2, color='k')
axd.axhline(1.0, color='gray', lw=0.8, ls='--')
axd.set_xlabel(r'$\alpha_{3N}$')
axd.set_ylabel('with FSI / without FSI')
axd.set_xlim(edges[0], edges[-1])
axd.grid(alpha=0.3, lw=0.4)

fig.tight_layout()
fig.savefig(OUT + '.pdf', dpi=300)
fig.savefig(OUT + '.png', dpi=200)
print('Saved', OUT + '.png')
