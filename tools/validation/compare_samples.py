#!/usr/bin/env python3
"""Compare two event samples statistically: the split pipeline vs the combined
generators (or any two samples of the same physics).

Reads either format:
  - the GCF event-file format (tree "events", from gcf_generate/gcf_fsi);
  - the combined generators' output (tree "genT", from genQE_3N_FSI or
    genQE_FSI): electron pe[3], struck nucleon p1_gen[3] (3N) or pLead[3]
    (2N), and the same part_* final-state list.

For each observable it does an unweighted two-sample chi^2 test (FSI does
not depend on the event weight, so the unweighted test has the full
statistical power of the sample) and reports the weighted cross section
per generation attempt of each sample. Exits 1 if any p-value is below
--alpha, so it can serve as a pass/fail flag.

  compare_samples.py A.root B.root [--attempts-a N] [--attempts-b N] [-o plot.png]
"""

import argparse
import sys

import awkward as ak
import numpy as np
import uproot
from scipy import stats

M_N = 0.93892
PROTON, NEUTRON = 2212, 2112


def load(path, attempts=None, beam=None):
    f = uproot.open(path)
    keys = {k.split(";")[0] for k in f.keys()}
    if "events" in keys:
        t = f["events"]
        a = t.arrays(["weight", "beam_E", "e_px", "e_py", "e_pz", "e_E", "nuc_role",
                      "nuc_px", "nuc_py", "nuc_pz", "part_pdg", "part_px", "part_py", "part_pz",
                      "part_rescatterCode"])
        beamE = np.asarray(a["beam_E"])
        e = np.stack([np.asarray(a[k]) for k in ("e_px", "e_py", "e_pz", "e_E")])
        struck = a["nuc_role"] == 0
        s = np.stack([np.asarray(ak.firsts(a[k][struck])) for k in ("nuc_px", "nuc_py", "nuc_pz")])
        fmt = "gcf"
    else:
        t = f["genT"]
        lead = "p1_gen" if "p1_gen" in t.keys() else "pLead"
        a = t.arrays(["weight", "pe", lead, "part_pdg", "part_px", "part_py", "part_pz",
                      "part_rescatterCode"])
        if beam is None:
            raise SystemExit(f"{path}: combined-generator file, give the beam energy with --beam")
        beamE = np.full(len(a), beam)
        pe = np.asarray(a["pe"]).T
        e = np.vstack([pe, np.sqrt((pe**2).sum(0))])
        s = np.asarray(a[lead]).T
        fmt = f"genT ({lead})"
    if attempts is None:
        if "hAttempts" not in keys:
            raise SystemExit(f"{path}: no hAttempts; give the number of attempts")
        attempts = float(f["hAttempts"].values()[0])
    w = np.asarray(a["weight"])

    q = np.stack([-e[0], -e[1], beamE - e[2]])
    nu = beamE - e[3]
    Q2 = (q**2).sum(0) - nu**2
    pdg = a["part_pdg"]
    p = np.sqrt(a["part_px"]**2 + a["part_py"]**2 + a["part_pz"]**2)
    pp = p[pdg == PROTON]
    fastest = ak.argmax(pp, axis=1, keepdims=True)
    pz_p = a["part_pz"][pdg == PROTON]
    p_fast = np.asarray(ak.fill_none(ak.firsts(pp[fastest]), np.nan))
    cos_fast = np.asarray(ak.fill_none(ak.firsts(pz_p[fastest]), np.nan)) / p_fast
    untouched = ak.all((a["part_rescatterCode"] == 1) | (a["part_rescatterCode"] == -1), axis=1)
    obs = {
        "Q2 [GeV^2]": (Q2, np.linspace(0, 6, 41)),
        "x_B": (Q2 / (2 * M_N * nu), np.linspace(0, 3, 41)),
        "struck nucleon |p|, generated [GeV/c]": (np.sqrt((s**2).sum(0)), np.linspace(0, 7, 41)),
        "protons after FSI": (np.asarray(ak.sum(pdg == PROTON, axis=1)), np.arange(-0.5, 12.5)),
        "neutrons after FSI": (np.asarray(ak.sum(pdg == NEUTRON, axis=1)), np.arange(-0.5, 12.5)),
        "charged pions after FSI": (np.asarray(ak.sum(abs(pdg) == 211, axis=1)), np.arange(-0.5, 5.5)),
        "fastest proton |p| [GeV/c]": (p_fast, np.linspace(0, 6, 41)),
        "fastest proton cos(theta)": (cos_fast, np.linspace(-1, 1, 41)),
        "no FSI interaction (all nucleons fate 1)": (np.asarray(untouched).astype(float), np.array([-0.5, 0.5, 1.5])),
    }
    return dict(path=path, fmt=fmt, n=len(w), w=w, attempts=attempts, obs=obs)


def chi2_unweighted(x1, x2, bins):
    n1, _ = np.histogram(x1[np.isfinite(x1)], bins=bins)
    n2, _ = np.histogram(x2[np.isfinite(x2)], bins=bins)
    N1, N2 = n1.sum(), n2.sum()
    use = (n1 + n2) > 0
    chi2 = np.sum((N2 * n1[use] - N1 * n2[use])**2 / (N1 * N2 * (n1[use] + n2[use])))
    ndf = max(int(use.sum()) - 1, 1)
    return chi2, ndf, stats.chi2.sf(chi2, ndf), n1, n2


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("a")
    ap.add_argument("b")
    ap.add_argument("--attempts-a", type=float)
    ap.add_argument("--attempts-b", type=float)
    ap.add_argument("--beam", type=float, help="beam energy for genT files (they don't store it)")
    ap.add_argument("--alpha", type=float, default=1e-3, help="fail if any p-value is below this")
    ap.add_argument("-o", "--output", help="PNG with the unweighted overlays")
    args = ap.parse_args()

    A = load(args.a, args.attempts_a, args.beam)
    B = load(args.b, args.attempts_b, args.beam)
    print(f"A: {A['path']} [{A['fmt']}], {A['n']} events, {A['attempts']:.0f} attempts")
    print(f"B: {B['path']} [{B['fmt']}], {B['n']} events, {B['attempts']:.0f} attempts\n")

    def xs(S):
        return S["w"].sum() / S["attempts"], np.sqrt((S["w"]**2).sum()) / S["attempts"]
    (xa, ea), (xb, eb) = xs(A), xs(B)
    pull = (xa - xb) / np.hypot(ea, eb)
    print(f"sum(w)/attempts:  A = {xa:.5g} +- {ea:.2g}   B = {xb:.5g} +- {eb:.2g}   "
          f"pull = {pull:+.2f}  (N_eff A = {A['w'].sum()**2 / (A['w']**2).sum():.0f}, "
          f"B = {B['w'].sum()**2 / (B['w']**2).sum():.0f})\n")
    print(f"{'observable (unweighted)':44s} {'chi2/ndf':>12s} {'p-value':>9s}  mean A / mean B")
    worst = 1.0
    results = []
    for name in A["obs"]:
        xa_, bins = A["obs"][name]
        xb_, _ = B["obs"][name]
        chi2, ndf, p, n1, n2 = chi2_unweighted(xa_, xb_, bins)
        worst = min(worst, p)
        results.append((name, bins, n1, n2))
        print(f"{name:44s} {chi2:7.1f}/{ndf:<4d} {p:9.3f}  {np.nanmean(xa_):.4f} / {np.nanmean(xb_):.4f}")
    ok = worst >= args.alpha and abs(pull) < 4
    print(f"\n{'PASS' if ok else 'FAIL'}: smallest p-value {worst:.3g} (alpha {args.alpha}), "
          f"cross-section pull {pull:+.2f}")

    if args.output:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        cols = 3
        rows = (len(results) + cols - 1) // cols
        fig, axes = plt.subplots(rows, cols, figsize=(4.4 * cols, 3.2 * rows), squeeze=False)
        for ax, (name, bins, n1, n2) in zip(axes.flat, results):
            centers = 0.5 * (bins[1:] + bins[:-1])
            for n, lab in ((n1, "A"), (n2, "B")):
                ax.errorbar(centers, n / n.sum(), yerr=np.sqrt(n) / n.sum(), fmt="o" if lab == "A" else "s",
                            ms=3, label=lab, alpha=0.8)
            ax.set_xlabel(name)
            ax.legend(fontsize=8)
        for ax in list(axes.flat)[len(results):]:
            ax.axis("off")
        fig.suptitle(f"A: {A['path']}   B: {B['path']}   (unweighted, unit area)", fontsize=9)
        fig.tight_layout()
        fig.savefig(args.output, dpi=100)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
