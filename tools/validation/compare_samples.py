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

An unweighted test can't see a change that only moves weights, such as a
different contact: the events are the same and only their weights differ.
--weighted decides on a weighted test instead. Its statistic is Pearson's chi^2 of the
normalized weighted histograms, with the pooled shape and each sample's effective number
of events in the denominator; each bin's own sum of w^2 there made a bin that came out low
look significant for that reason alone. The weights in this project are too heavy-tailed
for that statistic to follow a chi^2 distribution, and the samples can have different
weight distributions (say, with and without importance sampling), so swapping events
between them isn't valid either. Its p-value comes from a centred bootstrap: each
sample is resampled with its own weights, and the difference of each replica from its
own sample is compared with the difference between the samples. This tests the shapes;
the difference in scale is the pull of sum(w)/attempts. It starts with BOOTSTRAPS
replicas and goes to BOOTSTRAPS_MAX for a p-value below 0.05, so its smallest p-value
is about 1/BOOTSTRAPS_MAX. Samples made with the same seed share events, which makes
this test conservative.

Files without FSI (no part_* branches) get the generated observables only.

  compare_samples.py A.root B.root [--attempts-a N] [--attempts-b N] [-o plot.png]
                     [--weighted] [--max-events N] [--json OUT.json]
"""

import argparse
import json
import sys

import awkward as ak
import numpy as np
import uproot
from scipy import stats

M_N = 0.93892
PROTON, NEUTRON = 2212, 2112
BOOTSTRAPS = 200
BOOTSTRAPS_MAX = 1000
FSI_BRANCHES = ["part_pdg", "part_px", "part_py", "part_pz", "part_rescatterCode"]


def load(path, attempts=None, beam=None, max_events=None):
    """A sample's observables, from at most max_events events. Its weight sum and
    attempts always come from the whole file, so the normalization is unaffected."""
    f = uproot.open(path)
    keys = {k.split(";")[0] for k in f.keys()}
    if "events" in keys:
        t = f["events"]
        fsi = set(FSI_BRANCHES) <= set(t.keys())
        residual = [b for b in ("res_px", "res_py", "res_pz") if b in t.keys()]
        a = t.arrays(["weight", "beam_E", "e_px", "e_py", "e_pz", "e_E", "nuc_role",
                      "nuc_px", "nuc_py", "nuc_pz"] + residual + (FSI_BRANCHES if fsi else []), entry_stop=max_events)
        beamE = np.asarray(a["beam_E"])
        e = np.stack([np.asarray(a[k]) for k in ("e_px", "e_py", "e_pz", "e_E")])
        struck = a["nuc_role"] == 0
        s = np.stack([np.asarray(ak.firsts(a[k][struck])) for k in ("nuc_px", "nuc_py", "nuc_pz")])
        partner = a["nuc_role"] == 1
        recoil = np.stack([np.asarray(ak.fill_none(ak.firsts(a[k][partner]), np.nan), dtype=float)
                           for k in ("nuc_px", "nuc_py", "nuc_pz")])
        res = np.stack([np.asarray(a[k], dtype=float) for k in residual]) if len(residual) == 3 else None
        fmt = "gcf"
    else:
        t = f["genT"]
        fsi = True
        lead = "p1_gen" if "p1_gen" in t.keys() else "pLead"
        a = t.arrays(["weight", "pe", lead, "part_pdg", "part_px", "part_py", "part_pz",
                      "part_rescatterCode"], entry_stop=max_events)
        if beam is None:
            raise SystemExit(f"{path}: combined-generator file, give the beam energy with --beam")
        beamE = np.full(len(a), beam)
        pe = np.asarray(a["pe"]).T
        e = np.vstack([pe, np.sqrt((pe**2).sum(0))])
        s = np.asarray(a[lead]).T
        recoil = res = None
        fmt = f"genT ({lead})"
    if attempts is None:
        if "hAttempts" not in keys:
            raise SystemExit(f"{path}: no hAttempts; give the number of attempts")
        attempts = float(f["hAttempts"].values()[0])
    w = np.asarray(a["weight"], dtype=float)
    w_all = t["weight"].array(library="np").astype(float) if max_events and t.num_entries > max_events else w

    q = np.stack([-e[0], -e[1], beamE - e[2]])
    nu = beamE - e[3]
    Q2 = (q**2).sum(0) - nu**2
    obs = {
        "Q2 [GeV^2]": (Q2, np.linspace(0, 6, 41)),
        "x_B": (Q2 / (2 * M_N * nu), np.linspace(0, 3, 41)),
        "struck nucleon |p|, generated [GeV/c]": (np.sqrt((s**2).sum(0)), np.linspace(0, 7, 41)),
    }
    if recoil is not None:
        # The observables that respond to the nuclear inputs (sigma_CM, E*, contacts).
        p1 = s - q
        obs["p_miss = p_struck - q [GeV/c]"] = (np.sqrt((p1**2).sum(0)), np.linspace(0, 1.5, 31))
        obs["partner |p|, generated [GeV/c]"] = (np.sqrt((recoil**2).sum(0)), np.linspace(0, 1.5, 31))
        with np.errstate(invalid="ignore", divide="ignore"):
            cos_open = (p1 * recoil).sum(0) / np.sqrt((p1**2).sum(0) * (recoil**2).sum(0))
        obs["cos(p_miss, partner)"] = (cos_open, np.linspace(-1, 1, 41))
    if res is not None:
        obs["residual nucleus |p| [GeV/c]"] = (np.sqrt((res**2).sum(0)), np.linspace(0, 1, 31))
    if fsi:
        pdg = a["part_pdg"]
        p = np.sqrt(a["part_px"]**2 + a["part_py"]**2 + a["part_pz"]**2)
        pp = p[pdg == PROTON]
        fastest = ak.argmax(pp, axis=1, keepdims=True)
        pz_p = a["part_pz"][pdg == PROTON]
        p_fast = np.asarray(ak.fill_none(ak.firsts(pp[fastest]), np.nan))
        cos_fast = np.asarray(ak.fill_none(ak.firsts(pz_p[fastest]), np.nan)) / p_fast
        untouched = ak.all((a["part_rescatterCode"] == 1) | (a["part_rescatterCode"] == -1), axis=1)
        obs.update({
            "protons after FSI": (np.asarray(ak.sum(pdg == PROTON, axis=1)), np.arange(-0.5, 12.5)),
            "neutrons after FSI": (np.asarray(ak.sum(pdg == NEUTRON, axis=1)), np.arange(-0.5, 12.5)),
            "charged pions after FSI": (np.asarray(ak.sum(abs(pdg) == 211, axis=1)), np.arange(-0.5, 5.5)),
            "fastest proton |p| [GeV/c]": (p_fast, np.linspace(0, 6, 41)),
            "fastest proton cos(theta)": (cos_fast, np.linspace(-1, 1, 41)),
            "no FSI interaction (all nucleons fate 1)": (np.asarray(untouched).astype(float), np.array([-0.5, 0.5, 1.5])),
        })
    return dict(path=path, fmt=fmt, n=t.num_entries, used=len(w), w=w, w_all=w_all, attempts=attempts, obs=obs)


def chi2_unweighted(x1, x2, bins):
    n1, _ = np.histogram(x1[np.isfinite(x1)], bins=bins)
    n2, _ = np.histogram(x2[np.isfinite(x2)], bins=bins)
    N1, N2 = n1.sum(), n2.sum()
    use = (n1 + n2) > 0
    chi2 = np.sum((N2 * n1[use] - N1 * n2[use])**2 / (N1 * N2 * (n1[use] + n2[use])))
    ndf = max(int(use.sum()) - 1, 1)
    return chi2, ndf, stats.chi2.sf(chi2, ndf), n1, n2


def weighted_hist(x, w, bins):
    ok = np.isfinite(x)
    return (np.histogram(x[ok], bins=bins, weights=w[ok])[0], np.histogram(x[ok], bins=bins, weights=w[ok]**2)[0],
            np.histogram(x[ok], bins=bins)[0])


def chi2_weighted(x1, w1, x2, w2, bins, rng=None):
    """Whether two weighted samples have the same shape: Pearson's chi^2 of their normalized
    weighted histograms, with a p-value from a centred bootstrap. Returns (chi2, number
    of bins in use - 1, p); NaNs if there are too few bins."""
    rng = rng or np.random.default_rng(20261104)

    def prepare(x, w):
        ok = np.isfinite(x)
        x, w = x[ok], w[ok]
        nb = len(bins) - 1
        idx = np.searchsorted(bins, x, side="right") - 1
        idx[x == bins[-1]] = nb - 1
        inside = (idx >= 0) & (idx < nb)
        return idx[inside], w[inside]

    nb = len(bins) - 1
    idx1, w1 = prepare(x1, w1)
    idx2, w2 = prepare(x2, w2)
    if not len(idx1) or not len(idx2) or w1.sum() <= 0 or w2.sum() <= 0:
        return float("nan"), 0, float("nan")

    def shape(idx, w):
        """The normalized histogram of a sample and its effective number of events."""
        total = np.bincount(idx, w, nb)
        return total / total.sum(), w.sum()**2 / (w * w).sum()

    p1, n1 = shape(idx1, w1)
    p2, n2 = shape(idx2, w2)
    # Pearson's denominator, of the pooled shape: the same for the data and its replicas,
    # and hardly moved by a fluctuation in the bin. A bin's own sum of w^2 is neither, and
    # a bin that came out low then looks significant because it came out low.
    pooled = (n1 * p1 + n2 * p2) / (n1 + n2)
    use = pooled > 0
    if use.sum() < 2:
        return float("nan"), 0, float("nan")
    scale = 1 / (1 / n1 + 1 / n2)

    def statistic(d):
        return float(scale * np.sum(d[use]**2 / pooled[use]))

    observed = statistic(p1 - p2)

    def replica(idx, w, centre):
        # A Poisson bootstrap of one sample, with its own weights.
        total = np.bincount(idx, w * rng.poisson(1.0, len(w)), nb)
        norm = total.sum()
        return total / norm - centre if norm > 0 else np.zeros(nb)

    above = tried = 0
    while True:
        target = BOOTSTRAPS if tried < BOOTSTRAPS else BOOTSTRAPS_MAX
        for _ in range(target - tried):
            above += statistic(replica(idx1, w1, p1) - replica(idx2, w2, p2)) >= observed
        tried = target
        p = (1 + above) / (1 + tried)
        if tried >= BOOTSTRAPS_MAX or p >= 0.05:
            return observed, int(use.sum()) - 1, p


def sample_info(S):
    w = S["w_all"]
    sum_w, sum_w2 = float(w.sum()), float((w**2).sum())
    return dict(path=str(S["path"]), format=S["fmt"], events=int(S["n"]), events_used=int(S["used"]),
                attempts=float(S["attempts"]), yield_=sum_w / S["attempts"], yield_err=sum_w2**0.5 / S["attempts"],
                n_eff=sum_w**2 / sum_w2 if sum_w2 > 0 else 0.)


def plot(results, A, B, output, weighted):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    cols = 3
    rows = (len(results) + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(4.4 * cols, 3.2 * rows), squeeze=False)
    for ax, r in zip(axes.flat, results):
        name, bins = r["name"], r["bins"]
        centers = 0.5 * (bins[1:] + bins[:-1])
        for S, lab, marker in ((A, "A", "o"), (B, "B", "s")):
            x, _ = S["obs"][name]
            content, variance, _ = weighted_hist(x, S["w"] if weighted else np.ones_like(S["w"]), bins)
            total = content.sum()
            if total > 0:
                ax.errorbar(centers, content / total, yerr=np.sqrt(variance) / total, fmt=marker, ms=3, label=lab, alpha=0.8)
        ax.set_xlabel(name)
        ax.legend(fontsize=8)
    for ax in list(axes.flat)[len(results):]:
        ax.axis("off")
    kind = "weighted" if weighted else "unweighted"
    fig.suptitle(f"A: {A['path']}   B: {B['path']}   ({kind}, unit area)", fontsize=9)
    fig.tight_layout()
    fig.savefig(output, dpi=100)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("a")
    ap.add_argument("b")
    ap.add_argument("--attempts-a", type=float)
    ap.add_argument("--attempts-b", type=float)
    ap.add_argument("--beam", type=float, help="beam energy for genT files (they don't store it)")
    ap.add_argument("--alpha", type=float, default=1e-3, help="fail if any p-value is below this")
    ap.add_argument("--weighted", action="store_true", help="decide on the weighted test and plot weighted histograms")
    ap.add_argument("--max-events", type=int, help="use at most this many events of each file for the observables")
    ap.add_argument("--json", help="write the results to this file")
    ap.add_argument("-o", "--output", help="PNG with the overlays")
    args = ap.parse_args()

    A = load(args.a, args.attempts_a, args.beam, args.max_events)
    B = load(args.b, args.attempts_b, args.beam, args.max_events)
    print(f"A: {A['path']} [{A['fmt']}], {A['n']} events, {A['attempts']:.0f} attempts")
    print(f"B: {B['path']} [{B['fmt']}], {B['n']} events, {B['attempts']:.0f} attempts\n")

    info_a, info_b = sample_info(A), sample_info(B)
    xa, ea, xb, eb = info_a["yield_"], info_a["yield_err"], info_b["yield_"], info_b["yield_err"]
    pull = (xa - xb) / np.hypot(ea, eb)
    print(f"sum(w)/attempts:  A = {xa:.5g} +- {ea:.2g}   B = {xb:.5g} +- {eb:.2g}   "
          f"pull = {pull:+.2f}  (N_eff A = {info_a['n_eff']:.0f}, B = {info_b['n_eff']:.0f})\n")
    print(f"{'observable':44s} {'unweighted':>12s} {'p':>9s}  {'weighted':>12s} {'p':>9s}  mean A / mean B")
    worst = 1.0
    results = []
    for name in A["obs"]:
        if name not in B["obs"]:
            continue           # only one of the files has FSI
        xa_, bins = A["obs"][name]
        xb_, _ = B["obs"][name]
        chi2, ndf, p, _, _ = chi2_unweighted(xa_, xb_, bins)
        chi2w, ndfw, pw = chi2_weighted(xa_, A["w"], xb_, B["w"], bins)
        decide = pw if args.weighted else p
        if np.isfinite(decide):
            worst = min(worst, decide)
        results.append(dict(name=name, bins=bins, chi2=float(chi2), ndf=ndf, p=float(p),
                            chi2_weighted=chi2w, ndf_weighted=ndfw, p_weighted=pw,
                            mean_a=float(np.nanmean(xa_)), mean_b=float(np.nanmean(xb_))))
        print(f"{name:44s} {chi2:7.1f}/{ndf:<4d} {p:9.3f}  {chi2w:7.1f}/{ndfw:<4d} {pw:9.3f}  "
              f"{np.nanmean(xa_):.4f} / {np.nanmean(xb_):.4f}")
    ok = worst >= args.alpha and abs(pull) < 4
    print(f"\n{'PASS' if ok else 'FAIL'}: smallest {'weighted ' if args.weighted else ''}p-value {worst:.3g} "
          f"(alpha {args.alpha}), cross-section pull {pull:+.2f}")

    if args.json:
        out = dict(a=info_a, b=info_b, pull=float(pull), alpha=args.alpha, weighted=args.weighted,
                   smallest_p=worst, ok=bool(ok),
                   observables=[{k: v for k, v in r.items() if k != "bins"} for r in results])
        with open(args.json, "w") as f:
            # NaN is not JSON: a test with too few bins has no result.
            json.dump(json.loads(json.dumps(out), parse_constant=lambda _: None), f, indent=1)
    if args.output:
        plot(results, A, B, args.output, args.weighted)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
