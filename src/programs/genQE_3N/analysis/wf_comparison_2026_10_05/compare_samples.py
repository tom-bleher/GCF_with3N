#!/usr/bin/env python3
"""Compare generated 3N samples across wavefunction tables.

Samples (events/wf_test/): pwia_u{1,4,6,9}.root, fsi_u{1,4,6,9}.root and, if
present, the old-code AV8 baseline. u=1 AV18, u=4 CD-Bonn, u=6 pD AV18,
u=9 pD CD-Bonn.

Pre-FSI observables come from the generated triplet (p1_gen = struck nucleon,
p2_gen/p3_gen = the other two); post-FSI ones from the final-state particle list.
Everything is weighted by the event weight; yields are normalised per generation
attempt (hAttempts) so cross-section-level differences survive.
"""
import os, sys, glob
import numpy as np, awkward as ak, uproot
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

EV = "/home/ilayws/Ilay-Generators/src/programs/genQE_3N/events/wf_test"
OUT = os.path.dirname(os.path.abspath(__file__))
GeVfm = 0.1973269804
M_P = 0.93827
LABELS = {1: "AV18 (u=1)", 4: "CD-Bonn (u=4)", 6: "pD AV18 (u=6)", 9: "pD CD-Bonn (u=9)",
          "av8": "AV8 old code"}
COL = {1: "#3366bb", 4: "#118844", 6: "#cc3355", 9: "#ee9944", "av8": "#888888"}
LS = {1: "-", 4: "-", 6: "--", 9: "--", "av8": ":"}

OBS = [  # key, label, bins
    ("p1", r"$|p_1^{\rm init}| = |p_1-q|$ (struck, initial) [GeV/c]", np.linspace(0, 1.2, 61)),
    ("p_lead", r"$|p_{\rm lead}|$ (outgoing) [GeV/c]", np.linspace(0, 3.0, 61)),
    ("pmiss", r"$p_{\rm miss}=|p_1-q|$ [GeV/c]", np.linspace(0, 1.2, 61)),
    ("p2", r"$|p_2|$ [GeV/c]", np.linspace(0, 1.2, 61)),
    ("p3", r"$|p_3|$ [GeV/c]", np.linspace(0, 1.2, 61)),
    ("pb_tab", r"$p_b$ = half pp-relative (table axis) [GeV/c]", np.linspace(0, 1.0, 61)),
    ("pa_tab", r"$p_a$ = neutron momentum (table axis) [GeV/c]", np.linspace(0, 1.2, 61)),
    ("theta_ab_tab", r"$\theta_{ab}$ in table frame [deg]", np.linspace(0, 180, 61)),
    ("ptriplet", r"$|p_1+p_2+p_3|$ [GeV/c]", np.linspace(0, 1.2, 61)),
    ("Q2", r"$Q^2$ [GeV$^2$]", np.linspace(0, 6, 61)),
    ("xB", r"$x_B$", np.linspace(0.8, 3.0, 61)),
]


def analyse(path, fsi):
    f = uproot.open(path)
    att = float(f["hAttempts"].values()[0]) if "hAttempts" in [k.split(";")[0] for k in f.keys()] else np.nan
    t = f["genT"]
    have = set(t.keys())
    # The pre-replacement (AV8) code writes pLead_pre/p2_pre/p3_pre and has no
    # flat final-state particle list; the current code writes p*_gen + part_*.
    old_schema = "p1_gen" not in have
    gen = ("pLead_pre", "p2_pre", "p3_pre") if old_schema else ("p1_gen", "p2_gen", "p3_gen")
    br = ["weight", "pe", *gen, "N1_type", "N2_type", "N3_type"]
    has_list = fsi and not old_schema
    if has_list:
        br += ["part_pdg", "part_px", "part_py", "part_pz"]
    a = t.arrays(br, library="ak")
    w = ak.to_numpy(a.weight)
    pe = ak.to_numpy(a.pe)
    p1 = ak.to_numpy(a[gen[0]]); p2 = ak.to_numpy(a[gen[1]]); p3 = ak.to_numpy(a[gen[2]])
    Ebeam = 6.0
    q = np.column_stack([-pe[:, 0], -pe[:, 1], Ebeam - pe[:, 2]])
    nu = Ebeam - np.linalg.norm(pe, axis=1)
    Q2 = np.sum(q**2, axis=1) - nu**2
    # p1_gen / pLead_pre is the OUTGOING lead (it has already absorbed q);
    # the INITIAL struck momentum is p1 - q, i.e. the p_miss vector. All the
    # wavefunction (Jacobi) variables must use the initial momenta.
    p1_init = p1 - q
    d = dict(w=w, attempts=att)
    d["p_lead"] = np.linalg.norm(p1, axis=1)
    d["p1"] = np.linalg.norm(p1_init, axis=1)
    d["p2"] = np.linalg.norm(p2, axis=1)
    d["p3"] = np.linalg.norm(p3, axis=1)
    d["pmiss"] = np.linalg.norm(p1 - q, axis=1)
    # Jacobi variables, the SAME ones the density table is binned in
    # (Andrew's convention): p_a = p1 - P_cm/3 ("k_cm" axis of the table),
    # p_b = (p2 - p3)/2 ("k_rel" axis), theta_ab = angle(p_a, p_b).
    Pcm = p1_init + p2 + p3
    pa = p1_init - Pcm/3.0
    pb = (p2 - p3)/2.0
    # TABLE-frame Jacobi pair: get_rho() permutes so the NEUTRON sits in the
    # "a" slot and the two protons form the pair, then looks up
    # rho(|v_a'|, |v_b'|, angle). Reproduce that permutation here.
    t1 = ak.to_numpy(a.N1_type); t2 = ak.to_numpy(a.N2_type); t3 = ak.to_numpy(a.N3_type)
    NEUT = 2112
    q1, q2, q3 = p1_init.copy(), p2.copy(), p3.copy()
    m2 = (t2 == NEUT)          # neutron is N2 -> swap 1<->2
    q1[m2], q2[m2] = p2[m2], p1_init[m2]
    m3 = (t3 == NEUT)          # neutron is N3 -> swap 1<->3
    q1[m3], q3[m3] = p3[m3], p1_init[m3]
    pa_t = q1 - Pcm/3.0
    pb_t = (q2 - q3)/2.0
    d["pa_tab"] = np.linalg.norm(pa_t, axis=1)
    d["pb_tab"] = np.linalg.norm(pb_t, axis=1)
    cabt = np.sum(pa_t*pb_t, axis=1)/np.maximum(d["pa_tab"]*d["pb_tab"], 1e-12)
    d["theta_ab_tab"] = np.degrees(np.arccos(np.clip(cabt, -1, 1)))
    d["n_is_struck"] = (t1 == NEUT).astype(float)
    d["pa"] = np.linalg.norm(pa, axis=1)
    d["pb"] = np.linalg.norm(pb, axis=1)
    cab = np.sum(pa*pb, axis=1)/np.maximum(d["pa"]*d["pb"], 1e-12)
    d["theta_ab"] = np.degrees(np.arccos(np.clip(cab, -1, 1)))
    d["ptriplet"] = np.linalg.norm(Pcm, axis=1)
    d["Q2"] = Q2
    d["xB"] = Q2 / (2 * M_P * np.maximum(nu, 1e-9))
    if has_list:
        pm = np.sqrt(a.part_px**2 + a.part_py**2 + a.part_pz**2)
        isp = a.part_pdg == 2212
        d["n_p_hard"] = ak.to_numpy(ak.sum(isp & (pm > 0.35), axis=1))
        d["n_N_hard"] = ak.to_numpy(ak.sum(((a.part_pdg == 2212) | (a.part_pdg == 2112)) & (pm > 0.35), axis=1))
        d["n_part"] = ak.to_numpy(ak.num(a.part_pdg))
    elif fsi:
        n = len(w)
        d["n_p_hard"] = np.full(n, np.nan); d["n_N_hard"] = np.full(n, np.nan)
        d["n_part"] = np.full(n, np.nan)
    return d


samples = {}
for u in (1, 4, 6, 9):
    for mode, fsi in (("pwia", False), ("fsi", True)):
        p = f"{EV}/{mode}_u{u}.root"
        if os.path.exists(p):
            samples[(mode, u)] = analyse(p, fsi)
            print(f"loaded {mode} u={u}: {len(samples[(mode,u)]['w'])} events", flush=True)
for mode, fsi in (("pwia", False), ("fsi", True)):
    p = f"{EV}/{mode}_av8.root"
    if os.path.exists(p):
        samples[(mode, "av8")] = analyse(p, fsi)
        print(f"loaded {mode} av8: {len(samples[(mode,'av8')]['w'])} events", flush=True)

lines = []
lines.append("=== yields: sum(weight) per generation attempt (cross-section scale) ===")
lines.append(f"{'sample':22s} {'events':>10s} {'attempts':>12s} {'sum w':>12s} {'sum w / attempt':>16s} {'ratio to AV18':>14s}")
for mode in ("pwia", "fsi"):
    ref = None
    for u in (1, 4, 6, 9, "av8"):
        if (mode, u) not in samples:
            continue
        d = samples[(mode, u)]
        y = d["w"].sum() / d["attempts"]
        if u == 1:
            ref = y
        r = y / ref if ref else float("nan")
        lines.append(f"{mode+' '+LABELS[u]:22s} {len(d['w']):10d} {d['attempts']:12.0f} "
                     f"{d['w'].sum():12.5g} {y:16.6g} {r:14.4f}")
    lines.append("")

# ---------- distribution plots ----------
for mode in ("pwia", "fsi"):
    keys = [k for k in samples if k[0] == mode]
    if not keys:
        continue
    fig, axes = plt.subplots(2, 5, figsize=(22, 8))
    for ax, (key, xlabel, bins) in zip(axes.ravel(), OBS):
        for (m, u) in keys:
            d = samples[(m, u)]
            h, _ = np.histogram(d[key], bins=bins, weights=d["w"])
            if h.sum() > 0:
                ax.step(0.5*(bins[1:]+bins[:-1]), h/h.sum(), where="mid",
                        color=COL[u], ls=LS[u], lw=1.8, label=LABELS[u])
        ax.set_xlabel(xlabel); ax.set_ylim(bottom=0)
    axes[0, 0].set_ylabel("normalised"); axes[1, 0].set_ylabel("normalised")
    axes[0, 0].legend(frameon=False, fontsize=9)
    fig.suptitle(f"3N generated sample, {mode.upper()} — wavefunction comparison (weighted, area-normalised)", fontsize=13)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, f"{mode}_observables.png"), dpi=160); plt.close(fig)

    # log-scale tails for the momentum observables
    fig, axes = plt.subplots(1, 4, figsize=(18, 4.3))
    for ax, key in zip(axes, ("p1", "pmiss", "pb_tab", "pa_tab")):
        bins = dict(OBS and [(k, b) for k, _, b in OBS])[key]
        for (m, u) in keys:
            d = samples[(m, u)]
            h, _ = np.histogram(d[key], bins=bins, weights=d["w"])
            if h.sum() > 0:
                ax.semilogy(0.5*(bins[1:]+bins[:-1]), h/h.sum(), color=COL[u], ls=LS[u], lw=1.8,
                            label=LABELS[u])
        ax.set_xlabel([o[1] for o in OBS if o[0] == key][0])
    axes[0].set_ylabel("normalised"); axes[0].legend(frameon=False, fontsize=9)
    fig.suptitle(f"{mode.upper()} — high-momentum tails (log)", fontsize=13)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, f"{mode}_tails.png"), dpi=160); plt.close(fig)

# ---------- ratio plots: pD / full, same potential ----------
for mode in ("pwia", "fsi"):
    pairs = [(6, 1, "pD AV18 / AV18"), (9, 4, "pD CD-Bonn / CD-Bonn"), (4, 1, "CD-Bonn / AV18")]
    pairs = [(a, b, l) for a, b, l in pairs if (mode, a) in samples and (mode, b) in samples]
    if not pairs:
        continue
    fig, axes = plt.subplots(1, 4, figsize=(18, 4.3))
    for ax, key in zip(axes, ("p1", "pmiss", "pb_tab", "theta_ab_tab")):
        bins = [b for k, _, b in OBS if k == key][0]
        c = 0.5*(bins[1:]+bins[:-1])
        for num, den, lab in pairs:
            dn, dd = samples[(mode, num)], samples[(mode, den)]
            hn, _ = np.histogram(dn[key], bins=bins, weights=dn["w"])
            hd, _ = np.histogram(dd[key], bins=bins, weights=dd["w"])
            hn, hd = hn/hn.sum(), hd/hd.sum()
            r = np.divide(hn, hd, out=np.full_like(hn, np.nan), where=hd > 0)
            ax.plot(c, r, lw=1.8, label=lab)
        ax.axhline(1.0, color="k", lw=0.6, ls=":")
        ax.set_xlabel([o[1] for o in OBS if o[0] == key][0]); ax.set_ylim(0, 2.5)
    axes[0].set_ylabel("shape ratio"); axes[0].legend(frameon=False, fontsize=9)
    fig.suptitle(f"{mode.upper()} — shape ratios (area-normalised distributions)", fontsize=13)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, f"{mode}_shape_ratios.png"), dpi=160); plt.close(fig)

# ---------- numbers ----------
lines.append("=== weighted means / tail fractions ===")
hdr = f"{'sample':22s}" + "".join(f"{k:>14s}" for k in ("<p1>", "<pmiss>", "<p_b>tab", "<p_a>tab", "<th_ab>tab", "f(p1>0.4)"))
for mode in ("pwia", "fsi"):
    lines.append(f"-- {mode.upper()}")
    lines.append(hdr)
    for u in (1, 4, 6, 9, "av8"):
        if (mode, u) not in samples: continue
        d = samples[(mode, u)]; w = d["w"]; W = w.sum()
        vals = [np.sum(w*d[k])/W for k in ("p1", "pmiss", "pb_tab", "pa_tab", "theta_ab_tab")]
        vals.append(w[d["p1"] > 0.4].sum()/W)
        lines.append(f"{LABELS[u]:22s}" + "".join(f"{v:14.4f}" for v in vals))
    lines.append("")

if ("fsi", 1) in samples:
    lines.append("=== post-FSI final-state multiplicity (weighted) ===")
    lines.append(f"{'sample':22s} {'<n protons>0.35>':>18s} {'<n nucleons>0.35>':>19s} {'<n particles>':>14s}")
    for u in (1, 4, 6, 9, "av8"):
        if ("fsi", u) not in samples: continue
        d = samples[("fsi", u)]; w = d["w"]; W = w.sum()
        lines.append(f"{LABELS[u]:22s} {np.sum(w*d['n_p_hard'])/W:18.4f} "
                     f"{np.sum(w*d['n_N_hard'])/W:19.4f} {np.sum(w*d['n_part'])/W:14.4f}")

open(os.path.join(OUT, "samples_summary.txt"), "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
print("\nsaved plots + samples_summary.txt in", OUT)
