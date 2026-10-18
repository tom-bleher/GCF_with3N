"""A page of weighted control plots for an event file (PNG).

Uses uproot, awkward and matplotlib. Panels that need branches the file
doesn't have (e.g. no electron in a file from another generator, or no FSI
yet) are left out.
"""

from pathlib import Path

M_N = 0.93892  # average nucleon mass (GeV)
PROTON = 2212


def _need_modules():
    try:
        import awkward  # noqa: F401
        import matplotlib  # noqa: F401
        import numpy  # noqa: F401
        import uproot  # noqa: F401
    except ImportError as e:
        raise RuntimeError(f"quicklook needs uproot, awkward, numpy and matplotlib ({e})") from None


def make(path, output=None, tree="events"):
    _need_modules()
    import awkward as ak
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import uproot

    path = Path(path)
    output = Path(output) if output else path.with_suffix(".png")
    f = uproot.open(path)
    t = f[tree]
    have = set(t.keys())
    branches = [b for b in ["weight", "beam_E", "e_px", "e_py", "e_pz", "e_E",
                            "nuc_pdg", "nuc_role", "nuc_px", "nuc_py", "nuc_pz",
                            "part_pdg", "part_px", "part_py", "part_pz", "fsi_ran"] if b in have]
    a = t.arrays(branches)
    n = len(a)
    w = np.asarray(a["weight"]) if "weight" in have else np.ones(n)
    attempts = None
    if "hAttempts" in [k.split(";")[0] for k in f.keys()]:
        attempts = float(f["hAttempts"].values()[0])

    panels = []

    def hist(ax, values, weights, bins, label=None, **kw):
        values, weights = np.asarray(values), np.asarray(weights)
        ok = np.isfinite(values)
        if ok.sum() == 0:
            return
        ax.hist(values[ok], bins=bins, weights=weights[ok], density=True, histtype="step",
                linewidth=1.5, label=label, **kw)

    has_e = {"beam_E", "e_px", "e_py", "e_pz", "e_E"} <= have
    has_nuc = {"nuc_pdg", "nuc_px", "nuc_py", "nuc_pz"} <= have
    has_fsi = {"part_pdg", "part_px", "part_py", "part_pz"} <= have

    if has_e:
        beam = np.asarray(a["beam_E"])
        ex, ey, ez, eE = (np.asarray(a[k]) for k in ("e_px", "e_py", "e_pz", "e_E"))
        qx, qy, qz, nu = -ex, -ey, beam - ez, beam - eE
        Q2 = qx**2 + qy**2 + qz**2 - nu**2
        xB = Q2 / (2 * M_N * nu)
        theta = np.degrees(np.arctan2(np.hypot(ex, ey), ez))
        panels.append(("Q$^2$ [GeV$^2$]", lambda ax: hist(ax, Q2, w, 60)))
        panels.append(("x$_B$", lambda ax: hist(ax, xB, w, np.linspace(0, 3, 61))))
        panels.append(("electron $\\theta$ [deg]", lambda ax: hist(ax, theta, w, 60)))
        if has_nuc:
            role = a["nuc_role"] if "nuc_role" in have else ak.local_index(a["nuc_pdg"])
            struck = role == 0
            px = ak.fill_none(ak.firsts(a["nuc_px"][struck]), np.nan)
            py = ak.fill_none(ak.firsts(a["nuc_py"][struck]), np.nan)
            pz = ak.fill_none(ak.firsts(a["nuc_pz"][struck]), np.nan)
            pmiss = np.sqrt((np.asarray(px) - qx)**2 + (np.asarray(py) - qy)**2
                            + (np.asarray(pz) - qz)**2)
            panels.append(("p$_{miss}$ = p$_{struck}$ $-$ q, before FSI [GeV/c]",
                           lambda ax: hist(ax, pmiss, w, np.linspace(0, 1.5, 61))))

    def lead_proton(pdg, px, py, pz):
        p = np.sqrt(px**2 + py**2 + pz**2)
        p = p[pdg == PROTON]
        return np.asarray(ak.fill_none(ak.max(p, axis=1), np.nan))

    if has_nuc:
        npro_pre = np.asarray(ak.sum(a["nuc_pdg"] == PROTON, axis=1))
        lead_pre = lead_proton(a["nuc_pdg"], a["nuc_px"], a["nuc_py"], a["nuc_pz"])
        npro_post = lead_post = None
        if has_fsi:
            npro_post = np.asarray(ak.sum(a["part_pdg"] == PROTON, axis=1))
            lead_post = lead_proton(a["part_pdg"], a["part_px"], a["part_py"], a["part_pz"])

        def multiplicity(ax):
            top = int(max(npro_pre.max(initial=0),
                          npro_post.max(initial=0) if npro_post is not None else 0)) + 2
            bins = np.arange(-0.5, top + 0.5)
            hist(ax, npro_pre, w, bins, label="before FSI")
            if npro_post is not None:
                hist(ax, npro_post, w, bins, label="after FSI")
            ax.legend(fontsize=8)

        def leading(ax):
            bins = np.linspace(0, max(1.0, np.nanpercentile(lead_pre, 99.5) * 1.1), 61)
            hist(ax, lead_pre, w, bins, label="before FSI")
            if lead_post is not None:
                hist(ax, lead_post, w, bins, label="after FSI")
            ax.legend(fontsize=8)

        panels.append(("protons per event", multiplicity))
        panels.append(("fastest proton |p| [GeV/c]", leading))

    if not panels:
        raise RuntimeError(f"{path}: no branches to plot (see docs/PIPELINE.md for the format)")

    cols = 3
    rows = (len(panels) + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(4.2 * cols, 3.3 * rows), squeeze=False)
    for ax, (title, draw) in zip(axes.flat, panels):
        draw(ax)
        ax.set_xlabel(title)
        ax.set_yticks([])
    for ax in list(axes.flat)[len(panels):]:
        ax.axis("off")
    info = f"{path.name}: {n} events"
    if n and np.sum(w**2) > 0:
        # Effective number of events: few heavy events make spiky histograms.
        info += f" (N$_{{eff}}$ = {w.sum()**2 / np.sum(w**2):.0f})"
    if attempts:
        info += f", {attempts:.0f} attempts, $\\Sigma w$/attempts = {w.sum() / attempts:.4g}"
    if "fsi_ran" in have:
        info += f", FSI ran for {100 * np.mean(np.asarray(a['fsi_ran'])):.1f}%"
    fig.suptitle(info + "  (weighted, unit area)", fontsize=10)
    fig.tight_layout()
    fig.savefig(output, dpi=110)
    plt.close(fig)
    return output
