"""A page of weighted control plots for an event file (PNG).

Uses uproot, awkward and matplotlib. Reads the file a chunk at a time and fills
histograms, so its memory use does not grow with the file. Panels that need
branches the file doesn't have (e.g. no electron in a file from another
generator, or no FSI yet) are left out.

A selection can be applied: `cuts` maps q2, xb, theta_e, pmiss and protons to
[min, max], either of which may be None.
"""

from pathlib import Path

M_N = 0.93892  # average nucleon mass (GeV)
PROTON, NEUTRON = 2212, 2112
COLUMNS = 2
CUT_LABELS = {"q2": "Q²", "xb": "x_B", "theta_e": "θ_e", "pmiss": "p_miss", "protons": "protons"}
CUT_UNITS = {"q2": " GeV²", "xb": "", "theta_e": "°", "pmiss": " GeV/c", "protons": ""}


def _need_modules():
    try:
        import awkward  # noqa: F401
        import matplotlib  # noqa: F401
        import numpy  # noqa: F401
        import uproot  # noqa: F401
    except ImportError as e:
        raise RuntimeError(f"quicklook needs uproot, awkward, numpy and matplotlib ({e})") from None


def check_cuts(cuts):
    """Cuts with only the set limits, as floats; ValueError if they make no sense."""
    out = {}
    for key, limits in (cuts or {}).items():
        if key not in CUT_LABELS:
            raise ValueError(f"unknown cut {key!r}; allowed: {', '.join(CUT_LABELS)}")
        if not isinstance(limits, (list, tuple)) or len(limits) != 2:
            raise ValueError(f"cut {key}: expected [min, max], got {limits!r}")
        lo, hi = (None if x is None else float(x) for x in limits)
        if lo is None and hi is None:
            continue
        if lo is not None and hi is not None and hi < lo:
            raise ValueError(f"cut {key}: max must not be below min")
        out[key] = [lo, hi]
    return out


def describe_cuts(cuts):
    parts = []
    for key, (lo, hi) in cuts.items():
        name, unit = CUT_LABELS[key], CUT_UNITS[key]
        if lo is not None and hi is not None:
            parts.append(f"{lo:g} < {name} < {hi:g}{unit}")
        elif lo is not None:
            parts.append(f"{name} > {lo:g}{unit}")
        else:
            parts.append(f"{name} < {hi:g}{unit}")
    return ", ".join(parts)


class Panel:
    """A 1D plot of one or more weighted series with common bin edges."""

    def __init__(self, xlabel, edges, integer=False):
        self.xlabel, self.edges, self.series, self.integer = xlabel, edges, {}, integer

    def fill(self, label, values, weights):
        import numpy as np
        values = np.asarray(values, dtype=float)
        ok = np.isfinite(values)
        counts, _ = np.histogram(values[ok], bins=self.edges, weights=np.asarray(weights)[ok])
        if label in self.series:
            self.series[label] += counts
        else:
            self.series[label] = counts

    def empty(self):
        return not any(c.sum() > 0 for c in self.series.values())

    def draw(self, ax):
        import numpy as np
        edges = self.edges
        if self.integer:         # up to the last filled bin
            last = max((int(np.nonzero(c)[0][-1]) for c in self.series.values() if c.sum() > 0), default=0)
            edges = edges[:last + 2]
        for label, counts in self.series.items():
            total = counts.sum()
            if total <= 0:
                continue
            c = counts[:len(edges) - 1]
            ax.stairs(c / total / np.diff(edges), edges, linewidth=1.5, label=label or None)
        ax.set_xlabel(self.xlabel)
        ax.set_yticks([])
        if any(self.series):
            ax.legend(fontsize=9)


class Panel2D:
    """A 2D map of a weighted pair of values."""

    def __init__(self, xlabel, ylabel, edges):
        import numpy as np
        self.xlabel, self.ylabel, self.edges = xlabel, ylabel, edges
        self.counts = np.zeros((len(edges) - 1, len(edges) - 1))

    def fill(self, x, y, weights):
        import numpy as np
        x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
        ok = np.isfinite(x) & np.isfinite(y)
        c, _, _ = np.histogram2d(x[ok], y[ok], bins=[self.edges, self.edges], weights=np.asarray(weights)[ok])
        self.counts += c

    def empty(self):
        return self.counts.sum() <= 0

    def draw(self, ax):
        ax.pcolormesh(self.edges, self.edges, self.counts.T, cmap="viridis", shading="flat")
        ax.set_xlabel(self.xlabel)
        ax.set_ylabel(self.ylabel)


def _angle(ax_, ay, az, bx, by, bz):
    """Angle between two vectors, in degrees; NaN where one is missing or null."""
    import numpy as np
    with np.errstate(invalid="ignore", divide="ignore"):
        cos = (ax_ * bx + ay * by + az * bz) / np.sqrt((ax_**2 + ay**2 + az**2) * (bx**2 + by**2 + bz**2))
    return np.degrees(np.arccos(np.clip(cos, -1, 1)))


def _role(a, have, role):
    """px, py, pz of the nucleon with a role (NaN if the event has none)."""
    import awkward as ak
    import numpy as np
    roles = a["nuc_role"] if "nuc_role" in have else ak.local_index(a["nuc_pdg"])
    sel = roles == role
    return tuple(np.asarray(ak.fill_none(ak.firsts(a[k][sel]), np.nan), dtype=float)
                 for k in ("nuc_px", "nuc_py", "nuc_pz"))


def _fastest_proton(pdg, px, py, pz):
    import awkward as ak
    import numpy as np
    p = np.sqrt(px**2 + py**2 + pz**2)
    return np.asarray(ak.fill_none(ak.max(p[pdg == PROTON], axis=1), np.nan), dtype=float)


def make(path, output=None, tree="events", cuts=None, stats=None, step_size="200 MB"):
    """Draw the page; returns its path. `stats`, if given, is filled with the
    numbers of the selection and the filled histograms."""
    _need_modules()
    import awkward as ak
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import uproot

    path = Path(path)
    cuts = check_cuts(cuts)
    if output:
        output = Path(output)
    else:
        output = path.with_name(path.stem + ("_cuts.png" if cuts else ".png"))

    f = uproot.open(path)
    t = f[tree]
    have = set(t.keys())
    branches = [b for b in ["weight", "beam_E", "e_px", "e_py", "e_pz", "e_E",
                            "nuc_pdg", "nuc_role", "nuc_px", "nuc_py", "nuc_pz",
                            "part_pdg", "part_px", "part_py", "part_pz", "fsi_ran"] if b in have]
    has_e = {"beam_E", "e_px", "e_py", "e_pz", "e_E"} <= have
    has_nuc = {"nuc_pdg", "nuc_px", "nuc_py", "nuc_pz"} <= have
    has_fsi = {"part_pdg", "part_px", "part_py", "part_pz"} <= have
    if not has_e and not has_nuc:
        raise RuntimeError(f"{path}: no branches to plot (see docs/PIPELINE.md for the format)")
    for key in cuts:
        if key == "protons" and not (has_nuc or has_fsi):
            raise ValueError("cut protons: this file has no particle list")
        if key == "pmiss" and not (has_e and has_nuc):
            raise ValueError("cut pmiss: this file has no electron or struck nucleon")
        if key in ("q2", "xb", "theta_e") and not has_e:
            raise ValueError(f"cut {key}: this file has no electron")

    attempts = None
    if "hAttempts" in [k.split(";")[0] for k in f.keys()]:
        attempts = float(f["hAttempts"].values()[0])

    panels = None
    n_all = n_sel = 0
    sum_all = sum_sel = sum2_sel = fsi_ran = 0.
    integer_edges = np.arange(-0.5, 40.5)
    for a in t.iterate(branches, step_size=step_size):
        n = len(a)
        if not n:
            continue
        w = np.asarray(a["weight"], dtype=float) if "weight" in have else np.ones(n)
        v = {}
        if has_e:
            beam = np.asarray(a["beam_E"], dtype=float)
            ex, ey, ez, eE = (np.asarray(a[k], dtype=float) for k in ("e_px", "e_py", "e_pz", "e_E"))
            qx, qy, qz, nu = -ex, -ey, beam - ez, beam - eE
            v["q2"] = qx**2 + qy**2 + qz**2 - nu**2
            with np.errstate(invalid="ignore", divide="ignore"):
                v["xb"] = v["q2"] / (2 * M_N * nu)
            v["theta_e"] = np.degrees(np.arctan2(np.hypot(ex, ey), ez))
        if has_nuc:
            v["protons_pre"] = np.asarray(ak.sum(a["nuc_pdg"] == PROTON, axis=1))
            v["lead_pre"] = _fastest_proton(a["nuc_pdg"], a["nuc_px"], a["nuc_py"], a["nuc_pz"])
            p2, p3 = _role(a, have, 1), _role(a, have, 2)
            v["p2"] = np.sqrt(p2[0]**2 + p2[1]**2 + p2[2]**2)
            v["p3"] = np.sqrt(p3[0]**2 + p3[1]**2 + p3[2]**2)
            v["theta_23"] = _angle(*p2, *p3)
            if has_e:
                s = _role(a, have, 0)
                # The struck nucleon's momentum before the electron hit it.
                p1 = (s[0] - qx, s[1] - qy, s[2] - qz)
                v["pmiss"] = np.sqrt(p1[0]**2 + p1[1]**2 + p1[2]**2)
                v["theta_12"] = _angle(*p1, *p2)
                v["theta_13"] = _angle(*p1, *p3)
        if has_fsi:
            pdg = a["part_pdg"]
            v["protons_post"] = np.asarray(ak.sum(pdg == PROTON, axis=1))
            v["neutrons_post"] = np.asarray(ak.sum(pdg == NEUTRON, axis=1))
            v["pions_post"] = np.asarray(ak.sum(abs(pdg) == 211, axis=1))
            v["lead_post"] = _fastest_proton(pdg, a["part_px"], a["part_py"], a["part_pz"])

        keep = np.ones(n, dtype=bool)
        for key, (lo, hi) in cuts.items():
            x = v.get("protons_post", v.get("protons_pre")) if key == "protons" else v[key]
            x = np.asarray(x, dtype=float)
            if lo is not None:
                keep &= x >= lo
            if hi is not None:
                keep &= x <= hi
        n_all += n
        n_sel += int(keep.sum())
        sum_all += float(w.sum())
        sum_sel += float(w[keep].sum())
        sum2_sel += float((w[keep]**2).sum())
        if "fsi_ran" in have:
            fsi_ran += float(np.asarray(a["fsi_ran"], dtype=float)[keep].sum())

        if panels is None:       # the data-driven ranges, from the first chunk
            def top(x, q, floor):
                x = np.asarray(x, dtype=float)
                x = x[np.isfinite(x)]
                return max(floor, float(np.percentile(x, q)) * 1.1) if len(x) else floor
            panels = {}
            if has_e:
                panels["q2"] = Panel("Q$^2$ [GeV$^2$]", np.linspace(0, top(v["q2"], 99.9, 1.0), 61))
                panels["xb"] = Panel("x$_B$", np.linspace(0, 3, 61))
                panels["theta_e"] = Panel("electron $\\theta$ [deg]", np.linspace(0, min(180.0, top(v["theta_e"], 100, 5.0)), 61))
            if has_nuc:
                if has_e:
                    panels["pmiss"] = Panel("p$_{miss}$ = p$_{struck}$ $-$ q, before FSI [GeV/c]", np.linspace(0, 1.5, 61))
                panels["partners"] = Panel("partner |p|, before FSI [GeV/c]", np.linspace(0, 1.5, 61))
                if has_e:
                    panels["opening"] = Panel("angle to the struck nucleon's p$_{miss}$, before FSI [deg]", np.linspace(0, 180, 46))
                    panels["theta_map"] = Panel2D("$\\theta_{12}$ (struck, partner 1) [deg]", "$\\theta_{23}$ (partners) [deg]",
                                                  np.linspace(0, 180, 37))
                panels["protons"] = Panel("protons per event", integer_edges, integer=True)
                panels["lead"] = Panel("fastest proton |p| [GeV/c]", np.linspace(0, top(v["lead_pre"], 99.5, 1.0), 61))
            if has_fsi:
                panels["neutrons"] = Panel("neutrons per event, after FSI", integer_edges, integer=True)
                panels["pions"] = Panel("charged pions per event, after FSI", integer_edges, integer=True)

        ws = w[keep]
        if has_e:
            panels["q2"].fill("", v["q2"][keep], ws)
            panels["xb"].fill("", v["xb"][keep], ws)
            panels["theta_e"].fill("", v["theta_e"][keep], ws)
        if has_nuc:
            if has_e:
                panels["pmiss"].fill("", v["pmiss"][keep], ws)
                panels["opening"].fill("partner 1", v["theta_12"][keep], ws)
                panels["opening"].fill("partner 2", v["theta_13"][keep], ws)
                panels["theta_map"].fill(v["theta_12"][keep], v["theta_23"][keep], ws)
            panels["partners"].fill("partner 1", v["p2"][keep], ws)
            panels["partners"].fill("partner 2", v["p3"][keep], ws)
            panels["protons"].fill("before FSI", v["protons_pre"][keep], ws)
            panels["lead"].fill("before FSI", v["lead_pre"][keep], ws)
        if has_fsi:
            if has_nuc:
                panels["protons"].fill("after FSI", v["protons_post"][keep], ws)
                panels["lead"].fill("after FSI", v["lead_post"][keep], ws)
            panels["neutrons"].fill("", v["neutrons_post"][keep], ws)
            panels["pions"].fill("", v["pions_post"][keep], ws)

    if panels is None:
        raise RuntimeError(f"{path}: the file has no events")
    if n_sel == 0:
        raise ValueError("no event passes this selection")
    drawn = [p for p in panels.values() if not p.empty()]

    n_eff = sum_sel**2 / sum2_sel if sum2_sel > 0 else 0.
    rows = (len(drawn) + COLUMNS - 1) // COLUMNS
    with plt.rc_context({"font.size": 12}):
        fig, axes = plt.subplots(rows, COLUMNS, figsize=(4.6 * COLUMNS, 3.0 * rows + 0.8), squeeze=False)
        for ax, panel in zip(axes.flat, drawn):
            panel.draw(ax)
        for ax in list(axes.flat)[len(drawn):]:
            ax.axis("off")
        info = f"{path.name}: {n_sel} events (N$_{{eff}}$ = {n_eff:.0f})"
        if attempts:
            info += f", {attempts:.0f} attempts,\n$\\Sigma w$/attempts = {sum_sel / attempts:.4g}"
        if "fsi_ran" in have and n_sel:
            info += f", FSI ran for {100 * fsi_ran / n_sel:.1f}%"
        if cuts:
            info += f"\n{describe_cuts(cuts)}: {100 * sum_sel / sum_all:.3g}% of $\\Sigma w$"
        fig.suptitle(info + "  (weighted, unit area)", fontsize=11)
        fig.tight_layout()
        fig.savefig(output, dpi=110)
        plt.close(fig)
    if stats is not None:
        stats.update(events=n_sel, events_all=n_all, n_eff=n_eff, sum_w=sum_sel, sum_w_all=sum_all,
                     fraction=sum_sel / sum_all if sum_all > 0 else None,
                     yield_=sum_sel / attempts if attempts else None,
                     yield_err=sum2_sel**0.5 / attempts if attempts else None,
                     cuts=cuts, description=describe_cuts(cuts), panels=[p.xlabel for p in drawn],
                     histograms={p.xlabel: {k: c.tolist() for k, c in p.series.items()}
                                 for p in drawn if isinstance(p, Panel)})
    return output
