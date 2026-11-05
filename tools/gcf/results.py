"""The numbers of an event file: size, normalization, statistical quality, provenance.

Reads only the branches it needs, a chunk at a time, so a file of any size fits
in memory. Uses uproot and numpy.
"""

import json
import math
import re
from pathlib import Path

# Below these N_eff a sample gets a warning: its relative statistical error on
# the yield is 1/sqrt(N_eff), and its histograms are spiky.
NEFF_BAD = 100
NEFF_WARN = 1000
# A warning if one event carries more than this fraction of the summed weight.
MAX_WEIGHT_WARN = 0.05

# What Σw/attempts of each model can be compared with (docs/PHYSICS.md).
SCALE_NOTES = {
    "3N": "The absolute scale of the 3N channel is a placeholder (#73): compare Σw/attempts "
          "between 3N samples only, not with 2N or MF.",
    "MF": "The mean-field weight is not on the same absolute scale as the SRC weight (#34): "
          "compare Σw/attempts between MF samples only.",
}


def _need_modules():
    try:
        import numpy  # noqa: F401
        import uproot  # noqa: F401
    except ImportError as e:
        raise RuntimeError(f"the run summary needs uproot and numpy ({e})") from None


def _provenance(f):
    """Per stage: its code versions, seeds and first command, from the TNamed objects.
    A merged file holds one per job, as cycles of the same name."""
    out = {}
    for key in f.keys(cycle=True):
        name = key.split(";")[0]
        if f.classname_of(key).split("::")[-1] != "TNamed":
            continue
        fields = {}
        for line in str(f[key].member("fTitle")).splitlines():
            k, _, v = line.partition(":")
            if v:
                fields[k.strip()] = v.strip()
        if "version" not in fields:
            continue
        s = out.setdefault(name, {"versions": [], "seeds": [], "command": fields.get("command", ""),
                                  "date": fields.get("date", ""), "model": fields.get("model"),
                                  "inputs": fields.get("nuclear inputs")})
        if fields["version"] not in s["versions"]:
            s["versions"].append(fields["version"])
        if "seed" in fields:
            s["seeds"].append(fields["seed"])
    return out


def summarize(path, tree=None, attempts=None):
    """A dict of plain numbers for an event file; see the keys below. `tree` is the
    event-file format's or, without one, that of the suite's other generators;
    `attempts` is for files that don't record theirs."""
    _need_modules()
    import numpy as np
    import uproot

    path = Path(path)
    with uproot.open(path) as f:
        keys = {k.split(";")[0] for k in f.keys()}
        if tree is None:
            tree = "events" if "events" in keys or "genT" not in keys else "genT"
        have = set(f[tree].keys())
        attempts = float(f["hAttempts"].values()[0]) if "hAttempts" in keys else attempts
        provenance = _provenance(f)
        n = sum_w = sum_w2 = max_w = 0.
        fsi_ran = 0
        branches = [b for b in ("weight", "fsi_ran") if b in have]
        if not branches:
            n = f[tree].num_entries
        else:
            for a in f[tree].iterate(branches, step_size="100 MB", library="np"):
                w = a["weight"] if "weight" in a else np.ones(len(a["fsi_ran"]))
                n += len(w)
                sum_w += float(w.sum())
                sum_w2 += float((w * w).sum())
                if len(w):
                    max_w = max(max_w, float(w.max()))
                if "fsi_ran" in a:
                    fsi_ran += int(a["fsi_ran"].sum())
    n = int(n)
    has_w = "weight" in have
    out = {"events": n, "attempts": attempts, "has_weights": has_w,
           "sum_w": sum_w if has_w else None,
           "yield": sum_w / attempts if has_w and attempts else None,
           "yield_err": math.sqrt(sum_w2) / attempts if has_w and attempts else None,
           "n_eff": sum_w * sum_w / sum_w2 if has_w and sum_w2 > 0 else None,
           "efficiency": n / attempts if attempts else None,
           "max_w_fraction": max_w / sum_w if has_w and sum_w > 0 else None,
           "fsi_ran_fraction": fsi_ran / n if "fsi_ran" in have and n else None,
           "provenance": provenance,
           "model": next((s["model"] for s in provenance.values() if s.get("model")), None)}
    out["notes"] = notes(out)
    return out


def notes(s):
    """[(level, text)] about a summary: 'bad', 'warn' or 'info'."""
    out = []
    if s.get("events") == 0:
        out.append(("bad", "No attempt had a positive weight, so this file holds no events. "
                           "Check the beam energy against the reaction's threshold, and the target and limits."))
    n_eff = s.get("n_eff")
    if n_eff is not None:
        err = 100 / math.sqrt(n_eff) if n_eff > 0 else float("inf")
        if n_eff < NEFF_BAD:
            out.append(("bad", f"N_eff = {n_eff:.0f} of {s['events']} events: a few events dominate this "
                               f"sample, so its shapes are unreliable and its yield has a {err:.0f}% "
                               f"statistical error. Generate more events."))
        elif n_eff < NEFF_WARN:
            out.append(("warn", f"N_eff = {n_eff:.0f} of {s['events']} events: the yield has a {err:.0f}% "
                                f"statistical error."))
    mw = s.get("max_w_fraction")
    if mw is not None and mw > MAX_WEIGHT_WARN:
        out.append(("warn", f"The heaviest event carries {100 * mw:.0f}% of the summed weight."))
    fr = s.get("fsi_ran_fraction")
    if fr is not None and fr < 1:
        out.append(("info", f"The cascade ran for {100 * fr:.1f}% of the events. The others kept their "
                            f"input particles (residual A < 3, or nothing to transport)."))
    if s.get("attempts") is None:
        out.append(("info", "No hAttempts in this file, so it has no normalization per attempt."))
    versions = {v for st in s.get("provenance", {}).values() for v in st["versions"]}
    if any(v.endswith("-dirty") or v == "unknown" for v in versions):
        out.append(("warn", "Made by a build with uncommitted changes or an unknown version: "
                            + ", ".join(sorted(versions)) + "."))
    elif len(versions) > 1:
        out.append(("warn", "Its stages were built from different code versions: "
                            + ", ".join(sorted(versions)) + "."))
    if s.get("model") in SCALE_NOTES:
        out.append(("info", SCALE_NOTES[s["model"]]))
    return out


def cached(path, cache_file, key):
    """The summary stored under `key` in `cache_file`, if `path` has not changed since."""
    try:
        entry = json.loads(Path(cache_file).read_text()).get(key)
        st = Path(path).stat()
    except (OSError, ValueError):
        return None
    if entry and entry.get("_mtime") == st.st_mtime and entry.get("_size") == st.st_size:
        return entry
    return None


def store(path, cache_file, key, summary):
    """Add a summary to `cache_file`; returns it with the file's identity."""
    cache_file = Path(cache_file)
    try:
        data = json.loads(cache_file.read_text())
    except (OSError, ValueError):
        data = {}
    st = Path(path).stat()
    summary = dict(summary, _mtime=st.st_mtime, _size=st.st_size)
    data[key] = summary
    tmp = cache_file.with_suffix(".tmp")
    tmp.write_text(json.dumps(data, indent=1))
    tmp.replace(cache_file)
    return summary


def text(s):
    """A summary as lines of text, for the terminal."""
    lines = [f"events        {s['events']}"]
    if s["attempts"]:
        lines.append(f"attempts      {s['attempts']:.0f}  (efficiency {100 * s['efficiency']:.3g}%)")
    if s["yield"] is not None:
        lines.append(f"sum(w)/attempts  {s['yield']:.5g} +- {s['yield_err']:.2g}")
    if s["n_eff"] is not None:
        lines.append(f"N_eff         {s['n_eff']:.0f}")
    if s["max_w_fraction"] is not None:
        lines.append(f"heaviest event  {100 * s['max_w_fraction']:.2g}% of sum(w)")
    if s["fsi_ran_fraction"] is not None:
        lines.append(f"FSI ran for   {100 * s['fsi_ran_fraction']:.1f}%")
    for stage, p in s["provenance"].items():
        lines.append(f"{stage}: version {', '.join(p['versions'])}, {len(p['seeds'])} job(s), {p['date']}")
    lines += [f"[{level}] {re.sub(r'Σw', 'sum(w)', t)}" for level, t in s["notes"]]
    return "\n".join(lines)
