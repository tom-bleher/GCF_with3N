"""Run configurations: TOML in, validated dict out, and back to TOML.

A config has up to four tables, [run], [generate], [fsi] and [eic]; their keys
are listed in schema.py. Validation fills in defaults, so a saved (resolved)
config records every setting a run used.
"""

import os
import re
import shutil
import tomllib
from pathlib import Path

from . import schema


class ConfigError(Exception):
    pass


def _coerce(opt, value, where):
    kind = opt.kind
    bad = ConfigError(f"{where}: expected {kind}, got {value!r}")
    if kind == "float":
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise bad
        value = float(value)
    elif kind == "int":
        if isinstance(value, bool) or not isinstance(value, int):
            raise bad
    elif kind == "str":
        if not isinstance(value, str):
            raise bad
    elif kind == "bool":
        if not isinstance(value, bool):
            raise bad
    elif kind == "range":
        if (not isinstance(value, list) or len(value) != 2
                or not all(isinstance(v, (int, float)) and not isinstance(v, bool) for v in value)):
            raise ConfigError(f"{where}: expected [min, max], got {value!r}")
        value = [float(value[0]), float(value[1])]
        if value[1] <= value[0]:
            raise ConfigError(f"{where}: max must exceed min, got {value!r}")
    elif kind == "target":
        if (not isinstance(value, list) or len(value) != 2
                or not all(isinstance(v, int) and not isinstance(v, bool) for v in value)):
            raise ConfigError(f"{where}: expected [A, Z], got {value!r}")
        if value[0] < 1 or not 0 <= value[1] <= value[0]:
            raise ConfigError(f"{where}: invalid nucleus A={value[0]}, Z={value[1]}")
    elif kind == "slices":
        value = _coerce_slices(value, where)
    elif kind == "list":
        if not isinstance(value, list) or not all(isinstance(v, str) for v in value):
            raise ConfigError(f"{where}: expected a list of strings, got {value!r}")
    if opt.choices is not None:
        vals = value if kind == "list" else [value]
        for v in vals:
            if v not in opt.choices:
                raise ConfigError(f"{where}: {v!r} is not one of {opt.choices}")
    return value


def _coerce_slices(value, where):
    """[[q2 min, q2 max, attempts], ...]: disjoint Q^2 ranges, in ascending order in the result."""
    number = lambda x: isinstance(x, (int, float)) and not isinstance(x, bool)
    if (not isinstance(value, list) or not value
            or not all(isinstance(s, list) and len(s) == 3 and all(map(number, s)) for s in value)):
        raise ConfigError(f"{where}: expected [[q2 min, q2 max, attempts], ...], got {value!r}")
    out = []
    for lo, hi, n in value:
        if n != int(n) or n < 1:
            raise ConfigError(f"{where}: attempts must be a positive whole number, got {n!r}")
        if not 0 <= lo < hi < 1000:
            raise ConfigError(f"{where}: need 0 <= Q^2 min < Q^2 max, got {lo:g}:{hi:g}")
        out.append([float(lo), float(hi), int(n)])
    out.sort()
    for a, b in zip(out, out[1:]):
        if b[0] < a[1]:
            raise ConfigError(f"{where}: {a[0]:g}:{a[1]:g} and {b[0]:g}:{b[1]:g} overlap; slices add, "
                              f"so they must be disjoint")
    return out


def eic_software_defaults():
    """What an empty software path of [eic] resolves to on this machine."""
    home = Path.home()
    shell = os.environ.get("EIC_SHELL") or (str(home / "eic/eic-shell") if (home / "eic/eic-shell").exists()
                                            else shutil.which("eic-shell") or str(home / "eic/eic-shell"))
    return {"epic_install": os.environ.get("EPIC_INSTALL") or "/opt/detector/epic-main",
            "eicrecon_install": os.environ.get("EICRECON_INSTALL") or "/opt/local",
            "eic_shell": shell,
            "fieldmaps": os.environ.get("FIELDMAPS") or str(home / "eic/fieldmaps")}


def _validate_eic(cfg):
    run, eic = cfg["run"], cfg["eic"]
    if len(run["stages"]) > 1:
        raise ConfigError("[run] stages: eic is a chain of its own and does not combine with generate or fsi")
    if run["seed"] > schema.EIC_MAX_SEED:
        raise ConfigError(f"[run] seed must be at most {schema.EIC_MAX_SEED} for the eic stage "
                          f"(npsim is seeded with seed * 100000 + chunk)")
    if not eic["beams"]:
        raise ConfigError("[eic] beams: choose at least one beam setting")
    eic["beams"] = [b for b in schema.EIC_BEAMS if b in eic["beams"]]
    sample = schema.EIC_SAMPLES[eic["sample"]]
    if eic["ion"] not in sample.ions:
        raise ConfigError(f"[eic] ion: {eic['ion']} does not go with this sample "
                          f"(allowed: {', '.join(sample.ions)})")
    if eic["deexcite"] == "auto":
        eic["deexcite"] = "geant4" if "geant4" in sample.deexcite else "none"
    elif eic["deexcite"] not in sample.deexcite:
        raise ConfigError(f"[eic] deexcite: {eic['deexcite']} is not available for this sample "
                          f"(allowed: {', '.join(sample.deexcite)})")
    cpus = os.cpu_count() or 1
    for key in ("sim_jobs", "sim_threads", "reco_jobs", "chunk_size"):
        if eic[key] < 1:
            raise ConfigError(f"[eic] {key} must be at least 1")
    if eic["sim_jobs"] * eic["sim_threads"] > cpus:
        raise ConfigError(f"[eic] sim_jobs: {eic['sim_jobs']} x {eic['sim_threads']} threads is more than "
                          f"the {cpus} CPUs of this machine")
    if eic["reco_jobs"] > cpus:
        raise ConfigError(f"[eic] reco_jobs must be at most {cpus} (the CPUs of this machine)")
    if eic["compress"] < 0:
        raise ConfigError("[eic] compress must be >= 0")
    if eic["lumi_fb"] <= 0:
        raise ConfigError("[eic] lumi_fb must be positive")
    if eic["gen_args"] and not re.fullmatch(r"[-\w .:=+]*", eic["gen_args"]):
        raise ConfigError(f"[eic] gen_args: only generator flags and numbers, got {eic['gen_args']!r}")
    for key, value in eic_software_defaults().items():
        eic[key] = str(Path(eic[key]).expanduser()) if eic[key] else value


def validate(raw, base_dir="."):
    """Check a parsed config and return it with all defaults filled in."""
    unknown = set(raw) - set(schema.SECTIONS)
    if unknown:
        raise ConfigError(f"unknown table(s) {sorted(unknown)}; allowed: {sorted(schema.SECTIONS)}")
    cfg = {}
    for section, opts in schema.SECTIONS.items():
        given = raw.get(section, {})
        if not isinstance(given, dict):
            raise ConfigError(f"[{section}] must be a table")
        known = {o.key: o for o in opts}
        extra = set(given) - set(known)
        if extra:
            raise ConfigError(f"[{section}]: unknown option(s) {sorted(extra)}; "
                              f"allowed: {sorted(known)}")
        cfg[section] = {}
        for o in opts:
            if o.key in given:
                cfg[section][o.key] = _coerce(o, given[o.key], f"[{section}] {o.key}")
            else:
                cfg[section][o.key] = o.default
    run, gen, fsi = cfg["run"], cfg["generate"], cfg["fsi"]

    if not run["stages"]:
        raise ConfigError("[run] stages is empty")
    run["stages"] = [s for s in ["generate", "fsi", "eic"] if s in run["stages"]]
    max_jobs = max(1, os.cpu_count() or 1)
    if not 1 <= run["jobs"] <= max_jobs:
        raise ConfigError(f"[run] jobs must be between 1 and {max_jobs} (the CPUs of this machine)")
    if run["seed"] < 0:
        raise ConfigError("[run] seed must be >= 0")
    if not run["name"] or "/" in run["name"] or not run["name"].isprintable():
        raise ConfigError("[run] name must be a non-empty file name without control characters")
    out = Path(run["output_dir"]).expanduser()
    run["output_dir"] = str(out if out.is_absolute() else (Path(base_dir) / out).resolve())

    if "eic" in run["stages"]:
        _validate_eic(cfg)
    if "generate" in run["stages"]:
        if gen["beam"] <= 0:
            raise ConfigError("[generate] beam must be positive")
        if gen["events"] < 1:
            raise ConfigError("[generate] events must be >= 1")
        model = gen["model"]
        supported = [(a, z) for a, z, _ in schema.TARGETS[model]]
        if tuple(gen["target"]) not in supported:
            names = ", ".join(n for _, _, n in schema.TARGETS[model])
            raise ConfigError(f"[generate] target: A={gen['target'][0]}, Z={gen['target'][1]} is not "
                              f"in the {model} tables; supported: {names}")
        given = raw.get("generate", {})
        for o in schema.GENERATE:
            if o.key in given and model not in o.models:
                raise ConfigError(f"[generate] {o.key} applies to model(s) "
                                  f"{', '.join(o.models)}, not {model}")
    if "fsi" in run["stages"]:
        applicable = {c[0]: (c[2] if fsi["model"] == "hN" else c[3]) for c in schema.CHANNELS}
        for c in fsi["off"]:
            if not applicable[c]:
                raise ConfigError(f"[fsi] off: '{c}' is not a separate fate in the {fsi['model']} "
                                  f"cascade (see gcf_fsi --list-channels)")
        if fsi["model"] == "hA" and {"nucleon.inelastic", "nucleon.cex", "nucleon.absorption",
                                     "nucleon.piprod"} <= set(fsi["off"]):
            raise ConfigError("[fsi] off: in hA at least one of nucleon.inelastic/cex/absorption/"
                              "piprod must stay on")
        for key in ("nucleon_mfp_scale", "pion_mfp_scale"):
            if fsi[key] is not None and fsi[key] <= 0:
                raise ConfigError(f"[fsi] {key} must be positive")
    if "fsi" in run["stages"] and "generate" not in run["stages"]:
        if not fsi["input"]:
            raise ConfigError("[fsi] input is required when the generate stage is not run")
        inp = Path(fsi["input"]).expanduser()
        fsi["input"] = str(inp if inp.is_absolute() else (Path(base_dir) / inp).resolve())
    elif fsi["input"] and "generate" in run["stages"]:
        raise ConfigError("[fsi] input must not be set when the generate stage runs "
                          "(FSI then runs on the generated events)")
    return cfg


def load(path):
    path = Path(path)
    try:
        with open(path, "rb") as f:
            raw = tomllib.load(f)
    except tomllib.TOMLDecodeError as e:
        raise ConfigError(f"{path}: {e}") from None
    return validate(raw, base_dir=path.resolve().parent)


def _toml_value(v):
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, (int, float)):
        return repr(v)
    if isinstance(v, str):
        out = v.replace("\\", "\\\\").replace('"', '\\"')
        out = "".join(c if c.isprintable() or c == " " else "\\u%04X" % ord(c) for c in out)
        return '"' + out + '"'
    if isinstance(v, list):
        return "[" + ", ".join(_toml_value(x) for x in v) + "]"
    raise TypeError(f"cannot write {v!r} to TOML")


def dumps(cfg, comment=None):
    """Write a (validated) config as TOML; unset (None) options are omitted."""
    lines = []
    if comment:
        lines += [f"# {line}" for line in comment.splitlines()] + [""]
    stages = cfg["run"]["stages"]
    for section, opts in schema.SECTIONS.items():
        if section != "run" and section not in stages:
            continue
        lines.append(f"[{section}]")
        model = cfg.get("generate", {}).get("model")
        for o in opts:
            v = cfg[section].get(o.key)
            if v is None:
                continue
            if section == "generate" and model not in o.models:
                continue
            if section == "run" and "eic" in stages and o.key in ("jobs", "keep_chunks"):
                continue           # the chain has its own job counts and keeps its chunks
            lines.append(f"{o.key} = {_toml_value(v)}")
        lines.append("")
    return "\n".join(lines)
