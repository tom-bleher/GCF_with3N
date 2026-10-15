"""Run configurations: TOML in, validated dict out, and back to TOML.

A config has up to three tables, [run], [generate] and [fsi]; their keys are
listed in schema.py. Validation fills in defaults, so a saved (resolved)
config records every setting a run used.
"""

import os
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
    elif kind == "list":
        if not isinstance(value, list) or not all(isinstance(v, str) for v in value):
            raise ConfigError(f"{where}: expected a list of strings, got {value!r}")
    if opt.choices is not None:
        vals = value if kind == "list" else [value]
        for v in vals:
            if v not in opt.choices:
                raise ConfigError(f"{where}: {v!r} is not one of {opt.choices}")
    return value


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
    run["stages"] = [s for s in ["generate", "fsi"] if s in run["stages"]]
    max_jobs = max(1, os.cpu_count() or 1)
    if not 1 <= run["jobs"] <= max_jobs:
        raise ConfigError(f"[run] jobs must be between 1 and {max_jobs} (the CPUs of this machine)")
    if run["seed"] < 0:
        raise ConfigError("[run] seed must be >= 0")
    if not run["name"] or "/" in run["name"] or not run["name"].isprintable():
        raise ConfigError("[run] name must be a non-empty file name without control characters")
    out = Path(run["output_dir"]).expanduser()
    run["output_dir"] = str(out if out.is_absolute() else (Path(base_dir) / out).resolve())

    if "generate" in run["stages"]:
        if gen["beam"] <= 0:
            raise ConfigError("[generate] beam must be positive")
        if gen["events"] < 1:
            raise ConfigError("[generate] events must be >= 1")
        model = gen["model"]
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
            lines.append(f"{o.key} = {_toml_value(v)}")
        lines.append("")
    return "\n".join(lines)
