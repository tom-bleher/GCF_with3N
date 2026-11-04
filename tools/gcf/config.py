"""Run configurations: TOML in, validated dict out, and back to TOML.

A config has up to three tables, [run], [generate] and [fsi]; their keys are
listed in schema.py. Validation fills in defaults, so a saved (resolved)
config records every setting a run used.
"""

import math
import os
import tomllib
from pathlib import Path

from . import schema


class ConfigError(Exception):
    pass


def _check_limits(opt, value, where):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return
    unit = f" {opt.unit}" if opt.unit else ""
    if opt.ge is not None and value < opt.ge:
        raise ConfigError(f"{where}: must be at least {opt.ge:g}{unit}, got {value:g}")
    if opt.gt is not None and value <= opt.gt:
        raise ConfigError(f"{where}: must be above {opt.gt:g}{unit}, got {value:g}")
    if opt.le is not None and value > opt.le:
        raise ConfigError(f"{where}: must be at most {opt.le:g}{unit}, got {value:g}")
    if isinstance(value, float) and not math.isfinite(value):
        raise ConfigError(f"{where}: must be a finite number")


def _check_combinations(gen):
    """Options that are each valid but do nothing where they are set. The binary
    refuses them too (#123); here the form hears about it as you type."""
    model, (a, _z) = gen["model"], gen["target"]
    given = {o.key: gen.get(o.key) for o in schema.GENERATE if model in o.models}

    def refuse(key, why):
        raise ConfigError(f"[generate] {key}: {why}")
    if model == "2N" and a == 2:
        for key in ("sigma_cm", "estar", "sigma_e"):
            if given.get(key) is not None:
                refuse(key, "does not apply to the deuteron: a pair is the whole nucleus, with no residual")
        for key in ("cpp0", "cnn0", "cpn0", "cpn1"):
            if given.get(key) is not None:
                refuse(key, "does not apply to the deuteron: it uses a momentum distribution, not contacts")
        if given.get("randomize"):
            refuse("randomize", "does not apply to the deuteron: it uses a momentum distribution, not contacts")
    if model == "3N" and a == 3:
        if given.get("sigma_cm") is not None:
            refuse("sigma_cm", "does not apply to A = 3: the triplet is the whole nucleus, at rest")
        if given.get("estar") is not None:
            refuse("estar", "does not apply to A = 3: a triplet leaves no residual")
    if model == "2N" and given.get("onset") == "local" and given.get("kcut") is not None:
        refuse("kcut", "has no effect with the local SRC onset, which replaces the cut with k_F(r)")
    if model == "3N" and given.get("cm") is False and given.get("sigma_cm") is not None:
        refuse("sigma_cm", "has no effect with the CM motion switched off")
    if model == "MF" and given.get("fermi_gas") == "local" and given.get("kF_fg") is not None:
        refuse("kF_fg", "has no effect with the local Fermi gas, which takes k_F(r) from the density")
    if model == "2N" and given.get("light_cone"):
        if given.get("prel_is"):
            refuse("prel_is", "does not apply with the light-cone cross section")
        if given.get("onset") == "local":
            refuse("onset", "the light-cone cross section has no local onset")
    if model == "2N" and any(given.get(k) == 0 for k in ("cpp0", "cnn0", "cpn0", "cpn1")) \
            and all(given.get(k) == 0 for k in ("cpp0", "cnn0", "cpn0", "cpn1")):
        refuse("cpn1", "all contacts are zero, so every event would have zero weight")


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
        for end in value:
            _check_limits(opt, end, where)
    elif kind == "target":
        if (not isinstance(value, list) or len(value) != 2
                or not all(isinstance(v, int) and not isinstance(v, bool) for v in value)):
            raise ConfigError(f"{where}: expected [A, Z], got {value!r}")
        if value[0] < 1 or not 0 <= value[1] <= value[0]:
            raise ConfigError(f"{where}: invalid nucleus A={value[0]}, Z={value[1]}")
    elif kind == "list":
        if not isinstance(value, list) or not all(isinstance(v, str) for v in value):
            raise ConfigError(f"{where}: expected a list of strings, got {value!r}")
    elif kind == "scales":
        names = {c[0] for c in schema.CHANNELS}
        if not isinstance(value, dict):
            raise ConfigError(f"{where}: expected a table of channel = factor, got {value!r}")
        out = {}
        for name, factor in value.items():
            if name not in names:
                raise ConfigError(f"{where}: unknown channel {name!r}; allowed: {sorted(names)}")
            if name in schema.SWITCH_CHANNELS:
                raise ConfigError(f"{where}: {name} can only be switched off, not scaled")
            if isinstance(factor, bool) or not isinstance(factor, (int, float)) or not factor >= 0:
                raise ConfigError(f"{where}: the factor of {name} must be a number >= 0, got {factor!r}")
            out[name] = float(factor)
        value = out
    if kind in ("float", "int"):
        _check_limits(opt, value, where)
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
    run["stages"] = [s for s in ["generate", "fsi", "program"] if s in run["stages"]]
    max_jobs = max(1, os.cpu_count() or 1)
    if not 1 <= run["jobs"] <= max_jobs:
        raise ConfigError(f"[run] jobs must be between 1 and {max_jobs} (the CPUs of this machine)")
    if run["seed"] < 0:
        raise ConfigError("[run] seed must be >= 0")
    if run["chunks"] is not None and run["chunks"] < 1:
        raise ConfigError("[run] chunks must be at least 1")
    if not run["name"] or "/" in run["name"] or not run["name"].isprintable():
        raise ConfigError("[run] name must be a non-empty file name without control characters")
    out = Path(run["output_dir"]).expanduser()
    run["output_dir"] = str(out if out.is_absolute() else (Path(base_dir) / out).resolve())

    if "program" in run["stages"]:
        if len(run["stages"]) > 1:
            raise ConfigError("[run] stages: another generator runs on its own, not with generate or fsi "
                              "(it doesn't write the event-file format)")
        name = cfg["program"]["program"]
        given = raw.get("program", {})
        for o in schema.PROGRAM:
            if o.key in given and name not in o.models:
                raise ConfigError(f"[program] {o.key} does not apply to {name}")
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
        _check_combinations(gen)
    if "fsi" in run["stages"]:
        applicable = {c[0]: (c[2] if fsi["model"] == "hN" else c[3]) for c in schema.CHANNELS}
        for c in fsi["off"]:
            if not applicable[c]:
                raise ConfigError(f"[fsi] off: '{c}' is not a separate fate in the {fsi['model']} "
                                  f"cascade (see gcf_fsi --list-channels)")
        for c, factor in fsi["scale"].items():
            if not applicable[c]:
                raise ConfigError(f"[fsi] scale: '{c}' is not a separate fate in the {fsi['model']} "
                                  f"cascade (see gcf_fsi --list-channels)")
            if c in fsi["off"]:
                raise ConfigError(f"[fsi] scale: '{c}' is also switched off")
        # A factor of 0 is the same as switching the channel off: one spelling in the saved config.
        zero = [c for c, factor in fsi["scale"].items() if factor == 0]
        fsi["off"] = fsi["off"] + zero
        fsi["scale"] = {c: factor for c, factor in fsi["scale"].items() if factor != 0}
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
    if isinstance(v, dict):
        return "{" + ", ".join(f"{_toml_value(k)} = {_toml_value(x)}" for k, x in v.items()) + "}"
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
            if section == "program" and cfg["program"]["program"] not in o.models:
                continue
            lines.append(f"{o.key} = {_toml_value(v)}")
        lines.append("")
    return "\n".join(lines)
