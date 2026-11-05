"""Parameter scans: one setup, run for several values of one or more options.

Each point of a scan is an ordinary run with its own directory and config, so
anything that works on a run works on a point. All points share a base seed. For
an option that only changes the weights, such as a contact, they then hold the
same events, which is what makes their shapes directly comparable.
"""

import copy
import itertools
import re
import secrets

from . import schema

MAX_POINTS = 500
# Options that make no sense to scan, or whose scan is meaningless.
NOT_SCANNABLE = {("generate", "model"), ("generate", "target"), ("generate", "phase_space"), ("fsi", "input"),
                 ("fsi", "keep"), ("fsi", "tree"), ("generate", "events")}
NAME_UNSAFE = re.compile(r"[^A-Za-z0-9_.+-]")


def option(path):
    """(section, key, Opt, channel) of an axis path: generate.sigma_cm, run.seed or
    fsi.scale:pion.absorption."""
    section, _, key = path.partition(".")
    if section not in schema.SECTIONS:
        raise ValueError(f"unknown option {path!r}")
    channel = None
    if key.startswith("scale:"):
        channel = key[len("scale:"):]
        if section != "fsi" or channel not in {c[0] for c in schema.CHANNELS} or channel in schema.SWITCH_CHANNELS:
            raise ValueError(f"unknown channel in {path!r}")
        key = "scale"
    for o in schema.SECTIONS[section]:
        if o.key == key:
            if (section, key) in NOT_SCANNABLE or (section == "run" and key != "seed") \
                    or (o.kind not in ("float", "int", "str", "scales")):
                raise ValueError(f"{path} can't be scanned")
            return section, key, o, channel
    raise ValueError(f"unknown option {path!r}")


def parse_values(text, kind):
    """The values of an axis from what someone typed: 0.1, 0.15, 0.2 or start:stop:count."""
    if isinstance(text, list):
        items = text
    else:
        text = str(text).strip()
        m = re.fullmatch(r"([^:,\s]+):([^:,\s]+):(\d+)", text)
        if m and kind in ("float", "int", "scales"):
            try:
                start, stop, count = float(m.group(1)), float(m.group(2)), int(m.group(3))
            except ValueError:
                raise ValueError(f"expected start:stop:count, got {text!r}") from None
            if count < 2:
                raise ValueError("start:stop:count needs a count of at least 2")
            if count > MAX_POINTS:
                raise ValueError(f"at most {MAX_POINTS} values")
            items = [start + (stop - start) * i / (count - 1) for i in range(count)]
        else:
            items = [x for x in re.split(r"[,\s]+", text) if x]
    if not items:
        raise ValueError("give at least one value")
    out = []
    for x in items:
        try:
            if kind == "int":
                v = float(x)
                if v != int(v):
                    raise ValueError
                out.append(int(v))
            elif kind in ("float", "scales"):
                out.append(round(float(x), 12))
            else:
                out.append(str(x))
        except ValueError:
            raise ValueError(f"{x!r} is not a valid value") from None
    if len(set(out)) != len(out):
        raise ValueError("a value appears twice")
    return out


def axis_label(path):
    section, _, o, channel = option(path)
    if channel:
        return f"{channel} scale factor"
    return (o.label or o.key) + (f" ({o.unit})" if o.unit else "")


def expand(raw, axes):
    """The points of a scan: a list of (run name, values, raw config). `axes` is a list
    of {path, values}, and `raw` the config every point starts from."""
    if not axes:
        raise ValueError("add a parameter to scan")
    parsed = []
    for axis in axes:
        path = axis.get("path", "")
        _, _, o, _ = option(path)
        parsed.append((path, parse_values(axis.get("values", ""), o.kind)))
    if len({p for p, _ in parsed}) != len(parsed):
        raise ValueError("the same parameter is there twice")
    count = 1
    for _, values in parsed:
        count *= len(values)
    if count > MAX_POINTS:
        raise ValueError(f"{count} points; a scan can have at most {MAX_POINTS}")
    group = raw.get("run", {}).get("name", "")
    if not group:
        raise ValueError("the scan needs a name")
    seed = raw.get("run", {}).get("seed") or 1 + secrets.randbelow(2**31 - 2**21)
    points = []
    names = set()
    for combination in itertools.product(*(values for _, values in parsed)):
        cfg = copy.deepcopy(raw)
        cfg.setdefault("run", {}).update(seed=seed, scan=group)
        parts = []
        for (path, _), value in zip(parsed, combination):
            section, key, _, channel = option(path)
            target = cfg.setdefault(section, {})
            if channel:
                target.setdefault("scale", {})[channel] = value
                short = channel.replace(".", "_")
            else:
                target[key] = value
                short = key
            parts.append(f"{short}-{NAME_UNSAFE.sub('_', format(value, 'g') if isinstance(value, float) else str(value))}")
        name = f"{group}_" + "_".join(parts)
        if len(name) > 100:      # run names are limited; keep them unique
            name = f"{group[:60]}_p{len(points):03d}"
        if name in names:
            raise ValueError(f"two points would both be named {name}")
        names.add(name)
        cfg["run"]["name"] = name
        points.append((name, list(combination), cfg))
    return points
