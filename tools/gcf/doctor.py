"""Is this machine ready to use the pipeline? What is missing, and how to fix it."""

import importlib.util
import os
import re
import shutil
import subprocess
from pathlib import Path

from . import config, runner, schema

LOW_DISK_BYTES = 10 * 1000**3


def check(name, ok, message="", fix="", required=True):
    return {"name": name, "ok": bool(ok), "message": message, "fix": fix, "required": required}


def built_version():
    """The commit the stage programs were built from, if they say so."""
    header = runner.BIN_DIR / "gcf_version.h"
    try:
        m = re.search(r'"([^"]+)"', header.read_text())
    except OSError:
        return None
    return m.group(1) if m else None


def stale_binaries():
    """A message if the code under src/ is newer than the stage programs."""
    version = built_version()
    binary = runner.BIN_DIR / "gcf_generate"
    if not version or version == "unknown" or not binary.exists():
        return None
    try:
        out = subprocess.run(["git", "-C", str(runner.REPO), "diff", "--name-only", version.removesuffix("-dirty"), "--", "src",
                              ":(exclude)src/programs/genQE_FSI/Generator-R-3_06_02"],
                             capture_output=True, text=True)
    except OSError:              # no git: nothing to compare with
        return None
    if out.returncode != 0:
        return None
    built = binary.stat().st_mtime
    for name in out.stdout.split():
        path = runner.REPO / name
        # Deleted, or changed after the build. A dirty build's own changes are older than it.
        if not path.exists() or path.stat().st_mtime > built:
            return f"they were built from {version}, and the code has changed since"
    return None


def logout_kills_processes():
    """True if the system kills a user's processes at logout, False if not, None if unknown."""
    if not shutil.which("busctl"):
        return None
    out = subprocess.run(["busctl", "get-property", "org.freedesktop.login1", "/org/freedesktop/login1",
                          "org.freedesktop.login1.Manager", "KillUserProcesses"], capture_output=True, text=True)
    if out.returncode != 0:
        return None
    if out.stdout.strip() != "b true":
        return False
    linger = subprocess.run(["loginctl", "show-user", os.environ.get("USER", ""), "-p", "Linger"],
                            capture_output=True, text=True)
    return "Linger=yes" not in linger.stdout


def has_module(name):
    return importlib.util.find_spec(name) is not None


def checks(runs_dir=None):
    runs_dir = Path(runs_dir) if runs_dir else runner.default_runs_dir()
    out = []
    env = runner.environment()
    for name in ("gcf_generate", "gcf_fsi"):
        path = runner.BIN_DIR / name
        if not os.access(path, os.X_OK):
            out.append(check(name, False, f"not built ({path} is missing)", "./install.sh"))
            continue
        run = subprocess.run([str(path), "--help"], capture_output=True, text=True, env=env, timeout=60)
        out.append(check(name, run.returncode == 0,
                         "" if run.returncode == 0 else f"does not start: {(run.stderr or run.stdout).strip()[-200:]}",
                         "" if run.returncode == 0 else "./install.sh, or check build/gcf_env.sh"))
    stale = stale_binaries()
    out.append(check("stage programs match the code", not stale, stale or "", "./install.sh" if stale else "", required=False))
    missing_programs = []
    for name in schema.PROGRAMS:
        try:
            runner.program_binary(name)
        except runner.RunError:
            missing_programs.append(name)
    # Worth a line here, not a banner: only someone who uses them needs them.
    out.append(check("other generators", True, "not built: " + ", ".join(missing_programs) + "; ./install.sh builds them" if missing_programs else ""))
    hadd = shutil.which("hadd", path=env.get("PATH"))
    out.append(check("hadd (ROOT)", hadd, "" if hadd else "not found", "" if hadd else "./install.sh; it needs ROOT"))
    missing = [m for m in ("uproot", "awkward", "numpy", "matplotlib") if not has_module(m)]
    out.append(check("Python packages for plots and summaries", not missing,
                     "missing: " + ", ".join(missing) if missing else "",
                     "pip install -r tools/gcf/requirements.txt" if missing else "", required=False))
    scipy = has_module("scipy")
    out.append(check("Python package for comparing runs", scipy, "" if scipy else "missing: scipy",
                     "" if scipy else "pip install -r tools/gcf/requirements.txt", required=False))
    shell = config.eic_software_defaults()["eic_shell"]
    out.append(check("eic-shell, for the ePIC chain", True, "" if Path(shell).exists() else f"not found at {shell}; only the EIC / ePIC stage needs it"))
    base = runs_dir if runs_dir.exists() else next((p for p in runs_dir.parents if p.exists()), Path("/"))
    writable = os.access(base, os.W_OK)
    out.append(check(f"runs directory {runs_dir}", writable, "" if writable else "not writable",
                     "" if writable else "choose another with --runs-dir"))
    others_write = base.stat().st_uid == os.getuid() and runner.others_can_write(base)
    out.append(check("nobody else can write to the runs directory", not others_write,
                     f"others can write to {base}" if others_write else "",
                     f"chmod go-w {base}" if others_write else "", required=False))
    free = shutil.disk_usage(base).free
    out.append(check("disk space", free >= LOW_DISK_BYTES, f"{free / 1e9:.0f} GB free in {base}",
                     "" if free >= LOW_DISK_BYTES else "choose another --runs-dir, or free some space", required=False))
    kills = logout_kills_processes()
    out.append(check("runs survive logout", kills is not True,
                     "this system stops your processes when you log out" if kills else
                     ("could not tell" if kills is None else ""),
                     "ask the administrator for: loginctl enable-linger $USER" if kills else "", required=False))
    return out


def text(results):
    lines = []
    for c in results:
        mark = "ok  " if c["ok"] else ("FAIL" if c["required"] else "warn")
        lines.append(f"{mark} {c['name']}" + (f": {c['message']}" if c["message"] else ""))
        if c["fix"] and not c["ok"]:
            lines.append(f"     fix: {c['fix']}")
    return "\n".join(lines)
