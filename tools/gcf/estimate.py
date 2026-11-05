"""What a run will cost, measured with a short test of the real programs.

Runs a small number of events through the stages of a config and scales the
measured rates and file sizes up. Nothing it writes is kept.
"""

import math
import shutil
import subprocess
import tempfile
import time
from pathlib import Path

from . import runner

PILOT_EVENTS = 2000
PILOT_FSI_EVENTS = 300
PILOT_PROGRAM_ATTEMPTS = 2000
PILOT_TIMEOUT = 60        # seconds per stage; a slower setup is scaled from what it reached
# For a quick look at the disk space without running anything (2N and 3N, measured).
ROUGH_BYTES_PER_EVENT = {"generate": 200, "fsi": 480}


def timed_run(args, env):
    """Run a stage; returns (events/s between its first and last progress line, events done)."""
    proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, errors="replace", env=env)
    first = last = None
    deadline = time.time() + PILOT_TIMEOUT
    timed_out = False
    try:
        for line in proc.stdout:
            if line.startswith("Progress:"):
                try:
                    n = int(line.split()[1].split("/")[0])
                except (IndexError, ValueError):
                    continue
                now = time.time()
                if first is None:
                    first = (now, n)
                last = (now, n)
            if time.time() > deadline:
                timed_out = True
                proc.kill()
                break
    finally:
        proc.wait()
    if not timed_out and proc.returncode != 0:
        raise runner.RunError(f"the test run failed (exit {proc.returncode})")
    if first is None or last[0] <= first[0] or last[1] <= first[1]:
        return None, last[1] if last else 0
    return (last[1] - first[1]) / (last[0] - first[0]), last[1]


def pilot_statistics(path, tree="events"):
    """(events stored, N_eff per event) of a test's file; (None, None) without uproot."""
    try:
        import uproot
        with uproot.open(path) as f:
            w = f[tree]["weight"].array(library="np")
    except Exception:  # noqa: BLE001 - only an extra in the estimate
        return None, None
    fraction = float(w.sum()**2 / (w**2).sum() / len(w)) if len(w) and (w**2).sum() > 0 else None
    return len(w), fraction


def estimate(cfg, runs_dir):
    """Rates, time and disk space of a validated config."""
    run, gen, fsi = cfg["run"], cfg["generate"], cfg["fsi"]
    env = runner.environment()
    out = {"pilot_events": 0, "rates": {}, "bytes_per_event": {}, "n_eff_fraction": None}
    with tempfile.TemporaryDirectory(prefix="gcf-estimate-") as tmp:
        tmp = Path(tmp)
        if "program" in run["stages"]:
            # They print no progress below 100000 attempts, so time the whole test.
            events = cfg["program"]["events"]
            result = tmp / "gen.root"
            started = time.time()
            subprocess.run(runner.program_args(cfg["program"], PILOT_PROGRAM_ATTEMPTS, result, 1), capture_output=True, env=env,
                           timeout=PILOT_TIMEOUT, stdin=subprocess.DEVNULL, check=True)
            out["pilot_events"] = PILOT_PROGRAM_ATTEMPTS
            out["rates"]["program"] = PILOT_PROGRAM_ATTEMPTS / max(time.time() - started, 1e-3)
            if result.is_file():
                out["bytes_per_event"]["program"] = result.stat().st_size / PILOT_PROGRAM_ATTEMPTS
                out["pilot_stored"], out["n_eff_fraction"] = pilot_statistics(result, "genT")
            source = result
        elif "generate" in run["stages"]:
            events = gen["events"]
            pwia = tmp / "pwia.root"
            rate, done = timed_run(runner.generate_args(gen, PILOT_EVENTS, pwia, 1, 1), env)
            out["pilot_events"] = done
            out["rates"]["generate"] = rate
            if pwia.is_file() and done:
                out["bytes_per_event"]["generate"] = pwia.stat().st_size / done
                out["pilot_stored"], out["n_eff_fraction"] = pilot_statistics(pwia)
            source = pwia
        else:
            events = runner.count_entries(fsi["input"], fsi.get("tree"))
            source = Path(fsi["input"])
        if "fsi" in run["stages"] and source.is_file():
            result = tmp / "fsi.root"
            rate, done = timed_run(runner.fsi_args(fsi, source, result, 2, 0, PILOT_FSI_EVENTS), env)
            out["rates"]["fsi"] = rate
            if result.is_file() and done:
                out["bytes_per_event"]["fsi"] = result.stat().st_size / done
    out["events"] = events
    known = [r for r in out["rates"].values() if r]
    if len(known) == len(out["rates"]) and known:
        out["seconds"] = sum(events / r for r in known) / run["jobs"]
    else:
        out["seconds"] = None
    final = sum(out["bytes_per_event"].values()) * events
    out["bytes"] = final
    # Several jobs: the per-job files and the merged file exist together for a while.
    out["peak_bytes"] = 2 * final if runner.checkpoint_jobs(run, events) > 1 else final
    out["free_bytes"] = shutil.disk_usage(runs_dir).free
    if out["n_eff_fraction"] is not None:
        out["n_eff_at_most"] = out["n_eff_fraction"] * events
    return out


def rough_peak_bytes(cfg):
    """Disk space a run needs, without running anything; for refusing the impossible."""
    run, gen = cfg["run"], cfg["generate"]
    if "generate" not in run["stages"]:     # the other generators' sizes vary too much to guess
        return 0
    per_event = sum(ROUGH_BYTES_PER_EVENT[s] for s in run["stages"] if s in ROUGH_BYTES_PER_EVENT)
    return per_event * gen["events"] * (2 if runner.checkpoint_jobs(run, gen["events"]) > 1 else 1)


def human_bytes(n):
    for unit in ("B", "kB", "MB", "GB", "TB", "PB"):
        if n < 1000 or unit == "PB":
            return f"{n:.3g} {unit}"
        n /= 1000
    return ""


def human_time(seconds):
    if seconds < 90:
        return f"{math.ceil(seconds)} s"
    if seconds < 5400:
        return f"{seconds / 60:.0f} min"
    if seconds < 172800:
        return f"{seconds / 3600:.1f} h"
    return f"{seconds / 86400:.1f} days"
