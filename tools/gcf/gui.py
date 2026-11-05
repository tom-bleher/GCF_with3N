"""`gcf gui`: a small web interface for the pipeline.

Runs on the analysis machine and is opened in a browser through an SSH
tunnel. It only builds a config, starts `gcf run` with it, and shows what
that run writes (status.json, logs, outputs, quick-look plots), so anything
done here can be repeated from the command line with the saved config.

Standard library only. Listens on localhost; every request must carry the
random token printed at start-up, because other users on the same machine
can reach localhost ports too.
"""

import json
import mimetypes
import os
import re
import shutil
import tomllib
import secrets
import signal
import socket
import subprocess
import sys
import threading
import time
import urllib.parse
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

from . import cli, config, doctor, estimate, results, runner, scan, schema

STATIC = Path(__file__).resolve().parent / "static"
MODULE_RE = re.compile(r"[a-z_]+\.js\Z")
COMPARE_DIR = ".compare"       # hidden: not a run
SCANS_DIR = ".scans"           # hidden: what makes runs a scan
COMPARE_FILE_RE = re.compile(r"[A-Za-z0-9_.+-]+\.png\Z")
COMPARE_SCRIPT = runner.REPO / "tools/validation/compare_samples.py"
FIRST_PORT = 8765
LOG_TAIL_BYTES = 64 * 1024
MAX_LOGS_SHOWN = 12            # the first logs of a run with many jobs
SCHEDULE_SECONDS = 2.0
MAX_AUTO_RESUMES = 3           # a run that keeps dying waits for you after this
RUNNING_STATES = ("running", "merging", "starting")


def tail_lines(path, count, keep=None):
    """The last lines of a file, without reading all of it."""
    with open(path, "rb") as f:
        f.seek(0, os.SEEK_END)
        f.seek(max(0, f.tell() - LOG_TAIL_BYTES))
        lines = f.read().decode(errors="replace").splitlines()
    if keep:
        lines = [line for line in lines if keep(line)]
    return lines[-count:]


def safe_to_run(d):
    """What is in a run's directory gets run as you, including paths in its config. So
    only if the directory and its configs are yours, and nobody else can write to them."""
    try:
        for path in (d, d / "config.toml", d / "submitted.toml"):
            if not path.exists():
                continue
            if path.stat().st_uid != os.getuid() or runner.others_can_write(path):
                return False
        return True
    except OSError:
        return False


def config_dir():
    return Path(os.environ.get("XDG_CONFIG_HOME") or Path.home() / ".config") / "gcf"


def stable_token(new=False):
    """The page's token, kept in a file only you can read, so that a bookmark goes on
    working after a restart. A random one for this start if that file can't be written."""
    path = token_path()
    if not new:
        try:
            token = path.read_text().strip()
            if token and path.stat().st_mode & 0o077 == 0:
                return token
        except OSError:
            pass
    token = secrets.token_urlsafe(18)
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.unlink(missing_ok=True)
        fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
        with os.fdopen(fd, "w") as f:
            f.write(token + "\n")
    except OSError:
        pass
    return token


def token_path():
    return config_dir() / "gui-token"


def known_dirs_path():
    return config_dir() / "gui-known-dirs.json"


def server_record(port):
    return config_dir() / f"gui-{port}.json"


def free_port(first=FIRST_PORT, count=50):
    """The first port from `first` that nothing listens on: others may run this too."""
    for port in range(first, first + count):
        with socket.socket() as probe:
            probe.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
            try:
                probe.bind(("127.0.0.1", port))
                return port
            except OSError:
                continue
    raise runner.RunError(f"no free port from {first} to {first + count - 1}; choose one with --port")


def my_servers(port=None):
    """The servers of yours that are running, on one port or all."""
    if port:
        rec = read_server(port)
        return [rec] if rec else []
    found = []
    for path in sorted(config_dir().glob("gui-*.json")):
        try:
            rec = read_server(int(path.stem.split("-")[1]))
        except ValueError:
            continue
        if rec:
            found.append(rec)
    return found


def read_server(port):
    """What is known of a server of yours on this port; None if it isn't running."""
    try:
        rec = json.loads(server_record(port).read_text())
    except (OSError, ValueError):
        return None
    return rec if runner.proc_is(rec.get("pid"), rec.get("pid_start")) else None
NAME_RE = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.+-]{0,99}\Z")
PLACEHOLDER_STATE = "failed"


def valid_name(name):
    return bool(name) and NAME_RE.fullmatch(name) is not None


class App:
    def __init__(self, runs_dir, serve=True, token=None, max_cpus=None):
        self.runs_dir = Path(runs_dir).expanduser().resolve()
        self.known_dirs = [self.runs_dir]   # the scheduler looks after all of these
        if serve:
            self.runs_dir.mkdir(parents=True, exist_ok=True)
            self.known_dirs += [p for p in self.saved_dirs() if p != self.runs_dir]
        self.fixed_token = token or secrets.token_urlsafe(18)
        self.follows_token_file = False     # set by serve(), if it uses the one in the file
        self.max_cpus = max_cpus or os.cpu_count() or 1
        self.lock = threading.Lock()
        self.schedule_lock = threading.Lock()
        self.children = {}            # run directory -> Popen of `gcf run`
        self.plot_lock = threading.Lock()   # matplotlib is not thread-safe
        self.results_lock = threading.Lock()
        self.nuclear_cache = {}             # (model, A, Z, nn) -> the generator's table values

    @property
    def token(self):
        """The token in force: a running server follows gcf gui --new-token."""
        if self.follows_token_file:
            try:
                current = token_path().read_text().strip()
                if current and token_path().stat().st_mode & 0o077 == 0:
                    return current
            except OSError:
                pass
        return self.fixed_token

    @staticmethod
    def saved_dirs():
        try:
            return [Path(p) for p in json.loads(known_dirs_path().read_text()) if Path(p).is_dir()]
        except (OSError, ValueError, TypeError):
            return []

    def save_dirs(self):
        try:
            known_dirs_path().parent.mkdir(parents=True, exist_ok=True)
            known_dirs_path().write_text(json.dumps([str(p) for p in self.known_dirs]))
        except OSError:
            pass

    # -- runs ----------------------------------------------------------------
    def run_dir(self, name):
        if not valid_name(name):
            raise ValueError(f"invalid run name {name!r}")
        return self.runs_dir / name

    @staticmethod
    def live_steps(st):
        """Stage processes of a run that are still the recorded processes."""
        return [s for steps in st.get("jobs", []) for s in steps
                if s.get("pid") and runner.proc_is(s["pid"], s.get("pid_start"))]

    def run_summary(self, d):
        st = runner.read_status(d)
        state = st.get("state", "starting" if (d / "submitted.toml").exists() else "unknown")
        with self.lock:
            child = self.children.get(str(d))
            if child is not None and child.poll() is not None:
                del self.children[str(d)]   # reap and forget finished supervisors
                child = None
        # The supervisor is only "running" if its pid is still the recorded
        # process (pid + start time); a reused pid after a crash or reboot is
        # not ours and must never be signalled.
        supervisor_alive = runner.proc_is(st.get("pid"), st.get("pid_start"))
        orphans = []
        if state in ("running", "merging") and not supervisor_alive:
            state = "stopped"          # the supervisor is gone (killed, machine reboot)
            orphans = self.live_steps(st)
        if (d / "queued").exists() or (d / ".starting").exists():
            state = "queued"
        if state == "starting" and child is None and not runner.is_locked(d):
            # No supervisor we know of and no status.json: `gcf run` died
            # before writing one (bad config, missing hadd, ...).
            state = "failed"
            st["error"] = (d / "gcf.log").read_text(errors="replace")[-2000:].strip() \
                if (d / "gcf.log").exists() else "gcf run exited before starting"
            if not st["error"]:
                st["error"] = "gcf run exited before starting"
        outputs = {}
        for stage, path in st.get("outputs", {}).items():
            p = Path(path)
            png = p.with_suffix(".png")
            outputs[stage] = {"path": path, "exists": p.is_file(),
                              "size": p.stat().st_size if p.is_file() else 0,
                              # A plot older than its file shows another sample.
                              "plot": png.is_file() and p.is_file() and png.stat().st_mtime >= p.stat().st_mtime}
        jobs = st.get("jobs", [])
        finished = sum(all(step.get("state") == "done" for step in job) for job in jobs)
        stopped = state in ("stopped", "failed", "cancelled", "partial")
        has_config = (d / "config.toml").exists()
        any_step_done = any(step.get("state") == "done" for job in jobs for step in job)
        return {"name": d.name, "state": state, "error": st.get("error"), "hint": runner.failure_hint(st.get("error")),
                "started": st.get("started"), "finished": st.get("finished"),
                "stages": st.get("stages", {}), "outputs": outputs,
                "orphans": len(orphans), "jobs_total": len(jobs), "jobs_finished": finished,
                "workers": st.get("workers"), "partial": st.get("partial"),
                "resumable": bool(stopped and has_config and jobs),
                "mergeable": bool(state in ("stopped", "failed", "cancelled") and any_step_done
                                  and {"generate", "program"} & set(st.get("stages", {}))
                                  and "simulate" not in st.get("stages", {})),
                "note": self.read_config(d).get("run", {}).get("note"),
                "scan": self.read_config(d).get("run", {}).get("scan"),
                "setup": self.setup_summary(d), "mtime": d.stat().st_mtime}

    @staticmethod
    def read_config(d):
        for name in ("config.toml", "submitted.toml"):
            try:
                return tomllib.loads((d / name).read_text())
            except (OSError, ValueError):   # unreadable, not UTF-8, or not TOML
                continue
        return {}

    def setup_summary(self, d):
        """One line saying what a run is, from its config: '3N · 12C · 6 GeV · 40k · hN'."""
        cfg = d / "config.toml"
        if not cfg.exists():
            cfg = d / "submitted.toml"
        try:
            c = tomllib.loads(cfg.read_text())
        except (OSError, ValueError):
            return ""
        parts = []
        g, f, e, r, program = c.get("generate"), c.get("fsi"), c.get("eic"), c.get("run", {}), c.get("program")
        if program:
            parts.append(program.get("program", "?"))
            if "target" in program:
                a, z = program["target"]
                parts.append(f"A={a} Z={z}")
            if "beam" in program:
                parts.append(f"{program['beam']:g} GeV")
            if "events" in program:
                n = program["events"]
                parts.append((f"{n / 1e6:g}M" if n >= 1e6 else f"{n / 1e3:g}k" if n >= 1e3 else str(n)) + " attempts")
            if r.get("jobs", 1) > 1:
                parts.append(f"{r['jobs']} jobs")
            return " · ".join(parts)
        if e:
            sample = schema.EIC_SAMPLES.get(e.get("sample"))
            slices = e.get("slices") or [[0, 0, 0]]
            attempts = sum(x[2] for x in slices)
            parts += ["ePIC", sample.label.split(" (")[0] if sample else str(e.get("sample")),
                      schema.EIC_IONS.get(e.get("ion"), str(e.get("ion"))), " + ".join(e.get("beams", [])),
                      f"Q2 {slices[0][0]:g}-{slices[-1][1]:g}",
                      (f"{attempts / 1e6:g}M" if attempts >= 1e6 else f"{attempts / 1e3:g}k") + " attempts"]
            if e.get("deexcite") not in (None, "none", "auto"):
                parts.append(e["deexcite"])
            return " · ".join(parts)
        if g:
            parts.append(g.get("model", "?"))
            if "target" in g:
                a, z = g["target"]
                names = {(2, 1): "2H", (3, 1): "3H", (3, 2): "3He", (4, 2): "4He", (6, 3): "6Li", (12, 6): "12C",
                         (14, 7): "14N", (16, 8): "16O", (27, 13): "27Al", (40, 18): "40Ar", (40, 20): "40Ca",
                         (56, 26): "56Fe", (208, 82): "208Pb"}
                parts.append(names.get((a, z), f"A={a} Z={z}"))
            if "beam" in g:
                parts.append(f"{g['beam']:g} GeV")
            if "events" in g:
                n = g["events"]
                parts.append((f"{n / 1e6:g}M" if n >= 1e6 else f"{n / 1e3:g}k" if n >= 1e3 else str(n)) + " events")
            if g.get("wave_function") and g.get("wave_function") != "AV18":
                parts.append(g["wave_function"])
        if f:
            fs = f.get("model", "hN") + " FSI"
            if f.get("off"):
                fs += " (off: " + ", ".join(x.split(".")[-1] if False else x for x in f["off"]) + ")"
            if not g and f.get("input"):
                parts.append(Path(f["input"]).name)
            parts.append(fs)
        elif "fsi" not in r.get("stages", []):
            parts.append("PWIA only")
        if r.get("jobs", 1) > 1:
            parts.append(f"{r['jobs']} jobs")
        return " · ".join(parts)

    def is_run(self, d):
        return d.is_dir() and valid_name(d.name) and \
            ((d / "status.json").exists() or (d / "submitted.toml").exists())

    def list_runs(self):
        runs = []
        if not self.runs_dir.is_dir():
            return runs
        for d in self.runs_dir.iterdir():
            try:
                if self.is_run(d):
                    runs.append(self.run_summary(d))
            except OSError:
                continue           # vanished between listing and reading
            except Exception as e:  # noqa: BLE001 - one damaged run must not hide the others
                runs.append(self.broken_run(d, e))
        return sorted(runs, key=lambda r: r["mtime"], reverse=True)

    @staticmethod
    def broken_run(d, problem):
        return {"name": d.name, "state": "unknown", "error": f"this run's files can't be read: {problem}", "hint": None,
                "started": None, "finished": None, "stages": {}, "outputs": {}, "orphans": 0, "jobs_total": 0,
                "jobs_finished": 0, "workers": None, "partial": None, "resumable": False, "mergeable": False,
                "note": None, "scan": None, "setup": "", "mtime": 0}

    def run_details(self, name):
        d = self.run_dir(name)
        if not d.is_dir():
            raise FileNotFoundError(name)
        try:
            out = self.run_summary(d)
        except OSError:
            raise
        except Exception as e:  # noqa: BLE001
            out = self.broken_run(d, e)
        jobs = runner.read_status(d).get("jobs")
        out["jobs"] = jobs if isinstance(jobs, list) else []
        cfg = d / "config.toml"
        if not cfg.exists():
            cfg = d / "submitted.toml"
        out["config"] = cfg.read_text(errors="replace") if cfg.exists() else None
        out["config_data"] = None
        if cfg.exists():
            try:
                out["config_data"] = tomllib.loads(out["config"])
            except tomllib.TOMLDecodeError:
                pass
        logs = []
        order = {"generate": 0, "fsi": 1, "simulate": 2, "reconstruct": 3, "resume": 4, "analyze": 5, "merge": 6}
        files = sorted((d / "logs").glob("*.log"), key=lambda p: (order.get(p.name.split("_")[0], 9), p.name)) \
            if (d / "logs").is_dir() else []
        # With many jobs: those that failed or are running, before the others.
        interesting = {Path(step.get("log", "")).name: step["state"] for job in out["jobs"] for step in job
                       if step.get("state") in ("failed", "running")}
        rank = {"failed": 0, "running": 1}
        shown = sorted(files, key=lambda p: rank.get(interesting.get(p.name), 2))[:MAX_LOGS_SHOWN]
        # drop progress lines and GENIE's start-up banner / messenger chatter
        keep = (lambda l: l.strip() and not l.startswith("Progress:")
                and not l.lstrip().startswith("*") and "NOTICE Messenger" not in l)
        for log in shown:
            logs.append({"name": log.name, "file": f"logs/{log.name}", "size": log.stat().st_size,
                         "state": interesting.get(log.name), "tail": "\n".join(tail_lines(log, 40, keep))})
        gcf_log = d / "gcf.log"
        if gcf_log.exists() and gcf_log.stat().st_size:
            logs.insert(0, {"name": "gcf.log", "file": "gcf.log", "size": gcf_log.stat().st_size,
                            "tail": "\n".join(tail_lines(gcf_log, 40))})
        out["logs"] = logs
        out["logs_hidden"] = max(0, len(files) - MAX_LOGS_SHOWN)
        out["logs_dir"] = str(d / "logs")
        out["size"] = sum(f.stat().st_size for f in d.rglob("*") if f.is_file())
        out["results"] = {}
        for stage, o in out["outputs"].items():
            if o["exists"]:
                cached = results.cached(o["path"], d / "results.json", stage)
                if cached:
                    out["results"][stage] = cached
        if out["config_data"] and "eic" in out["config_data"]:
            out["eic"] = self.eic_details(d, out["config_data"]["eic"])
            # while the container works the step logs are silent; each slice's run.log is not
            for s in out["eic"]["slices"]:
                run_log = d / s["label"] / "run.log"
                if run_log.is_file():
                    logs.append({"name": f"{s['label']}/run.log", "tail": "\n".join(tail_lines(run_log, 15))})
        return out

    def eic_details(self, d, eic):
        """Per slice: events and chunks; the software that made them; the figures."""
        slices, provenance = [], None
        for label, beams, lo, hi, attempts in runner.eic_slices(eic):
            s = d / label
            info = {"label": label, "beams": beams, "q2": [lo, hi], "attempts": attempts,
                    "path": str(s), "events": None, "jobs": None}
            gen = runner._read_json(s / "gen.hepmc3.json")
            if gen:
                info["events"] = gen.get("events_written")
            rec = runner._read_json(s / "run.json")
            if rec:
                info["jobs"] = rec.get("jobs")
                if provenance is None and rec.get("container"):
                    git = lambda part: (rec.get(part) or {}).get("git") or {}
                    provenance = {"gcf": " ".join(rec.get("gcf") or []), "generator": " ".join(rec.get("generator") or []),
                                  "epic": rec.get("epic", {}).get("install"), "epic_config": rec.get("epic", {}).get("config"),
                                  "epic_commit": git("epic").get("commit"), "eicrecon": rec.get("eicrecon", {}).get("install"),
                                  "eicrecon_commit": git("eicrecon").get("commit"), "abconv_preset": rec.get("abconv_preset"),
                                  "container_image": rec["container"].get("container_image"),
                                  "container_from": rec["container"].get("container_from"),
                                  "container_build_date": rec["container"].get("container_build_date"),
                                  "detector_path": rec["container"].get("detector_path")}
            slices.append(info)
        figures, plots = [], d / "plots"
        if plots.is_dir():
            figures = [{"name": p.stem, "png": p.name, "pdf": p.with_suffix(".pdf").name, "mtime": p.stat().st_mtime}
                       for p in sorted(plots.glob("*.png"))]
        return {"slices": slices, "provenance": provenance, "figures": figures,
                "summary": (plots / "summary.json").is_file(), "plots": str(plots)}

    @staticmethod
    def with_runs_dir(raw, runs_dir):
        """A copy of a request config with [run] output_dir forced to runs_dir."""
        if not isinstance(raw, dict) or any(not isinstance(v, dict) for v in raw.values()):
            raise ValueError("the config must be an object of tables ([run], [generate], [fsi])")
        raw = dict(raw)
        raw["run"] = dict(raw.get("run", {}), output_dir=str(runs_dir))
        return raw

    def check(self, raw):
        """Validate a config and run the stages' setup checks; returns the
        resolved TOML, a list of problems and, for the ePIC chain, what it will cost."""
        raw = self.with_runs_dir(raw, self.runs_dir)
        cfg = config.validate(raw, base_dir=self.runs_dir)
        problems = []
        name = cfg["run"]["name"]
        if not valid_name(name):
            problems.append(f"invalid run name {name!r}: letters, digits, . _ + - only")
        elif (self.runs_dir / name).exists():
            problems.append(f"a run named {name!r} already exists; choose another name")
        problems += runner.check_setup(cfg)
        problems += self.resource_problems(cfg)
        cost = None
        eic = "eic" in cfg["run"]["stages"]
        if not problems and not eic:
            try:
                cost = estimate.estimate(cfg, self.runs_dir)
            except (runner.RunError, OSError, subprocess.SubprocessError) as e:
                print(f"gcf gui: no estimate: {e}", file=sys.stderr)
        if cost and cost["peak_bytes"] > cost["free_bytes"]:
            problems.append(f"this run needs about {estimate.human_bytes(cost['peak_bytes'])} of disk space, and "
                            f"{self.runs_dir} has {estimate.human_bytes(cost['free_bytes'])} free")
        out = {"toml": config.dumps(cfg, comment="Run with:  gcf run <this file>"),
               "problems": problems, "ok": not problems, "estimate": cost,
               "estimate_text": estimate_text(cost, cfg) if cost else None}
        if eic:
            out["cost"] = dict(zip(("events", "gb", "core_hours"), runner.eic_cost(cfg)))
        return out

    def resource_problems(self, cfg):
        problems = []
        if "eic" in cfg["run"]["stages"]:
            needed = self.cpus_of_config(cfg)
            if needed > self.max_cpus:
                problems.append(f"[eic] sim_jobs: this chain needs {needed} CPUs at a time, more than the {self.max_cpus} "
                                f"this interface may use (gcf gui --max-cpus)")
            return problems
        if cfg["run"]["jobs"] > self.max_cpus:
            problems.append(f"[run] jobs: {cfg['run']['jobs']} is more than the {self.max_cpus} CPUs "
                            f"this interface may use (gcf gui --max-cpus)")
        need, free = estimate.rough_peak_bytes(cfg), shutil.disk_usage(self.runs_dir).free
        if need / 2 > free:       # a rough number, so only for what clearly can't fit
            problems.append(f"[generate] events: this needs roughly {estimate.human_bytes(need)} of disk space, and "
                            f"{self.runs_dir} has {estimate.human_bytes(free)} free")
        return problems

    def delete_run(self, name):
        d = self.run_dir(name)
        if not self.is_run(d):
            raise FileNotFoundError(name)
        st = self.run_summary(d)
        if st["state"] in RUNNING_STATES or st["orphans"]:
            raise ValueError("the run is still running; cancel it first")
        shutil.rmtree(d)
        with self.lock:
            self.children.pop(str(d), None)

    def start_run(self, raw):
        raw = self.with_runs_dir(raw, self.runs_dir)
        name = raw["run"].get("name", "")
        d = self.run_dir(name)
        if d.exists():
            raise ValueError(f"a run named {name!r} already exists; choose another name")
        cfg = config.validate(raw, base_dir=self.runs_dir)   # report errors now
        problems = runner.check_setup(cfg) + self.resource_problems(cfg)
        if problems:
            raise ValueError("; ".join(problems))
        self.submit(d, cfg)
        self.schedule()
        return name

    @staticmethod
    def submit(d, cfg):
        d.mkdir(parents=True, mode=0o755)
        runner.private(d)
        saved = d / "submitted.toml"
        saved.write_text(config.dumps(cfg, comment="Submitted from gcf gui"))
        runner.private(saved)
        (d / "queued").write_text("new")

    # -- scans ---------------------------------------------------------------
    def scan_points(self, raw, axes):
        """The points of a scan, each validated: a mistake in one starts none of them."""
        if not isinstance(raw, dict) or not isinstance(axes, list) or not all(isinstance(a, dict) for a in axes):
            raise ValueError("expected a config and a list of parameters")
        points = scan.expand(raw, axes)
        checked = []
        for name, values, point in points:
            if not valid_name(name):
                raise ValueError(f"invalid run name {name!r}: letters, digits, . _ + - only")
            if (self.runs_dir / name).exists():
                raise ValueError(f"a run named {name!r} already exists; choose another scan name")
            cfg = config.validate(self.with_runs_dir(point, self.runs_dir), base_dir=self.runs_dir)
            checked.append((name, values, cfg))
        return checked

    def preview_scan(self, raw, axes):
        points = self.scan_points(raw, axes)
        generate = [cfg["generate"]["events"] for _, _, cfg in points if "generate" in cfg["run"]["stages"]]
        return {"points": len(points), "first": points[0][0], "events": sum(generate)}

    def start_scan(self, raw, axes):
        points = self.scan_points(raw, axes)
        group = raw["run"]["name"]
        path = self.runs_dir / SCANS_DIR / f"{group}.json"
        if path.exists():
            raise ValueError(f"a scan named {group!r} already exists; choose another name")
        problems, need = [], 0
        for name, values, cfg in points:
            for problem in runner.check_setup(cfg) + self.resource_problems(cfg):
                problems.append(f"{name}: {problem}" if len(points) > 1 else problem)
            need += estimate.rough_peak_bytes(cfg)
        free = shutil.disk_usage(self.runs_dir).free
        if need / 2 > free:
            problems.append(f"this scan needs roughly {estimate.human_bytes(need)} of disk space, and "
                            f"{self.runs_dir} has {estimate.human_bytes(free)} free")
        if problems:
            raise ValueError("; ".join(problems[:4]) + (f"; and {len(problems) - 4} more" if len(problems) > 4 else ""))
        path.parent.mkdir(exist_ok=True)
        path.write_text(json.dumps({
            "name": group, "created": time.time(), "seed": points[0][2]["run"]["seed"],
            "axes": [{"path": a["path"], "label": scan.axis_label(a["path"])} for a in axes],
            "points": [{"name": name, "values": values} for name, values, _ in points]}, indent=1))
        for name, _, cfg in points:
            self.submit(self.runs_dir / name, cfg)
        self.schedule()
        return {"name": group, "points": len(points)}

    def read_scan(self, name):
        if not valid_name(name):
            raise ValueError(f"invalid scan name {name!r}")
        try:
            return json.loads((self.runs_dir / SCANS_DIR / f"{name}.json").read_text())
        except (OSError, ValueError):
            raise FileNotFoundError(name) from None

    def list_scans(self):
        out = []
        try:
            files = sorted((self.runs_dir / SCANS_DIR).glob("*.json"), key=lambda p: p.stat().st_mtime, reverse=True)
        except OSError:
            return out
        for path in files:
            try:
                s = json.loads(path.read_text())
            except (OSError, ValueError):
                continue
            states = {}
            for p in s["points"]:
                d = self.runs_dir / p["name"]
                try:
                    state = self.run_summary(d)["state"] if self.is_run(d) else "deleted"
                except Exception:  # noqa: BLE001 - one damaged run must not hide the rest
                    state = "unknown"
                states[state] = states.get(state, 0) + 1
            out.append({"name": s["name"], "points": len(s["points"]), "states": states,
                        "axes": [a["label"] for a in s["axes"]]})
        return out

    def scan_table(self, name):
        """One row per point: its state and, once it's done, its numbers."""
        s = self.read_scan(name)
        rows = []
        models = set()
        for p in s["points"]:
            d = self.runs_dir / p["name"]
            row = {"name": p["name"], "values": p["values"], "state": "deleted", "events": None, "attempts": None,
                   "n_eff": None, "yield": None, "yield_err": None}
            if self.is_run(d):
                run = self.run_summary(d)
                row["state"] = run["state"]
                models.add(self.read_config(d).get("generate", {}).get("model"))
                stage = "fsi" if "fsi" in run["outputs"] else "generate"
                o = run["outputs"].get(stage)
                if run["state"] in ("done", "partial") and o and o["exists"]:
                    try:
                        r = self.results(p["name"], stage)
                        row.update(events=r["events"], attempts=r["attempts"], n_eff=r["n_eff"],
                                   **{"yield": r["yield"], "yield_err": r["yield_err"]})
                    except (RuntimeError, OSError, ValueError) as e:
                        row["error"] = str(e)
            rows.append(row)
        notes = []
        for model in models:
            if model in results.SCALE_NOTES:
                notes.append(results.SCALE_NOTES[model])
        low = [r["n_eff"] for r in rows if r["n_eff"] is not None and r["n_eff"] < results.NEFF_WARN]
        if low:
            notes.append(f"{len(low)} point{'s have' if len(low) > 1 else ' has'} N_eff below {results.NEFF_WARN}, "
                         f"so differences between them may be statistical.")
        if any(a["path"] != "run.seed" for a in s["axes"]):
            notes.append("All points share a seed. For a parameter that only changes the weights, they hold the same "
                         "events, so the differences between them are smaller than the errors shown.")
        return {"name": s["name"], "axes": s["axes"], "points": rows, "notes": notes}

    def delete_scan(self, name):
        s = self.read_scan(name)
        runs = [self.runs_dir / p["name"] for p in s["points"]]
        active = [d.name for d in runs if self.is_run(d)
                  and (self.run_summary(d)["state"] in RUNNING_STATES or self.run_summary(d)["orphans"])]
        if active:
            raise ValueError(f"{len(active)} of its runs are still running; cancel them first")
        for d in runs:
            if self.is_run(d):
                (d / "queued").unlink(missing_ok=True)
                shutil.rmtree(d)
        (self.runs_dir / SCANS_DIR / f"{name}.json").unlink(missing_ok=True)

    # -- scheduling ----------------------------------------------------------
    def launch(self, d):
        """Start a queued run's supervisor. It is on its own from here on: no
        terminal, and this server may stop while it runs."""
        marker = d / "queued"
        claimed = d / ".starting"
        try:
            # Only one server can win this, if two look after one directory.
            os.rename(marker, claimed)
        except OSError:
            return False
        try:
            resume = claimed.read_text().strip() == "resume"
            if not safe_to_run(d):
                raise ValueError(f"not started: others can write to this run's directory or config, and what is in "
                                 f"them would run as you. If you trust them: chmod -R go-w {d}")
            if resume:
                # What the page shows until the new supervisor has written its own status.
                st = runner.read_status(d)
                st.update(state="starting", error=None)
                runner.write_status_file(d, st)
            self.spawn(d, resume)
        except (OSError, ValueError) as e:
            st = runner.read_status(d)
            st.update(state="failed", name=d.name, error=str(e))
            try:
                runner.write_status_file(d, st)
            except OSError:
                pass
            claimed.unlink(missing_ok=True)
            return False
        claimed.unlink(missing_ok=True)
        return True

    def recover_claims(self):
        """A server that died between claiming a queued run and starting it left a marker;
        back in the queue with it."""
        for base in self.known_dirs:
            try:
                entries = list(base.iterdir())
            except OSError:
                continue
            for d in entries:
                claimed = d / ".starting"
                try:
                    if time.time() - claimed.stat().st_mtime > 60 and not runner.is_locked(d):
                        os.rename(claimed, d / "queued")
                except OSError:
                    continue

    def spawn(self, d, resume):
        command = [sys.executable, str(runner.REPO / "gcf"), "run", "-q"]
        command += ["--resume", str(d / "config.toml")] if resume else [str(d / "submitted.toml")]
        with open(d / "gcf.log", "a" if resume else "w") as log:
            proc = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
                                    cwd=str(d.parent), start_new_session=True)
        with self.lock:
            self.children[str(d)] = proc

    def queued_runs(self):
        """Run directories waiting to start, the oldest first."""
        found = []
        for base in self.known_dirs:
            try:
                entries = list(base.iterdir())
            except OSError:
                continue
            for d in entries:
                try:
                    found.append((os.stat(d / "queued").st_mtime, str(d), d))
                except OSError:
                    continue
        return [d for _, _, d in sorted(found)]

    @staticmethod
    def cpus_of_config(cfg):
        """CPUs a config uses at a time: its jobs, or the larger stage of the ePIC chain."""
        eic = cfg.get("eic")
        if eic and "eic" in cfg.get("run", {}).get("stages", ["eic"]):
            return max(eic.get("sim_jobs", 1) * eic.get("sim_threads", 1), eic.get("reco_jobs", 1))
        jobs = cfg.get("run", {}).get("jobs", 1)
        return jobs if isinstance(jobs, int) else 1

    def run_cpus(self, d):
        """CPUs a run may use; a run larger than the limit still gets its turn, alone."""
        return max(1, min(self.cpus_of_config(self.read_config(d)), self.max_cpus))

    def cpus_in_use(self):
        """CPUs reserved by the runs in the directories this server looks after,
        including runs started from a terminal or by an earlier server."""
        used = 0
        for base in self.known_dirs:
            try:
                entries = list(base.iterdir())
            except OSError:
                continue
            for d in entries:
                st = runner.read_status(d)
                if st.get("state") in ("running", "merging") and runner.proc_is(st.get("pid"), st.get("pid_start")):
                    used += self.run_cpus(d)
                elif not st:
                    with self.lock:
                        child = self.children.get(str(d))
                    if child is not None and child.poll() is None:
                        used += self.run_cpus(d)
        return used

    def schedule(self):
        """Start queued runs, in order, as far as CPUs are free."""
        with self.schedule_lock:
            free = self.max_cpus - self.cpus_in_use()
            for d in self.queued_runs():
                try:
                    need = self.run_cpus(d)
                    if need > free:
                        break      # strictly in order, so a large run isn't overtaken forever
                    if self.launch(d):
                        free -= need
                except Exception as e:  # noqa: BLE001 - a problem with one run must not stop the queue
                    print(f"gcf gui: run {d.name}: {e!r}", file=sys.stderr, flush=True)

    def reattach(self):
        """Continue the runs that were interrupted rather than stopped: their status
        says they are running, but their supervisor is gone (it was killed, or the
        machine restarted). Stage processes still running get a new supervisor at
        once, or nothing merges their results; the others wait their turn for CPUs."""
        for base in self.known_dirs:
            try:
                entries = list(base.iterdir())
            except OSError:
                continue
            for d in entries:
                try:             # a problem with one run must not stop the others
                    self.reattach_run(d)
                except Exception as e:  # noqa: BLE001
                    print(f"gcf gui: run {d.name}: {e!r}", file=sys.stderr, flush=True)

    def reattach_run(self, d):
        st = runner.read_status(d)
        if st.get("state") not in ("running", "merging") or runner.proc_is(st.get("pid"), st.get("pid_start")):
            return
        if not (d / "config.toml").exists() or runner.is_locked(d) or (d / "queued").exists() \
                or not safe_to_run(d):
            return
        count_file = d / ".auto_resumes"
        try:
            count = int(count_file.read_text())
        except (OSError, ValueError):
            count = 0
        if count >= MAX_AUTO_RESUMES:
            return             # it shows as stopped, and a person decides
        count_file.write_text(str(count + 1))
        print(f"gcf gui: run {d.name} was interrupted; resuming it", file=sys.stderr, flush=True)
        if self.live_steps(st):
            self.spawn(d, True)
        else:
            (d / "queued").write_text("resume")

    def scheduler(self):
        while True:
            try:
                self.recover_claims()
                self.reattach()
                self.schedule()
            except Exception as e:  # noqa: BLE001 - the page must outlive any one problem
                print(f"gcf gui: scheduler: {e!r}", file=sys.stderr)
            time.sleep(SCHEDULE_SECONDS)

    def machine(self):
        """What the page shows about this machine, next to the runs."""
        try:
            load = os.getloadavg()[0]
        except OSError:
            load = None
        return {"cpus": os.cpu_count() or 1, "max_cpus": self.max_cpus, "in_use": self.cpus_in_use(),
                "queued": len(self.queued_runs()), "load": load,
                "disk_free": shutil.disk_usage(self.runs_dir).free, "runs_dir": str(self.runs_dir)}

    def resume_run(self, name):
        """Queue a run that stopped: finished jobs stay, the rest run again."""
        d = self.run_dir(name)
        if not self.is_run(d):
            raise FileNotFoundError(name)
        st = self.run_summary(d)
        if st["state"] in RUNNING_STATES or st["state"] == "queued" or runner.is_locked(d):
            raise ValueError("the run is still running")
        if not st["resumable"]:
            raise ValueError("this run has nothing to resume")
        (d / ".auto_resumes").unlink(missing_ok=True)
        (d / "queued").write_text("resume")
        self.schedule()

    def merge_finished(self, name):
        d = self.run_dir(name)
        if not self.is_run(d):
            raise FileNotFoundError(name)
        st = self.run_summary(d)
        if st["state"] in RUNNING_STATES or st["orphans"]:
            raise ValueError("the run is still running")
        return runner.merge_finished(d)     # takes the run's lock itself

    def cancel_run(self, name):
        d = self.run_dir(name)
        marker = d / "queued"
        if marker.exists():
            marker.unlink(missing_ok=True)
            st = runner.read_status(d)
            st.update(state="cancelled", name=name)
            runner.write_status_file(d, st)
            return
        try:
            st = json.loads((d / "status.json").read_text())
        except (OSError, json.JSONDecodeError):
            raise ValueError("the run is not running") from None
        if st.get("state") in ("running", "merging") and \
                runner.proc_is(st.get("pid"), st.get("pid_start")):
            os.kill(st["pid"], signal.SIGTERM)      # the supervisor cancels its jobs
            return
        # Supervisor gone but stage processes still ours: stop them directly
        # (each runs in its own session, so kill the whole group).
        orphans = self.live_steps(st)
        if not orphans:
            raise ValueError("the run is not running")
        for s in orphans:
            try:
                os.killpg(s["pid"], signal.SIGTERM)
            except (ProcessLookupError, PermissionError):
                pass

    def output_path(self, name, stage):
        d = self.run_dir(name)
        st = json.loads((d / "status.json").read_text())
        path = st.get("outputs", {}).get(stage)
        if not path or not Path(path).is_file():
            raise ValueError(f"no {stage} output yet")
        return d, path

    def quicklook(self, name, stage, cuts=None):
        from . import quicklook
        if stage == "program":
            raise ValueError("plots are for the event-file format, and this generator writes its own")
        _, path = self.output_path(name, stage)
        stats = {}
        with self.plot_lock:
            out = quicklook.make(path, cuts=cuts, stats=stats)
        return {"plot": out.name, "stats": stats}

    def results(self, name, stage):
        """The summary of an output, computed once per version of the file."""
        d, path = self.output_path(name, stage)
        cache = d / "results.json"
        attempts = None
        if stage == "program" and runner.read_status(d).get("state") == "done":
            # These generators store only the attempts with a positive weight.
            attempts = self.read_config(d).get("program", {}).get("events")
        with self.results_lock:
            return results.cached(path, cache, stage) or \
                results.store(path, cache, stage, results.summarize(path, attempts=attempts))

    def nuclear(self, model, a, z, nn):
        """The table values of a setup, for the form. They come from the binary, so
        they are the ones a run uses; one lookup per setup per server."""
        if model not in schema.MODELS:
            raise ValueError(f"unknown model {model!r}")
        if nn and nn not in schema.NN_INTERACTIONS:
            raise ValueError(f"unknown interaction {nn!r}")
        target = (int(a), int(z))
        key = (model, *target, nn or "")
        with self.lock:
            if key in self.nuclear_cache:
                return self.nuclear_cache[key]
        params = runner.nuclear_params(model, target, nn or None)
        with self.lock:
            self.nuclear_cache[key] = params
        return params

    def run_config(self, name):
        try:
            return tomllib.loads((self.run_dir(name) / "config.toml").read_text())
        except (OSError, tomllib.TOMLDecodeError):
            return {}

    def compare(self, a, b, max_events=300000):
        """tools/validation/compare_samples.py on two runs' outputs; `a` and `b` are
        {run, stage}. Its results stay in <runs>/.compare until an output changes."""
        if not isinstance(a, dict) or not isinstance(b, dict):
            raise ValueError("expected samples a and b, each with a run and a stage")
        if not isinstance(max_events, int) or isinstance(max_events, bool) or max_events < 1000:
            raise ValueError("max_events must be at least 1000")
        if "program" in (a.get("stage"), b.get("stage")):
            raise ValueError("comparing is for the event-file format, and this generator writes its own")
        _, path_a = self.output_path(a.get("run"), a.get("stage"))
        _, path_b = self.output_path(b.get("run"), b.get("stage"))
        if path_a == path_b:
            raise ValueError("choose two different samples")
        out_dir = self.runs_dir / COMPARE_DIR
        out_dir.mkdir(exist_ok=True)
        stem = f"{a['run']}.{a['stage']}--{b['run']}.{b['stage']}--{max_events}"
        js, png = out_dir / f"{stem}.json", out_dir / f"{stem}.png"
        newest = max(Path(path_a).stat().st_mtime, Path(path_b).stat().st_mtime)
        if not (js.is_file() and png.is_file() and js.stat().st_mtime > newest):
            proc = subprocess.run([sys.executable, str(COMPARE_SCRIPT), path_a, path_b, "--weighted",
                                   "--max-events", str(max_events), "--json", str(js), "-o", str(png)],
                                  capture_output=True, text=True, timeout=3600, cwd=str(self.runs_dir))
            if proc.returncode not in (0, 1) or not js.is_file():
                missing = re.search(r"No module named '(\w+)'", proc.stderr)
                raise RuntimeError(f"comparing needs the Python package {missing.group(1)}" if missing
                                   else f"the comparison failed: {proc.stderr.strip()[-600:]}")
        out = json.loads(js.read_text())
        out["plot"] = png.name
        out["a"]["name"] = f"{a['run']} ({'PWIA' if a['stage'] == 'generate' else 'with FSI'})"
        out["b"]["name"] = f"{b['run']} ({'PWIA' if b['stage'] == 'generate' else 'with FSI'})"
        ca, cb = self.run_config(a["run"]), self.run_config(b["run"])
        ma, mb = ca.get("generate", {}).get("model"), cb.get("generate", {}).get("model")
        notes = []
        if ma and mb and ma != mb:
            notes.append(f"These samples are from different models ({ma}, {mb}), whose Σw/attempts are not on a "
                         "common scale (docs/PHYSICS.md), so ignore the pull.")
        if ca.get("run", {}).get("seed") and ca["run"]["seed"] == cb.get("run", {}).get("seed") and a["run"] != b["run"]:
            notes.append("These runs used the same seed, so they may share events. The test assumes independent "
                         "samples and is conservative for shared ones.")
        out["notes"] = notes
        return out

    def parse_config(self, text):
        """A TOML config as data, for the form to load; errors say what's wrong with it."""
        if not isinstance(text, str) or len(text) > 200000:
            raise ValueError("this is not a config file")
        try:
            raw = tomllib.loads(text)
        except tomllib.TOMLDecodeError as e:
            raise ValueError(f"this is not valid TOML: {e}") from None
        raw.setdefault("run", {})
        config.validate(self.with_runs_dir(raw, self.runs_dir), base_dir=self.runs_dir)
        return raw

    def set_runs_dir(self, path, create=False):
        """Look at, and start runs in, another directory. Runs in the earlier one go
        on, and its queue is still looked after."""
        p = Path(path or "").expanduser()
        if not path or not p.is_absolute():
            raise ValueError("give a full path")
        p = p.resolve()
        if not p.exists():
            if not create:
                raise ValueError(f"{p} does not exist")
            p.mkdir(parents=True)
        if not p.is_dir() or not os.access(p, os.W_OK):
            raise ValueError(f"{p} is not a directory you can write to")
        with self.lock:
            self.runs_dir = p
            if p not in self.known_dirs:
                self.known_dirs.append(p)
        self.save_dirs()
        return str(p)

    def preview(self, raw):
        raw = self.with_runs_dir(raw, self.runs_dir)
        cfg = config.validate(raw, base_dir=self.runs_dir)
        return config.dumps(cfg, comment="Run with:  gcf run <this file>")

    def list_files(self, path):
        p = Path(path or self.runs_dir).expanduser().resolve()
        if p.is_file():
            p = p.parent
        if not p.is_dir():
            raise FileNotFoundError(str(p))
        dirs, files = [], []
        for e in sorted(p.iterdir(), key=lambda e: e.name.lower()):
            if e.name.startswith("."):
                continue
            try:
                if e.is_dir():
                    dirs.append(e.name)
                elif e.suffix == ".root":
                    files.append({"name": e.name, "size": e.stat().st_size})
            except OSError:
                continue
        return {"path": str(p), "parent": str(p.parent), "dirs": dirs, "files": files}


def estimate_text(c, cfg):
    """A sentence or two about what a run will cost."""
    parts = []
    if c.get("seconds"):
        jobs = cfg["run"]["jobs"]
        parts.append(f"About {estimate.human_time(c['seconds'])} on {jobs} job{'s' if jobs > 1 else ''} at a time")
    if c.get("bytes"):
        size = f"{estimate.human_bytes(c['bytes'])} of files"
        if c["peak_bytes"] > c["bytes"]:
            size += f" ({estimate.human_bytes(c['peak_bytes'])} while merging)"
        parts.append(size + f"; {estimate.human_bytes(c['free_bytes'])} free")
    text = ", ".join(parts) + "." if parts else ""
    if c.get("pilot_stored") == 0:
        text += (f" A test of {c['pilot_events']} attempts stored no events, so this run may too: check the beam "
                 f"energy against the reaction's threshold.")
    if c.get("n_eff_fraction") is not None:
        stored = c.get("pilot_stored") or c["pilot_events"]
        text += (f" A test with {stored} events had N_eff = {c['n_eff_fraction'] * stored:.0f}, "
                 f"so expect at most N_eff of about {c['n_eff_at_most'] * stored / c['pilot_events']:.3g} for this run; "
                 f"with heavy-tailed weights it grows more slowly than the event count.")
    return text.strip()


def make_handler(app):
    class Handler(BaseHTTPRequestHandler):
        server_version = "gcf-gui"

        def log_message(self, fmt, *args):      # keep the terminal quiet
            pass

        # -- helpers ---------------------------------------------------------
        def authorized(self, query):
            token = self.headers.get("X-GCF-Token") or query.get("token", [""])[0]
            if not token.isascii():
                return False
            return secrets.compare_digest(token, app.token)

        def end_headers(self):
            self.send_header("Referrer-Policy", "no-referrer")
            self.send_header("X-Content-Type-Options", "nosniff")
            super().end_headers()

        def send_json(self, obj, status=HTTPStatus.OK):
            body = json.dumps(obj).encode()
            self.send_response(status)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Cache-Control", "no-store")
            self.end_headers()
            self.wfile.write(body)

        def send_error_json(self, msg, status=HTTPStatus.BAD_REQUEST):
            self.send_json({"error": msg}, status)

        def send_file(self, path, download=False):
            """Send a file of any size, a piece at a time, honouring a simple Range
            request so that an interrupted download of an event file can continue."""
            path = Path(path)
            size = path.stat().st_size
            if path.suffix in (".toml", ".log", ".md"):
                ctype = "text/plain"
            else:
                ctype = mimetypes.guess_type(str(path))[0] or "application/octet-stream"
            if ctype.startswith("text/") or ctype == "application/json":
                ctype += "; charset=utf-8"
            start, end = 0, size - 1
            m = re.fullmatch(r"bytes=(\d+)-(\d*)", self.headers.get("Range", ""))
            partial = bool(m and int(m.group(1)) < size)
            if partial:
                start = int(m.group(1))
                end = min(int(m.group(2)), size - 1) if m.group(2) else size - 1
            self.send_response(HTTPStatus.PARTIAL_CONTENT if partial else HTTPStatus.OK)
            self.send_header("Content-Type", ctype)
            self.send_header("Content-Length", str(max(0, end - start + 1)))
            self.send_header("Accept-Ranges", "bytes")
            self.send_header("Cache-Control", "no-store")
            if partial:
                self.send_header("Content-Range", f"bytes {start}-{end}/{size}")
            if download:
                self.send_header("Content-Disposition", f'attachment; filename="{path.name}"')
            self.end_headers()
            with open(path, "rb") as f:
                f.seek(start)
                left = end - start + 1
                while left > 0:
                    chunk = f.read(min(1024 * 1024, left))
                    if not chunk:
                        break
                    self.wfile.write(chunk)
                    left -= len(chunk)

        def read_json(self):
            n = int(self.headers.get("Content-Length", 0))
            if n > 1_000_000:
                raise ValueError("request too large")
            body = json.loads(self.rfile.read(n) or b"{}")
            if not isinstance(body, dict):
                raise ValueError("the request body must be a JSON object")
            return body

        # -- routes ----------------------------------------------------------
        def do_GET(self):
            url = urllib.parse.urlparse(self.path)
            query = urllib.parse.parse_qs(url.query)
            if not self.authorized(query):
                body = (b"<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\"><title>GCF_with3N</title></head>"
                        b"<body><h1>GCF_with3N</h1><p>Open the full URL printed by <code>gcf gui</code>; "
                        b"it ends in <code>?token=...</code>.</p></body></html>")
                self.send_response(HTTPStatus.FORBIDDEN)
                self.send_header("Content-Type", "text/html")
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                self.wfile.write(body)
                return
            parts = [urllib.parse.unquote(p) for p in url.path.strip("/").split("/") if p]
            try:
                if not parts:
                    # No cookie: browsers send a localhost cookie to every port, so the
                    # token stays in the address and in the page's own requests.
                    return self.send_file(STATIC / "index.html")
                if len(parts) == 2 and parts[0] == "static" and MODULE_RE.fullmatch(parts[1]):
                    return self.send_file(STATIC / parts[1])
                if len(parts) == 3 and parts[:2] == ["api", "compare"] and COMPARE_FILE_RE.fullmatch(parts[2]):
                    return self.send_file(app.runs_dir / COMPARE_DIR / parts[2])
                if parts == ["api", "schema"]:
                    return self.send_json({"schema": schema.to_json(), "examples": cli.EXAMPLES,
                                           "doctor": doctor.checks(app.runs_dir), "version": doctor.built_version(),
                                           "eic_software": config.eic_software_defaults(),
                                           "runs_dir": str(app.runs_dir), "cpus": app.max_cpus,
                                           "binaries": all(os.access(runner.BIN_DIR / b, os.X_OK)
                                                           for b in ("gcf_generate", "gcf_fsi"))})
                if parts == ["api", "runs"]:
                    return self.send_json({"runs": app.list_runs(), "machine": app.machine(), "scans": app.list_scans()})
                if len(parts) == 3 and parts[:2] == ["api", "scans"]:
                    return self.send_json(app.scan_table(parts[2]))
                if len(parts) == 3 and parts[:2] == ["api", "runs"]:
                    return self.send_json(app.run_details(parts[2]))
                if len(parts) == 4 and parts[:2] == ["api", "runs"]:
                    d = app.run_dir(parts[2])
                    f = (d / parts[3]).resolve()
                    if f.parent != d.resolve() or f.suffix not in (".png", ".toml", ".log", ".json", ".root") \
                            or not f.is_file():
                        raise FileNotFoundError(parts[3])
                    return self.send_file(f, download="download" in query or f.suffix == ".root")
                if len(parts) == 5 and parts[:2] == ["api", "runs"] and parts[3] == "logs":
                    logs = app.run_dir(parts[2]) / "logs"
                    f = (logs / parts[4]).resolve()
                    if f.parent != logs.resolve() or f.suffix != ".log" or not f.is_file():
                        raise FileNotFoundError(parts[4])
                    return self.send_file(f, download="download" in query)
                if parts == ["api", "nuclear"]:
                    q = {k: query.get(k, [""])[0] for k in ("model", "A", "Z", "nn")}
                    return self.send_json(app.nuclear(q["model"], q["A"], q["Z"], q["nn"]))
                if len(parts) == 5 and parts[:2] == ["api", "runs"] and parts[3] == "plots":
                    plots = app.run_dir(parts[2]) / "plots"
                    f = (plots / parts[4]).resolve()
                    if f.parent != plots.resolve() or f.suffix not in (".png", ".pdf", ".json", ".txt") \
                            or not f.is_file():
                        raise FileNotFoundError(parts[4])
                    return self.send_file(f, download="download" in query)
                if parts == ["api", "files"]:
                    return self.send_json(app.list_files(query.get("path", [""])[0]))
                if parts == ["docs"]:
                    return self.send_file(STATIC / "docs.html")
                if parts == ["api", "docs"]:
                    name = query.get("name", ["PIPELINE"])[0]
                    if name not in ("GUI", "PIPELINE", "PHYSICS"):
                        raise FileNotFoundError(name)
                    return self.send_file(runner.REPO / "docs" / f"{name}.md")
                if len(parts) == 3 and parts[:2] == ["api", "docs"] and re.fullmatch(r"[a-z0-9_-]+\.png", parts[2]):
                    return self.send_file(runner.REPO / "docs" / "images" / parts[2])
                self.send_error_json("not found", HTTPStatus.NOT_FOUND)
            except FileNotFoundError as e:
                self.send_error_json(f"not found: {e}", HTTPStatus.NOT_FOUND)
            except (ValueError, OSError, runner.RunError, subprocess.TimeoutExpired) as e:
                self.send_error_json(str(e))
            except Exception as e:   # noqa: BLE001 - never drop the connection
                print(f"gcf gui: internal error on {self.path}: {e!r}", file=sys.stderr)
                self.send_error_json("internal error (see the gcf gui terminal)", HTTPStatus.INTERNAL_SERVER_ERROR)

        def do_POST(self):
            url = urllib.parse.urlparse(self.path)
            query = urllib.parse.parse_qs(url.query)
            if not self.authorized(query):
                return self.send_error_json("forbidden", HTTPStatus.FORBIDDEN)
            # Only the page's own requests: no other site, and JSON, which a plain form can't send.
            origin = self.headers.get("Origin")
            if origin and urllib.parse.urlparse(origin).netloc != self.headers.get("Host"):
                return self.send_error_json("forbidden: request from another site", HTTPStatus.FORBIDDEN)
            if self.headers.get("Content-Type", "").split(";")[0].strip() != "application/json":
                return self.send_error_json("expected JSON", HTTPStatus.UNSUPPORTED_MEDIA_TYPE)
            parts = [urllib.parse.unquote(p) for p in url.path.strip("/").split("/") if p]
            # A tab still showing another runs folder must not act on a same-named run here.
            claimed = self.headers.get("X-GCF-Runs-Dir")
            if claimed and parts != ["api", "runs-dir"] and urllib.parse.unquote(claimed) != str(app.runs_dir):
                return self.send_error_json("the runs folder was changed in another tab; try again once this "
                                            "page has caught up", HTTPStatus.CONFLICT)
            try:
                body = self.read_json()
                if parts == ["api", "compare"]:
                    return self.send_json(app.compare(body.get("a"), body.get("b"), body.get("max_events", 300000)))
                if parts == ["api", "preview"]:
                    return self.send_json({"toml": app.preview(body)})
                if parts == ["api", "check"]:
                    return self.send_json(app.check(body))
                if len(parts) == 4 and parts[:2] == ["api", "runs"] and parts[3] == "delete":
                    app.delete_run(parts[2])
                    return self.send_json({"ok": True})
                if parts == ["api", "runs"]:
                    return self.send_json({"name": app.start_run(body)})
                if len(parts) == 4 and parts[:2] == ["api", "runs"] and parts[3] == "cancel":
                    app.cancel_run(parts[2])
                    return self.send_json({"ok": True})
                if len(parts) == 4 and parts[:2] == ["api", "runs"] and parts[3] == "resume":
                    app.resume_run(parts[2])
                    return self.send_json({"ok": True})
                if len(parts) == 4 and parts[:2] == ["api", "runs"] and parts[3] == "merge":
                    return self.send_json({"partial": app.merge_finished(parts[2])})
                if parts == ["api", "scans", "preview"]:
                    return self.send_json(app.preview_scan(body.get("config"), body.get("axes")))
                if parts == ["api", "scans"]:
                    return self.send_json(app.start_scan(body.get("config"), body.get("axes")))
                if len(parts) == 4 and parts[:2] == ["api", "scans"] and parts[3] == "delete":
                    app.delete_scan(parts[2])
                    return self.send_json({"ok": True})
                if parts == ["api", "runs-dir"]:
                    return self.send_json({"runs_dir": app.set_runs_dir(body.get("path"), bool(body.get("create")))})
                if parts == ["api", "parse-config"]:
                    return self.send_json(app.parse_config(body.get("text", "")))
                if len(parts) == 4 and parts[:2] == ["api", "runs"] and parts[3] == "quicklook":
                    return self.send_json(app.quicklook(parts[2], body.get("stage", "fsi"), body.get("cuts")))
                if len(parts) == 4 and parts[:2] == ["api", "runs"] and parts[3] == "results":
                    return self.send_json(app.results(parts[2], body.get("stage", "fsi")))
                self.send_error_json("not found", HTTPStatus.NOT_FOUND)
            except (config.ConfigError, runner.RunError, ValueError, RuntimeError,
                    subprocess.TimeoutExpired) as e:
                self.send_error_json(str(e))
            except FileNotFoundError as e:
                self.send_error_json(f"not found: {e}", HTTPStatus.NOT_FOUND)
            except (OSError, json.JSONDecodeError) as e:
                self.send_error_json(str(e))
            except Exception as e:   # noqa: BLE001
                print(f"gcf gui: internal error on {self.path}: {e!r}", file=sys.stderr)
                self.send_error_json("internal error (see the gcf gui terminal)", HTTPStatus.INTERNAL_SERVER_ERROR)

    return Handler


def banner(app, host, port):
    user = os.environ.get("USER", "you")
    machine = os.uname().nodename
    url = f"http://localhost:{port}/?token={app.token}"
    lines = [f"gcf gui: runs are written to {app.runs_dir}", "",
             "On your laptop, open an SSH tunnel to this machine:",
             f"    ssh -N -L {port}:localhost:{port} {user}@{machine}",
             "(VS Code Remote-SSH forwards the port for you), then browse to",
             f"    {url}", "",
             "The token is the password: don't share this URL, and don't redirect this output to a file others can read."]
    if host not in ("127.0.0.1", "localhost", "::1"):
        lines += ["", f"WARNING: listening on {host}, over plain HTTP: others on the network can read the token. "
                      "Use the SSH tunnel instead unless this network is private."]
    return "\n".join(lines)


def serve(host, port, runs_dir, token=None, max_cpus=None, quiet=False):
    app = App(runs_dir, token=token, max_cpus=max_cpus)
    try:
        app.follows_token_file = token_path().read_text().strip() == token
    except OSError:
        pass
    try:
        server = ThreadingHTTPServer((host, port), make_handler(app))
    except OSError as e:
        raise runner.RunError(f"cannot listen on {host}:{port} ({e}); try --port") from None
    threading.Thread(target=app.scheduler, daemon=True).start()
    record = server_record(port)
    try:
        record.parent.mkdir(parents=True, exist_ok=True)
        record.write_text(json.dumps({"pid": os.getpid(), "pid_start": runner.proc_start(os.getpid()), "port": port,
                                      "host": host, "runs_dir": str(app.runs_dir)}))
    except OSError:
        pass
    if not quiet:
        print(banner(app, host, port) + "\n\nCtrl-C or closing this terminal stops the interface; runs already started "
              "keep going. To keep it up after you log out, use: gcf gui --background", flush=True)
        problems = [c for c in doctor.checks(app.runs_dir) if not c["ok"]]
        if problems:
            print("\n" + doctor.text(problems), flush=True)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()
        record.unlink(missing_ok=True)
