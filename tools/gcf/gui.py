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
import subprocess
import sys
import threading
import time
import urllib.parse
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

from . import cli, config, runner, schema

STATIC = Path(__file__).resolve().parent / "static"
NAME_RE = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.+-]{0,99}\Z")
PLACEHOLDER_STATE = "failed"


def valid_name(name):
    return bool(name) and NAME_RE.fullmatch(name) is not None


class App:
    def __init__(self, runs_dir):
        self.runs_dir = Path(runs_dir).resolve()
        self.runs_dir.mkdir(parents=True, exist_ok=True)
        self.token = secrets.token_urlsafe(18)
        self.lock = threading.Lock()
        self.children = {}            # run name -> Popen of `gcf run`
        self.plot_lock = threading.Lock()   # matplotlib is not thread-safe

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
        st = {}
        status = d / "status.json"
        if status.is_file():
            try:
                st = json.loads(status.read_text())
            except (OSError, json.JSONDecodeError):
                st = {}
        state = st.get("state", "starting" if (d / "submitted.toml").exists() else "unknown")
        with self.lock:
            child = self.children.get(d.name)
            if child is not None and child.poll() is not None:
                del self.children[d.name]   # reap and forget finished supervisors
                child = None
        # The supervisor is only "running" if its pid is still the recorded
        # process (pid + start time); a reused pid after a crash or reboot is
        # not ours and must never be signalled.
        supervisor_alive = runner.proc_is(st.get("pid"), st.get("pid_start"))
        orphans = []
        if state in ("running", "merging") and not supervisor_alive:
            state = "stopped"          # the supervisor is gone (killed, machine reboot)
            orphans = self.live_steps(st)
        if state == "starting" and child is None:
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
            outputs[stage] = {"path": path, "exists": p.is_file(),
                              "size": p.stat().st_size if p.is_file() else 0,
                              "plot": p.with_suffix(".png").is_file()}
        return {"name": d.name, "state": state, "error": st.get("error"),
                "started": st.get("started"), "finished": st.get("finished"),
                "stages": st.get("stages", {}), "outputs": outputs,
                "orphans": len(orphans),
                "setup": self.setup_summary(d), "mtime": d.stat().st_mtime}

    def setup_summary(self, d):
        """One line saying what a run is, from its config: '3N · 12C · 6 GeV · 40k · hN'."""
        cfg = d / "config.toml"
        if not cfg.exists():
            cfg = d / "submitted.toml"
        try:
            c = tomllib.loads(cfg.read_text())
        except (OSError, tomllib.TOMLDecodeError):
            return ""
        parts = []
        g, f, r = c.get("generate"), c.get("fsi"), c.get("run", {})
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
        for d in self.runs_dir.iterdir():
            try:
                if self.is_run(d):
                    runs.append(self.run_summary(d))
            except OSError:
                continue           # vanished between listing and reading
        return sorted(runs, key=lambda r: r["mtime"], reverse=True)

    def run_details(self, name):
        d = self.run_dir(name)
        if not d.is_dir():
            raise FileNotFoundError(name)
        out = self.run_summary(d)
        try:
            out["jobs"] = json.loads((d / "status.json").read_text()).get("jobs", [])
        except (OSError, json.JSONDecodeError):
            out["jobs"] = []
        cfg = d / "config.toml"
        if not cfg.exists():
            cfg = d / "submitted.toml"
        out["config"] = cfg.read_text() if cfg.exists() else None
        out["config_data"] = None
        if cfg.exists():
            try:
                out["config_data"] = tomllib.loads(out["config"])
            except tomllib.TOMLDecodeError:
                pass
        logs = []
        order = {"generate": 0, "fsi": 1, "merge": 2}
        files = sorted((d / "logs").glob("*.log"), key=lambda p: (order.get(p.name.split("_")[0], 9), p.name)) \
            if (d / "logs").is_dir() else []
        for log in files:
            lines = log.read_text(errors="replace").splitlines()
            # drop progress lines and GENIE's start-up banner / messenger chatter
            lines = [l for l in lines if l.strip() and not l.startswith("Progress:")
                     and not l.lstrip().startswith("*") and "NOTICE Messenger" not in l]
            logs.append({"name": log.name, "tail": "\n".join(lines[-40:])})
        gcf_log = d / "gcf.log"
        if gcf_log.exists() and gcf_log.stat().st_size:
            logs.insert(0, {"name": "gcf.log", "tail": gcf_log.read_text(errors="replace")[-3000:]})
        out["logs"] = logs
        return out

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
        resolved TOML and a list of problems."""
        raw = self.with_runs_dir(raw, self.runs_dir)
        cfg = config.validate(raw, base_dir=self.runs_dir)
        problems = []
        name = cfg["run"]["name"]
        if not valid_name(name):
            problems.append(f"invalid run name {name!r}: letters, digits, . _ + - only")
        elif (self.runs_dir / name).exists():
            problems.append(f"a run named {name!r} already exists; choose another name")
        problems += runner.check_setup(cfg)
        return {"toml": config.dumps(cfg, comment="Run with:  gcf run <this file>"),
                "problems": problems, "ok": not problems}

    def delete_run(self, name):
        d = self.run_dir(name)
        if not self.is_run(d):
            raise FileNotFoundError(name)
        st = self.run_summary(d)
        if st["state"] in ("running", "merging", "starting") or st["orphans"]:
            raise ValueError("the run is still running; cancel it first")
        shutil.rmtree(d)
        with self.lock:
            self.children.pop(name, None)

    def start_run(self, raw):
        raw = self.with_runs_dir(raw, self.runs_dir)
        name = raw["run"].get("name", "")
        d = self.run_dir(name)
        if d.exists():
            raise ValueError(f"a run named {name!r} already exists; choose another name")
        cfg = config.validate(raw, base_dir=self.runs_dir)   # report errors now
        problems = runner.check_setup(cfg)
        if problems:
            raise ValueError("; ".join(problems))
        d.mkdir(parents=True)
        submitted = d / "submitted.toml"
        submitted.write_text(config.dumps(cfg, comment="Submitted from gcf gui"))
        log = open(d / "gcf.log", "w")
        proc = subprocess.Popen([sys.executable, str(runner.REPO / "gcf"), "run", "-q",
                                 str(submitted)], stdout=log, stderr=subprocess.STDOUT,
                                cwd=str(self.runs_dir), start_new_session=True)
        log.close()
        with self.lock:
            self.children[name] = proc
        return name

    def cancel_run(self, name):
        d = self.run_dir(name)
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

    def quicklook(self, name, stage):
        from . import quicklook
        d = self.run_dir(name)
        st = json.loads((d / "status.json").read_text())
        path = st.get("outputs", {}).get(stage)
        if not path or not Path(path).is_file():
            raise ValueError(f"no {stage} output yet")
        with self.plot_lock:
            out = quicklook.make(path)
        return out.name

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


def make_handler(app):
    class Handler(BaseHTTPRequestHandler):
        server_version = "gcf-gui"

        def log_message(self, fmt, *args):      # keep the terminal quiet
            pass

        # -- helpers ---------------------------------------------------------
        def authorized(self, query):
            token = self.headers.get("X-GCF-Token") or query.get("token", [""])[0]
            if not token:
                cookie = self.headers.get("Cookie", "")
                for part in cookie.split(";"):
                    k, _, v = part.strip().partition("=")
                    if k == "gcf_token":
                        token = v
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

        def send_file(self, path, download=False, cookie=None):
            data = Path(path).read_bytes()
            if path.suffix in (".toml", ".log", ".md"):
                ctype = "text/plain"
            else:
                ctype = mimetypes.guess_type(str(path))[0] or "application/octet-stream"
            if ctype.startswith("text/") or ctype == "application/json":
                ctype += "; charset=utf-8"
            self.send_response(HTTPStatus.OK)
            self.send_header("Content-Type", ctype)
            self.send_header("Content-Length", str(len(data)))
            self.send_header("Cache-Control", "no-store")
            if download:
                self.send_header("Content-Disposition", f'attachment; filename="{path.name}"')
            if cookie:
                self.send_header("Set-Cookie", cookie)
            self.end_headers()
            self.wfile.write(data)

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
                    # The page keeps using the token from its URL; the cookie
                    # lets a bookmark without the token keep working.
                    return self.send_file(STATIC / "index.html",
                                          cookie=f"gcf_token={app.token}; HttpOnly; SameSite=Strict; Path=/")
                if parts == ["api", "schema"]:
                    return self.send_json({"schema": schema.to_json(), "examples": cli.EXAMPLES,
                                           "runs_dir": str(app.runs_dir), "cpus": os.cpu_count() or 1,
                                           "binaries": all(os.access(runner.BIN_DIR / b, os.X_OK)
                                                           for b in ("gcf_generate", "gcf_fsi"))})
                if parts == ["api", "runs"]:
                    return self.send_json({"runs": app.list_runs()})
                if len(parts) == 3 and parts[:2] == ["api", "runs"]:
                    return self.send_json(app.run_details(parts[2]))
                if len(parts) == 4 and parts[:2] == ["api", "runs"]:
                    d = app.run_dir(parts[2])
                    f = (d / parts[3]).resolve()
                    if f.parent != d.resolve() or f.suffix not in (".png", ".toml", ".log", ".json") \
                            or not f.is_file():
                        raise FileNotFoundError(parts[3])
                    return self.send_file(f, download="download" in query)
                if parts == ["api", "files"]:
                    return self.send_json(app.list_files(query.get("path", [""])[0]))
                if parts == ["docs"]:
                    return self.send_file(STATIC / "docs.html")
                if parts == ["api", "docs"]:
                    return self.send_file(runner.REPO / "docs" / "PIPELINE.md")
                self.send_error_json("not found", HTTPStatus.NOT_FOUND)
            except FileNotFoundError as e:
                self.send_error_json(f"not found: {e}", HTTPStatus.NOT_FOUND)
            except (ValueError, OSError) as e:
                self.send_error_json(str(e))
            except Exception as e:   # noqa: BLE001 - never drop the connection
                print(f"gcf gui: internal error on {self.path}: {e!r}", file=sys.stderr)
                self.send_error_json("internal error (see the gcf gui terminal)", HTTPStatus.INTERNAL_SERVER_ERROR)

        def do_POST(self):
            url = urllib.parse.urlparse(self.path)
            query = urllib.parse.parse_qs(url.query)
            if not self.authorized(query):
                return self.send_error_json("forbidden", HTTPStatus.FORBIDDEN)
            parts = [urllib.parse.unquote(p) for p in url.path.strip("/").split("/") if p]
            try:
                body = self.read_json()
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
                if len(parts) == 4 and parts[:2] == ["api", "runs"] and parts[3] == "quicklook":
                    return self.send_json({"plot": app.quicklook(parts[2], body.get("stage", "fsi"))})
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


def serve(host, port, runs_dir):
    app = App(runs_dir)
    try:
        server = ThreadingHTTPServer((host, port), make_handler(app))
    except OSError as e:
        raise runner.RunError(f"cannot listen on {host}:{port} ({e}); try --port") from None
    user = os.environ.get("USER", "you")
    machine = os.uname().nodename
    url = f"http://localhost:{port}/?token={app.token}"
    print(f"gcf gui: runs are written to {app.runs_dir}\n"
          f"\n"
          f"On your laptop, open an SSH tunnel to this machine:\n"
          f"    ssh -N -L {port}:localhost:{port} {user}@{machine}\n"
          f"(VS Code Remote-SSH forwards the port for you), then browse to\n"
          f"    {url}\n"
          f"\nThe token is the password: don't redirect this output to a file others can read.\n"
          f"Ctrl-C stops the interface; runs already started keep going.", flush=True)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()
