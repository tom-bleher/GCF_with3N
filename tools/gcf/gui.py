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
NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.+-]{0,99}$")


def pid_alive(pid):
    try:
        os.kill(pid, 0)
        return True
    except (ProcessLookupError, TypeError):
        return False
    except PermissionError:
        return True


class App:
    def __init__(self, runs_dir):
        self.runs_dir = Path(runs_dir).resolve()
        self.runs_dir.mkdir(parents=True, exist_ok=True)
        self.token = secrets.token_urlsafe(18)
        self.lock = threading.Lock()
        self.children = {}            # run name -> Popen of `gcf run`

    # -- runs ----------------------------------------------------------------
    def run_dir(self, name):
        if not NAME_RE.match(name or ""):
            raise ValueError(f"invalid run name {name!r}")
        return self.runs_dir / name

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
        if child is not None:
            child.poll()               # reap it once it has exited
        if state in ("running", "merging") and not pid_alive(st.get("pid")):
            state = "stopped"          # the run's process is gone (killed, machine reboot)
        if state == "starting" and child is not None and child.poll() is not None:
            state = "failed"
            st["error"] = (d / "gcf.log").read_text(errors="replace")[-2000:] \
                if (d / "gcf.log").exists() else "gcf run exited"
        outputs = {}
        for stage, path in st.get("outputs", {}).items():
            p = Path(path)
            outputs[stage] = {"path": path, "exists": p.is_file(),
                              "size": p.stat().st_size if p.is_file() else 0,
                              "plot": p.with_suffix(".png").is_file()}
        return {"name": d.name, "state": state, "error": st.get("error"),
                "started": st.get("started"), "finished": st.get("finished"),
                "stages": st.get("stages", {}), "outputs": outputs,
                "mtime": d.stat().st_mtime}

    def list_runs(self):
        runs = [self.run_summary(d) for d in self.runs_dir.iterdir()
                if d.is_dir() and ((d / "status.json").exists() or (d / "submitted.toml").exists())]
        return sorted(runs, key=lambda r: r["mtime"], reverse=True)

    def run_details(self, name):
        d = self.run_dir(name)
        if not d.is_dir():
            raise FileNotFoundError(name)
        out = self.run_summary(d)
        cfg = d / "config.toml"
        out["config"] = cfg.read_text() if cfg.exists() else None
        logs = []
        for log in sorted((d / "logs").glob("*.log")) if (d / "logs").is_dir() else []:
            lines = log.read_text(errors="replace").splitlines()
            lines = [l for l in lines if not l.startswith("Progress:")]
            logs.append({"name": log.name, "tail": "\n".join(lines[-15:])})
        out["logs"] = logs
        return out

    def start_run(self, raw):
        raw = dict(raw)
        run = dict(raw.get("run", {}))
        name = run.get("name", "")
        d = self.run_dir(name)
        if d.exists():
            raise ValueError(f"a run named {name!r} already exists; choose another name")
        run["output_dir"] = str(self.runs_dir)
        raw["run"] = run
        cfg = config.validate(raw, base_dir=self.runs_dir)   # report errors now
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
        st = json.loads((d / "status.json").read_text())
        pid = st.get("pid")
        if st.get("state") not in ("running", "merging") or not pid_alive(pid):
            raise ValueError("the run is not running")
        os.kill(pid, signal.SIGTERM)

    def quicklook(self, name, stage):
        from . import quicklook
        d = self.run_dir(name)
        st = json.loads((d / "status.json").read_text())
        path = st.get("outputs", {}).get(stage)
        if not path or not Path(path).is_file():
            raise ValueError(f"no {stage} output yet")
        out = quicklook.make(path)
        return out.name

    def preview(self, raw):
        raw = dict(raw)
        raw["run"] = dict(raw.get("run", {}), output_dir=str(self.runs_dir))
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
            return secrets.compare_digest(token, app.token)

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
            data = Path(path).read_bytes()
            if path.suffix in (".toml", ".log"):
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
            self.end_headers()
            self.wfile.write(data)

        def read_json(self):
            n = int(self.headers.get("Content-Length", 0))
            if n > 1_000_000:
                raise ValueError("request too large")
            return json.loads(self.rfile.read(n) or b"{}")

        # -- routes ----------------------------------------------------------
        def do_GET(self):
            url = urllib.parse.urlparse(self.path)
            query = urllib.parse.parse_qs(url.query)
            if not self.authorized(query):
                body = (b"<h1>gcf gui</h1><p>Open the full URL printed by <code>gcf gui</code>; "
                        b"it ends in <code>?token=...</code>.</p>")
                self.send_response(HTTPStatus.FORBIDDEN)
                self.send_header("Content-Type", "text/html")
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                self.wfile.write(body)
                return
            parts = [urllib.parse.unquote(p) for p in url.path.strip("/").split("/") if p]
            try:
                if not parts:
                    return self.send_file(STATIC / "index.html")
                if parts == ["api", "schema"]:
                    return self.send_json({"schema": schema.to_json(), "examples": cli.EXAMPLES,
                                           "runs_dir": str(app.runs_dir),
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
                self.send_error_json("not found", HTTPStatus.NOT_FOUND)
            except FileNotFoundError as e:
                self.send_error_json(f"not found: {e}", HTTPStatus.NOT_FOUND)
            except (ValueError, OSError) as e:
                self.send_error_json(str(e))

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
                if parts == ["api", "runs"]:
                    return self.send_json({"name": app.start_run(body)})
                if len(parts) == 4 and parts[:2] == ["api", "runs"] and parts[3] == "cancel":
                    app.cancel_run(parts[2])
                    return self.send_json({"ok": True})
                if len(parts) == 4 and parts[:2] == ["api", "runs"] and parts[3] == "quicklook":
                    return self.send_json({"plot": app.quicklook(parts[2], body.get("stage", "fsi"))})
                self.send_error_json("not found", HTTPStatus.NOT_FOUND)
            except (config.ConfigError, runner.RunError, ValueError, RuntimeError) as e:
                self.send_error_json(str(e))
            except FileNotFoundError as e:
                self.send_error_json(f"not found: {e}", HTTPStatus.NOT_FOUND)
            except (OSError, json.JSONDecodeError) as e:
                self.send_error_json(str(e))

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
          f"\nCtrl-C stops the interface; runs already started keep going.", flush=True)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()
