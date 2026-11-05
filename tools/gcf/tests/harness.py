"""Helpers for the front end's tests: a real `gcf gui` on a free port, its API,
and, when selenium and Firefox are there, a headless browser on its page.

The tests run the real stage binaries. Build them first (./install.sh), or point
GCF_BIN_DIR at them.
"""

import json
import os
import re
import shutil
import socket
import subprocess
import sys
import tempfile
import time
import unittest
import urllib.error
import urllib.request
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "tools"))

from gcf import runner  # noqa: E402


def need_binaries():
    for name in ("gcf_generate", "gcf_fsi"):
        if not os.access(runner.BIN_DIR / name, os.X_OK):
            raise unittest.SkipTest(f"{name} is not built in {runner.BIN_DIR}")


def free_port():
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


class Gui:
    """`gcf gui` as a subprocess; stops it, and removes its directory, on close."""

    def __init__(self, runs_dir=None, args=()):
        self.tmp = tempfile.mkdtemp(prefix="gcf-gui-test-")
        self.runs_dir = Path(runs_dir or Path(self.tmp) / "runs")
        self.port = free_port()
        self.log = Path(self.tmp) / "server.log"
        self.proc = None
        self.start(args)

    def start(self, args=()):
        out = open(self.log, "w")
        # The token and the record of the server stay out of the user's own settings.
        self.env = dict(os.environ, XDG_CONFIG_HOME=str(Path(self.tmp) / "config"))
        self.proc = subprocess.Popen(
            [sys.executable, str(REPO / "gcf"), "gui", "--port", str(self.port),
             "--runs-dir", str(self.runs_dir), *args],
            stdout=out, stderr=subprocess.STDOUT, start_new_session=True, env=self.env)
        out.close()
        deadline = time.time() + 20
        while time.time() < deadline:
            m = re.search(r"token=(\S+)", self.log.read_text())
            if m:
                self.token = m.group(1)
                break
            if self.proc.poll() is not None:
                raise RuntimeError(f"gcf gui exited: {self.log.read_text()}")
            time.sleep(0.05)
        else:
            raise RuntimeError("gcf gui printed no URL")
        while time.time() < deadline:
            try:
                socket.create_connection(("127.0.0.1", self.port), timeout=1).close()
                return
            except OSError:
                time.sleep(0.05)
        raise RuntimeError("gcf gui does not listen")

    def stop(self):
        if self.proc and self.proc.poll() is None:
            self.proc.terminate()
            self.proc.wait(timeout=10)

    def close(self):
        self.stop()
        shutil.rmtree(self.tmp, ignore_errors=True)

    @property
    def url(self):
        return f"http://localhost:{self.port}/?token={self.token}"

    def request(self, path, body=None, token=None):
        """(status, parsed JSON or raw bytes) of an API call."""
        req = urllib.request.Request(
            f"http://localhost:{self.port}/{path}",
            data=None if body is None else json.dumps(body).encode(),
            headers={"X-GCF-Token": token or self.token, "Content-Type": "application/json"})
        try:
            with urllib.request.urlopen(req, timeout=300) as r:
                status, raw = r.status, r.read()
        except urllib.error.HTTPError as e:
            status, raw = e.code, e.read()
        try:
            return status, json.loads(raw)
        except ValueError:
            return status, raw

    def api(self, path, body=None):
        status, out = self.request(path, body)
        if status != 200:
            raise AssertionError(f"{path}: HTTP {status}: {out}")
        return out

    def wait(self, name, timeout=280):
        """Details of a run once it stopped running."""
        deadline = time.time() + timeout
        while time.time() < deadline:
            d = self.api(f"api/runs/{name}")
            if d["state"] not in ("running", "merging", "starting", "queued"):
                return d
            time.sleep(0.5)
        raise AssertionError(f"run {name} still {d['state']} after {timeout} s")


def small_run(name, model="2N", events=2000, jobs=2, fsi="hN", seed=1234, **generate):
    """Config of a small run; `fsi=None` makes it generate-only."""
    gen = {"model": model, "beam": 6.0, "events": events, "target": [12, 6], **generate}
    cfg = {"run": {"name": name, "stages": ["generate", "fsi"] if fsi else ["generate"],
                   "jobs": jobs, "seed": seed}, "generate": gen}
    if fsi:
        cfg["fsi"] = {"model": fsi}
    return cfg


def browser(width=1600, height=1050):
    """A headless Firefox, or SkipTest if selenium or Firefox is missing."""
    try:
        from selenium import webdriver
        from selenium.webdriver.firefox.service import Service
    except ImportError:
        raise unittest.SkipTest("selenium is not installed") from None
    geckodriver = shutil.which("geckodriver")
    if not geckodriver or not shutil.which("firefox"):
        raise unittest.SkipTest("geckodriver or firefox is missing")
    if Path("/snap/bin/firefox").exists():
        # The snap can't see the system's /tmp, where the profile would go.
        tmp = Path.home() / "snap/firefox/common/tmp"
        tmp.mkdir(parents=True, exist_ok=True)
        os.environ["TMPDIR"] = str(tmp)
    options = webdriver.FirefoxOptions()
    options.add_argument("-headless")
    options.add_argument(f"--width={width}")
    options.add_argument(f"--height={height}")
    return webdriver.Firefox(options=options, service=Service(geckodriver))
