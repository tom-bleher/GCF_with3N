#!/usr/bin/env python3
"""Each of these is a way a run could go wrong without saying so, found by review and
reproduced first: results made with other settings kept, a job that couldn't start
ending as a finished run, and so on."""

import json
import os
import shutil
import signal
import subprocess
import sys
import tempfile
import threading
import time
import tomllib
import unittest
import urllib.error
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import harness  # noqa: E402
from harness import REPO, Gui, small_run  # noqa: E402
from gcf import config, runner  # noqa: E402

GCF = str(REPO / "gcf")


def gcf(*args, env=None, timeout=280):
    return subprocess.run([sys.executable, GCF, *map(str, args)], capture_output=True, text=True,
                          timeout=timeout, env=env)


def raw(gui, path, data=None, headers=None, token=True):
    """(status, response headers, body) of a request with exactly these headers."""
    send = dict(headers or {})
    if token:
        send["X-GCF-Token"] = gui.token
    req = urllib.request.Request(f"http://localhost:{gui.port}/{path}", data=data, headers=send)
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            return r.status, r.headers, r.read()
    except urllib.error.HTTPError as e:
        return e.code, e.headers, e.read()


class Case(unittest.TestCase):
    def setUp(self):
        harness.need_binaries()
        self.tmp = Path(tempfile.mkdtemp(prefix="gcf-review-"))
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)
        self.runs = self.tmp / "runs"

    def config(self, cfg):
        """Validated, with its runs in this test's directory."""
        cfg["run"]["output_dir"] = str(self.runs)
        return config.validate(cfg)

    def write(self, cfg, name="c.toml"):
        path = self.tmp / name
        path.write_text(config.dumps(self.config(cfg)))
        return path

    def failed_run(self, cfg, blocked):
        """Run it with one job unable to write its file, which leaves the other jobs'
        results in place, as an interrupted run does."""
        name = cfg["run"]["name"]
        chunks = self.runs / name / "chunks"
        (chunks / blocked).mkdir(parents=True)
        (chunks / blocked / "in_the_way").write_text("")
        out = gcf("run", self.write(cfg), "-q")
        self.assertNotEqual(out.returncode, 0, out.stdout)
        shutil.rmtree(chunks / blocked)
        return self.runs / name

    def gen_only(self, name, events=4000, chunks=4, seed=11, **generate):
        cfg = small_run(name, events=events, jobs=1, fsi=None, seed=seed, **generate)
        cfg["run"]["chunks"] = chunks
        return cfg

    def entries(self, path, tree="events"):
        import uproot
        with uproot.open(path) as f:
            return f[tree].num_entries

    def done_steps(self, cfg, run_dir):
        plan = runner.plan_from_config(self.config(cfg), run_dir, resume=True)
        return sum(s.state == "done" for steps in plan.jobs for s in steps)


class TestResume(Case):
    def test_finished_jobs_are_kept(self):
        cfg = self.gen_only("keep")
        d = self.failed_run(cfg, "pwia_003.root")
        self.assertEqual(self.done_steps(cfg, d), 3)
        kept = (d / "chunks/pwia_000.root").stat().st_mtime
        out = gcf("run", "--resume", d / "config.toml", "-q")
        self.assertEqual(out.returncode, 0, out.stderr)
        self.assertEqual(self.entries(d / "keep_pwia.root"), 4000)
        self.assertFalse((d / "chunks/pwia_000.root").exists() and (d / "chunks/pwia_000.root").stat().st_mtime != kept)

    def test_an_edited_config_keeps_nothing(self):
        cfg = self.gen_only("edited")
        d = self.failed_run(cfg, "pwia_003.root")
        self.assertEqual(self.done_steps(self.gen_only("edited", sigma_cm=0.2), d), 0)
        self.assertEqual(self.done_steps(self.gen_only("edited", seed=12), d), 0)

    def test_fsi_only_with_another_number_of_jobs_does_not_double_count(self):
        """Other jobs cover other event ranges; keeping the old ones gave 2500 of 2000."""
        import uproot
        out = gcf("run", self.write(small_run("input", events=2000, jobs=1, fsi=None, seed=5)), "-q")
        self.assertEqual(out.returncode, 0, out.stderr)
        source = self.runs / "input/input_pwia.root"
        cfg = {"run": {"name": "only", "stages": ["fsi"], "jobs": 1, "seed": 7, "chunks": 4},
               "fsi": {"model": "hN", "input": str(source)}}
        d = self.failed_run(cfg, "fsi_003.root")
        saved = (d / "config.toml").read_text()
        self.assertIn("chunks = 4", saved)
        (d / "config.toml").write_text(saved.replace("chunks = 4", "chunks = 5"))
        out = gcf("run", "--resume", d / "config.toml", "-q")
        self.assertEqual(out.returncode, 0, out.stderr)
        self.assertEqual(self.entries(d / "only_fsi.root"), 2000)
        with uproot.open(source) as a, uproot.open(d / "only_fsi.root") as b:
            self.assertEqual(a["hAttempts"].values().sum(), b["hAttempts"].values().sum())

    def test_fsi_only_records_its_number_of_jobs(self):
        out = gcf("run", self.write(small_run("input2", events=500, jobs=1, fsi=None, seed=5)), "-q")
        self.assertEqual(out.returncode, 0, out.stderr)
        cfg = {"run": {"name": "only2", "stages": ["fsi"], "jobs": 2, "seed": 7},
               "fsi": {"model": "hN", "input": str(self.runs / "input2/input2_pwia.root")}}
        out = gcf("run", self.write(cfg), "-q")
        self.assertEqual(out.returncode, 0, out.stderr)
        saved = tomllib.loads((self.runs / "only2/config.toml").read_text())
        self.assertEqual(saved["run"]["chunks"], 2)

    def test_force_then_interrupted_then_resumed_uses_nothing_of_the_old_run(self):
        d = self.failed_run(self.gen_only("force", events=8000, seed=1), "pwia_003.root")
        self.assertTrue((d / "chunks/pwia_001.root").is_file())
        path = self.write(self.gen_only("force", events=120000, seed=2), "new.toml")
        before = time.time()
        proc = subprocess.Popen([sys.executable, GCF, "run", str(path), "--force", "-q"],
                                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        record = d / "logs/generate_000.pid"
        deadline = time.time() + 60
        while time.time() < deadline and not (record.is_file() and record.stat().st_mtime >= before):
            time.sleep(0.02)
        os.kill(proc.pid, signal.SIGKILL)
        proc.wait()
        out = gcf("run", "--resume", d / "config.toml", "-q")
        self.assertEqual(out.returncode, 0, out.stderr)
        self.assertEqual(self.entries(d / "force_pwia.root"), 120000)

    def test_a_job_that_cannot_start_does_not_end_as_a_finished_run(self):
        """A worker that died left a run reported as done, with no file."""
        plan = runner.plan_from_config(self.config(self.gen_only("crash", events=400, chunks=2)), self.runs / "crash")
        r = runner.Runner(plan, quiet=True)
        original, first = r.run_step, plan.jobs[0][0]

        def run_step(step):
            if step is first:
                raise RuntimeError("no such thing")
            return original(step)
        r.run_step = run_step
        r.run()
        self.assertEqual(r.state, "failed")
        self.assertIn("no such thing", r.error)
        self.assertFalse((self.runs / "crash/crash_pwia.root").exists())

    def test_a_status_that_cannot_be_written_does_not_stop_the_jobs(self):
        d = self.runs / "nostatus"
        plan = runner.plan_from_config(self.config(small_run("nostatus", events=300, jobs=1, fsi=None)), d)
        (d / "status.json").mkdir(parents=True)
        r = runner.Runner(plan, quiet=True)
        r.run()
        self.assertEqual(r.state, "done", r.error)
        self.assertEqual(self.entries(d / "nostatus_pwia.root"), 300)

    def test_a_merge_of_finished_jobs_survives_a_resume_until_it_is_replaced(self):
        cfg = self.gen_only("part")
        d = self.failed_run(cfg, "pwia_003.root")
        partial = runner.merge_finished(d)
        self.assertTrue(partial)
        plan = runner.plan_from_config(self.config(cfg), d, resume=True)
        self.assertEqual(runner.Runner(plan, quiet=True).status()["partial"], partial)

    def test_a_merge_waits_for_the_supervisor_and_for_running_jobs(self):
        d = self.failed_run(self.gen_only("locked"), "pwia_003.root")
        lock = runner.RunLock(d)
        with self.assertRaises(runner.RunError):
            runner.merge_finished(d)
        lock.file.close()
        base = d / "logs/generate_003"
        record = base.with_suffix(".pid")
        record.unlink()
        proc = subprocess.Popen(["/bin/sh", "-c", runner.STAGE_WRAPPER, str(base), "sleep", "60"], start_new_session=True)
        self.addCleanup(lambda: (os.killpg(proc.pid, signal.SIGKILL), proc.wait()))
        deadline = time.time() + 5
        while time.time() < deadline and not (record.exists() and record.read_text().strip()):
            time.sleep(0.02)
        with self.assertRaises(runner.RunError):
            runner.merge_finished(d)

    def test_status_writers_do_not_get_in_each_others_way(self):
        d = self.tmp / "many"
        d.mkdir()
        errors = []

        def write():
            for i in range(200):
                try:
                    runner.write_status_file(d, {"i": i})
                except OSError as e:
                    errors.append(e)
        threads = [threading.Thread(target=write) for _ in range(8)]
        for t in threads:
            t.start()
        for t in threads:
            t.join()
        self.assertEqual(errors, [])
        self.assertEqual(sorted(p.name for p in d.iterdir()), ["status.json"])


class TestSeeds(Case):
    def seeds(self, base, **kw):
        cfg = self.config(small_run("s", events=4000, jobs=4, seed=base, **kw))
        plan = runner.plan_from_config(cfg, self.runs / "s")
        out = {}
        for steps in plan.jobs:
            for s in steps:
                args = list(map(str, s.args))
                out.setdefault(s.stage, []).append(int(args[args.index("--seed") + 1]))
        return out

    def test_neighbouring_base_seeds_share_no_jobs(self):
        """With seed + i, runs with base seeds 100 and 101 share all but one job, and
        look like independent replicas when they aren't."""
        for stage in ("generate", "fsi"):
            a, b = self.seeds(100)[stage], self.seeds(101)[stage]
            self.assertEqual(set(a) & set(b), set(), stage)
            self.assertEqual(len(set(a)), 4)

    def test_a_seed_makes_a_run_reproducible(self):
        self.assertEqual(self.seeds(100), self.seeds(100))

    def test_seeds_fit_the_stages(self):
        for seeds in self.seeds(7).values():
            for s in seeds:
                self.assertTrue(1 <= s < 2**31 - 1)


class TestScanCommand(Case):
    def scan(self, path, values="0.1,0.2", *extra):
        return gcf("scan", path, "--set", f"generate.sigma_cm={values}", "-q", *extra)

    def test_points_are_kept_only_if_they_are_this_scan(self):
        cfg = small_run("sc", events=400, jobs=1, fsi=None, seed=0)
        path = self.write(cfg)
        out = self.scan(path)
        self.assertEqual(out.returncode, 0, out.stderr)
        point = self.runs / "sc_sigma_cm-0.1/sc_sigma_cm-0.1_pwia.root"
        made = point.stat().st_mtime
        seed = json.loads((self.runs / ".scans/sc.json").read_text())["seed"]
        self.assertTrue(seed)

        out = self.scan(path)                    # again: nothing to do
        self.assertEqual(out.returncode, 0, out.stderr)
        self.assertEqual(point.stat().st_mtime, made)

        out = self.scan(path, "0.1,0.2,0.3")     # one more point: the same seed
        self.assertEqual(out.returncode, 0, out.stderr)
        self.assertEqual(point.stat().st_mtime, made)
        third = tomllib.loads((self.runs / "sc_sigma_cm-0.3/config.toml").read_text())
        self.assertEqual(third["run"]["seed"], seed)

        changed = self.write(small_run("sc", events=800, jobs=1, fsi=None, seed=0), "changed.toml")
        out = self.scan(changed)                 # other settings under the same name: refused
        self.assertNotEqual(out.returncode, 0)
        self.assertIn("other settings", out.stderr)
        self.assertEqual(point.stat().st_mtime, made)


class TestLauncher(unittest.TestCase):
    def launcher(self, wrapper=True, cfg=False):
        tmp = Path(tempfile.mkdtemp(prefix="gcf-launcher-"))
        self.addCleanup(shutil.rmtree, tmp, ignore_errors=True)
        shutil.copy(GCF, tmp / "gcf")
        (tmp / "tools").symlink_to(REPO / "tools")
        (tmp / ".venv/bin").mkdir(parents=True)
        python = tmp / ".venv/bin/python3"
        python.write_text(f'#!/bin/sh\nexec {sys.executable} "$@"\n')
        python.chmod(0o755)
        if cfg:
            (tmp / ".venv/pyvenv.cfg").write_text("home = /usr/bin\n")
        return tmp / "gcf"

    def run_it(self, path):
        env = {k: v for k, v in os.environ.items() if k not in ("VIRTUAL_ENV", "GCF_IN_VENV")}
        return subprocess.run([sys.executable, str(path), "--help"], capture_output=True, text=True, timeout=20, env=env)

    def test_a_directory_that_is_no_venv_is_left_alone(self):
        self.assertEqual(self.run_it(self.launcher()).returncode, 0)

    def test_a_venv_that_does_not_take_cannot_loop(self):
        self.assertEqual(self.run_it(self.launcher(cfg=True)).returncode, 0)


class TestDoctor(unittest.TestCase):
    def test_works_without_git(self):
        tmp = tempfile.mkdtemp(prefix="gcf-doctor-")
        self.addCleanup(shutil.rmtree, tmp, ignore_errors=True)
        for tool in ("bash", "sh", "env", "dirname", "uname", "cat", "grep", "sed", "root-config"):
            found = shutil.which(tool)
            if found:
                os.symlink(found, Path(tmp) / tool)
        code = f"import sys; sys.path.insert(0, {str(REPO / 'tools')!r}); from gcf import doctor; print(len(doctor.checks({tmp!r})))"
        out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=120,
                             env=dict(os.environ, PATH=tmp))
        self.assertEqual(out.returncode, 0, out.stderr)

    def test_a_group_of_your_own_is_nobody_else(self):
        import grp
        import pwd
        tmp = Path(tempfile.mkdtemp(prefix="gcf-doctor-"))
        self.addCleanup(shutil.rmtree, tmp, ignore_errors=True)
        tmp.chmod(0o777)
        self.assertTrue(runner.others_can_write(tmp))
        tmp.chmod(0o755)
        self.assertFalse(runner.others_can_write(tmp))
        tmp.chmod(0o775)
        me = pwd.getpwuid(os.getuid()).pw_name
        group = grp.getgrgid(os.stat(tmp).st_gid)
        self.assertEqual(runner.others_can_write(tmp), group.gr_name != me or any(m != me for m in group.gr_mem))


class TestComparison(unittest.TestCase):
    """A difference in weights' spread alone isn't a difference in shape."""

    def p(self, seed, n1, s1, n2, s2, power=1.0):
        import numpy as np
        sys.path.insert(0, str(REPO / "tools/validation"))
        import compare_samples
        rng = np.random.default_rng(seed)
        x1, x2 = rng.random(n1), rng.random(n2) ** power
        w1, w2 = rng.lognormal(0, s1, n1), rng.lognormal(0, s2, n2)
        return compare_samples.chi2_weighted(x1, w1, x2, w2, np.linspace(0, 1, 26))[2]

    def test_identical_shapes_with_different_weights_are_not_flagged(self):
        ps = [self.p(i, 5000, 0.3, 5000, 2.5) for i in range(12)]
        ps += [self.p(100 + i, 500, 1.0, 20000, 1.0) for i in range(6)]
        self.assertLessEqual(sum(p < 0.01 for p in ps), 1, ps)

    def test_a_clear_change_under_heavy_weights_is_flagged(self):
        """A dominant event left out of a bootstrap replica used to inflate every replica,
        so a change of 50% in a width, chi2/ndf of 40, came out as p = 0.35."""
        import numpy as np
        sys.path.insert(0, str(REPO / "tools/validation"))
        import compare_samples
        bins = np.linspace(0, 1, 31)
        ps = []
        for seed in range(4):
            rng = np.random.default_rng(seed)
            x1, x2 = np.abs(rng.normal(0, 0.24, 30000)), np.abs(rng.normal(0, 0.16, 30000))
            w1, w2 = rng.lognormal(0, 1.8, 30000), rng.lognormal(0, 1.8, 30000)
            ps.append(compare_samples.chi2_weighted(x1, w1, x2, w2, bins)[2])
        self.assertTrue(all(p < 0.01 for p in ps), ps)

    def test_a_different_shape_is_flagged(self):
        ps = [self.p(i, 20000, 1.0, 20000, 1.0, power=1.2) for i in range(4)]
        self.assertTrue(all(p < 0.01 for p in ps), ps)


class TestServer(unittest.TestCase):
    def setUp(self):
        harness.need_binaries()
        self.gui = Gui()
        self.addCleanup(self.gui.close)

    def submit_by_hand(self, name, cfg=None, marker=b"new"):
        d = self.gui.runs_dir / name
        d.mkdir(parents=True)
        cfg = cfg or small_run(name, events=300, jobs=1, fsi=None, seed=3)
        cfg["run"]["output_dir"] = str(self.gui.runs_dir)
        (d / "submitted.toml").write_text(config.dumps(config.validate(cfg)))
        if marker is not None:
            (d / "queued").write_bytes(marker)
        return d

    def wait_state(self, d, want, timeout=60):
        deadline = time.time() + timeout
        state = None
        while time.time() < deadline:
            state = runner.read_status(d).get("state")
            if state == want:
                return
            time.sleep(0.2)
        self.fail(f"{d.name}: {state}, not {want}")

    def test_a_new_token_takes_effect_at_once(self):
        old = self.gui.token
        out = gcf("gui", "--new-token", env=self.gui.env, timeout=30)
        self.assertEqual(out.returncode, 0, out.stderr)
        new = out.stdout.split("token=")[1].split()[0]
        self.assertNotEqual(new, old)
        self.assertEqual(self.gui.request("api/runs", token=new)[0], 200)
        self.assertEqual(self.gui.request("api/runs", token=old)[0], 403)

    def test_the_token_is_not_in_a_cookie(self):
        """Browsers send a cookie for localhost to all its ports, so to any other
        service there, and the page opened from one could use it."""
        status, headers, _ = raw(self.gui, f"?token={self.gui.token}", token=False)
        self.assertEqual(status, 200)
        self.assertIsNone(headers.get("Set-Cookie"))
        status, _, _ = raw(self.gui, "api/runs", token=False,
                           headers={"Cookie": f"gcf_token_{self.gui.port}={self.gui.token}"})
        self.assertEqual(status, 403)

    def test_only_the_page_itself_can_make_changes(self):
        path = "api/runs/nothing/cancel"
        as_json = {"Content-Type": "application/json"}
        self.assertEqual(raw(self.gui, path, b"{}", {**as_json, "Origin": "http://localhost:1"})[0], 403)
        self.assertEqual(raw(self.gui, path, b"{}", {"Content-Type": "text/plain"})[0], 415)
        self.assertNotIn(raw(self.gui, path, b"{}", {**as_json, "Origin": f"http://localhost:{self.gui.port}"})[0], (403, 415))

    def test_a_tab_showing_another_folder_cannot_act_here(self):
        path = "api/runs/nothing/cancel"
        as_json = {"Content-Type": "application/json"}
        self.assertEqual(raw(self.gui, path, b"{}", {**as_json, "X-GCF-Runs-Dir": "/somewhere/else"})[0], 409)
        here = urllib.request.quote(str(self.gui.runs_dir), safe="")
        self.assertNotEqual(raw(self.gui, path, b"{}", {**as_json, "X-GCF-Runs-Dir": here})[0], 409)

    def test_a_plot_older_than_its_file_is_not_shown(self):
        self.gui.api("api/runs", small_run("plot", events=300, jobs=1, fsi=None, seed=2))
        d = self.gui.wait("plot")
        self.assertEqual(d["state"], "done", d.get("error"))
        root = self.gui.runs_dir / "plot/plot_pwia.root"
        png = root.with_suffix(".png")
        png.write_bytes(b"x")
        os.utime(png, (time.time() - 100, time.time() - 100))
        self.assertFalse(self.gui.api("api/runs/plot")["outputs"]["generate"]["plot"])
        os.utime(png)
        self.assertTrue(self.gui.api("api/runs/plot")["outputs"]["generate"]["plot"])

    def test_a_directory_others_can_write_to_is_not_run(self):
        d = self.submit_by_hand("unsafe", marker=None)
        d.chmod(0o777)
        (d / "queued").write_text("new")
        self.wait_state(d, "failed")
        self.assertIn("not started", runner.read_status(d)["error"])
        self.assertFalse((d / "unsafe_pwia.root").exists())

    def test_one_damaged_run_hides_nothing_and_stops_nothing(self):
        self.gui.api("api/runs", small_run("fine", events=300, jobs=1, fsi=None, seed=2))
        self.gui.wait("fine")
        for name, files in {"binary": {"config.toml": b"\xff\xfe\x00", "status.json": b"{}"},
                            "shapes": {"status.json": b'{"state": 7, "stages": 5, "jobs": 3, "outputs": 1}'}}.items():
            d = self.gui.runs_dir / name
            d.mkdir()
            for f, content in files.items():
                (d / f).write_bytes(content)
        listed = self.gui.api("api/runs")["runs"]
        self.assertEqual({"fine", "binary", "shapes"}, {r["name"] for r in listed})
        for name in ("binary", "shapes"):
            status, out = self.gui.request(f"api/runs/{name}")
            self.assertEqual(status, 200, out)

        first = self.submit_by_hand("first", marker=b"\xff\xfe\x00")
        time.sleep(0.05)
        second = self.submit_by_hand("second")
        self.wait_state(second, "done")
        self.assertTrue(first.is_dir())

    def test_a_claim_left_by_a_dead_server_goes_back_in_the_queue(self):
        d = self.submit_by_hand("claimed", marker=None)
        (d / ".starting").write_text("new")
        os.utime(d / ".starting", (time.time() - 300, time.time() - 300))
        self.wait_state(d, "done")

    def test_all_folders_used_are_still_looked_after_after_a_restart(self):
        other = Path(self.gui.tmp) / "other"
        other.mkdir()
        self.gui.api("api/runs-dir", {"path": str(other)})
        self.gui.api("api/runs-dir", {"path": str(self.gui.runs_dir)})
        self.gui.stop()
        self.gui.start()
        d = other / "there"
        d.mkdir()
        cfg = small_run("there", events=300, jobs=1, fsi=None, seed=3)
        cfg["run"]["output_dir"] = str(other)
        (d / "submitted.toml").write_text(config.dumps(config.validate(cfg)))
        (d / "queued").write_text("new")
        self.wait_state(d, "done")


if __name__ == "__main__":
    unittest.main()
