#!/usr/bin/env python3
"""Tests of the web interface, against a real `gcf gui` and the real stages.

    python3 -m unittest discover -s tools/gcf/tests -v

Browser tests are skipped without selenium and Firefox.
"""

import json
import os
import re
import shutil
import signal
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import harness  # noqa: E402
from harness import REPO, Gui, small_run  # noqa: E402


class GuiCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        harness.need_binaries()
        cls.gui = Gui()
        cls.addClassCleanup(cls.gui.close)


class TestResults(GuiCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.gui.api("api/runs", small_run("res", events=3000, jobs=2))
        cls.details = cls.gui.wait("res")
        cls.pwia = cls.details["outputs"]["generate"]["path"]
        cls.fsi = cls.details["outputs"]["fsi"]["path"]

    def test_summary_matches_the_file(self):
        import numpy as np
        import uproot
        s = self.gui.api("api/runs/res/results", {"stage": "generate"})
        with uproot.open(self.pwia) as f:
            w = f["events"]["weight"].array(library="np")
            attempts = float(f["hAttempts"].values()[0])
        self.assertEqual(s["events"], 3000)
        self.assertEqual(s["attempts"], attempts)
        self.assertAlmostEqual(s["yield"], w.sum() / attempts)
        self.assertAlmostEqual(s["yield_err"], float(np.sqrt((w**2).sum())) / attempts)
        self.assertAlmostEqual(s["n_eff"], float(w.sum()**2 / (w**2).sum()))
        self.assertAlmostEqual(s["max_w_fraction"], float(w.max() / w.sum()))
        self.assertEqual(len(s["provenance"]["gcf_generate"]["seeds"]), 2)
        self.assertEqual(s["model"], "2N")

    def test_summary_is_cached_and_follows_the_file(self):
        self.gui.api("api/runs/res/results", {"stage": "fsi"})
        d = self.gui.api("api/runs/res")
        self.assertIn("fsi", d["results"])
        self.assertAlmostEqual(d["results"]["fsi"]["fsi_ran_fraction"], 1.0)

    def test_weights_do_not_change_in_fsi(self):
        a = self.gui.api("api/runs/res/results", {"stage": "generate"})
        b = self.gui.api("api/runs/res/results", {"stage": "fsi"})
        self.assertAlmostEqual(a["yield"], b["yield"])
        self.assertAlmostEqual(a["n_eff"], b["n_eff"])

    def test_selection_matches_a_direct_calculation(self):
        import numpy as np
        import uproot
        r = self.gui.api("api/runs/res/quicklook", {"stage": "generate", "cuts": {"q2": [2, 5], "xb": [1.0, None]}})
        with uproot.open(self.pwia) as f:
            a = f["events"].arrays(["weight", "beam_E", "e_px", "e_py", "e_pz", "e_E"], library="np")
        qz, nu = a["beam_E"] - a["e_pz"], a["beam_E"] - a["e_E"]
        q2 = a["e_px"]**2 + a["e_py"]**2 + qz**2 - nu**2
        keep = (q2 >= 2) & (q2 <= 5) & (q2 / (2 * 0.93892 * nu) >= 1.0)
        w = a["weight"]
        self.assertEqual(r["stats"]["events"], int(keep.sum()))
        self.assertAlmostEqual(r["stats"]["fraction"], float(w[keep].sum() / w.sum()))
        self.assertEqual(r["plot"], "res_pwia_cuts.png")
        self.assertEqual(self.gui.request(f"api/runs/res/{r['plot']}")[0], 200)

    def test_selection_with_no_events_says_so(self):
        status, out = self.gui.request("api/runs/res/quicklook", {"stage": "generate", "cuts": {"q2": [1000, None]}})
        self.assertEqual(status, 400)
        self.assertIn("no event passes", out["error"])

    def test_unknown_cut_is_refused(self):
        status, out = self.gui.request("api/runs/res/quicklook", {"stage": "generate", "cuts": {"bogus": [1, 2]}})
        self.assertEqual(status, 400)
        self.assertIn("unknown cut", out["error"])

    def test_chunks_do_not_change_the_histograms(self):
        """Filling chunk by chunk gives the same bins as filling in one go."""
        from gcf import quicklook
        whole, parts = {}, {}
        quicklook.make(self.fsi, Path(self.gui.tmp) / "a.png", stats=whole)
        quicklook.make(self.fsi, Path(self.gui.tmp) / "b.png", stats=parts, step_size=250)
        self.assertEqual(whole["events_all"], 3000)
        self.assertEqual(parts["events_all"], 3000)
        for panel in ("x$_B$", "protons per event", "neutrons per event, after FSI",
                      "partner |p|, before FSI [GeV/c]"):
            for label, counts in whole["histograms"][panel].items():
                for x, y in zip(counts, parts["histograms"][panel][label]):
                    self.assertAlmostEqual(x, y)

    def test_panels_follow_the_model(self):
        """A pair has one partner, a triplet two, and only a triplet gets the angle map."""
        from gcf import quicklook
        two, three = {}, {}
        self.gui.api("api/runs", small_run("three", model="3N", events=500, jobs=1))
        d = self.gui.wait("three")
        quicklook.make(self.pwia, Path(self.gui.tmp) / "2.png", stats=two)
        quicklook.make(d["outputs"]["fsi"]["path"], Path(self.gui.tmp) / "3.png", stats=three)
        angle_map = "$\\theta_{12}$ (struck, partner 1) [deg]"
        self.assertNotIn(angle_map, two["panels"])
        self.assertIn(angle_map, three["panels"])
        partners = "partner |p|, before FSI [GeV/c]"
        self.assertEqual(sum(1 for c in two["histograms"][partners].values() if sum(c) > 0), 1)
        self.assertEqual(sum(1 for c in three["histograms"][partners].values() if sum(c) > 0), 2)
        self.assertIn("neutrons per event, after FSI", three["panels"])


class TestCompare(GuiCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        for name, seed, fsi, extra in (("cmp_a", 11, "hN", {}), ("cmp_b", 22, "hN", {}), ("cmp_ha", 11, "hA", {}),
                                       ("cmp_sig", 33, "hN", {"sigma_cm": 0.05})):
            cls.gui.api("api/runs", small_run(name, events=30000, jobs=4, seed=seed, fsi=fsi, **extra))
        for name in ("cmp_a", "cmp_b", "cmp_ha", "cmp_sig"):
            cls.gui.wait(name)

    def compare(self, a, b, stage="fsi"):
        return self.gui.api("api/compare", {"a": {"run": a, "stage": stage}, "b": {"run": b, "stage": stage}})

    def test_same_settings_are_not_flagged(self):
        r = self.compare("cmp_a", "cmp_b")
        self.assertTrue(r["ok"], r["observables"])
        self.assertGreater(len(r["observables"]), 9)
        self.assertEqual(r["notes"], [])
        self.assertEqual(self.gui.request(f"api/compare/{r['plot']}")[0], 200)

    def test_another_cascade_model_is_flagged(self):
        """hN and hA on the same seed: the same generated events, other final states."""
        r = self.compare("cmp_a", "cmp_ha")
        by_name = {o["name"]: o for o in r["observables"]}
        self.assertLess(by_name["protons after FSI"]["p"], 1e-6)
        self.assertEqual(by_name["Q2 [GeV^2]"]["chi2"], 0)
        self.assertAlmostEqual(r["pull"], 0)
        self.assertTrue(any("same seed" in n for n in r["notes"]))

    def test_sigma_cm_is_seen_in_the_observables_that_respond_to_it(self):
        r = self.compare("cmp_a", "cmp_sig", stage="generate")
        by_name = {o["name"]: o for o in r["observables"]}
        self.assertLessEqual(by_name["residual nucleus |p| [GeV/c]"]["p_weighted"], 1e-3)
        self.assertFalse(r["ok"])

    def test_files_without_fsi_get_the_generated_observables(self):
        r = self.compare("cmp_a", "cmp_b", stage="generate")
        names = [o["name"] for o in r["observables"]]
        self.assertIn("partner |p|, generated [GeV/c]", names)
        self.assertNotIn("protons after FSI", names)

    def test_a_sample_with_itself_is_refused(self):
        status, out = self.gui.request("api/compare", {"a": {"run": "cmp_a", "stage": "fsi"}, "b": {"run": "cmp_a", "stage": "fsi"}})
        self.assertEqual(status, 400)
        self.assertIn("two different samples", out["error"])

    def test_weighted_test_sees_a_change_in_the_weights_only(self):
        """The same events, with other weights: invisible to the unweighted test."""
        import numpy as np
        sys.path.insert(0, str(REPO / "tools/validation"))
        import compare_samples as C
        rng = np.random.default_rng(1)
        x = rng.normal(0, 1, 40000)
        w_a, w_b = np.ones_like(x), np.exp(0.3 * x)
        bins = np.linspace(-4, 4, 21)
        self.assertEqual(C.chi2_unweighted(x, x, bins)[0], 0)
        self.assertLess(C.chi2_weighted(x, w_a, x, w_b, bins)[2], 0.01)


class TestPhysicsFlags(unittest.TestCase):
    """What the nuclear-input options of gcf_generate do to the events."""

    @classmethod
    def setUpClass(cls):
        harness.need_binaries()
        cls.tmp = tempfile.mkdtemp(prefix="gcf-flags-")
        cls.addClassCleanup(shutil.rmtree, cls.tmp, ignore_errors=True)
        cls.env = harness.runner.environment()
        cls.binary = harness.runner.binary("gcf_generate")
        cls.count = 0

    def generate(self, *args, model="2N", events=3000, seed=4242):
        import uproot
        type(self).count += 1
        out = str(Path(self.tmp) / f"e{self.count}.root")
        run = subprocess.run([self.binary, "--model", model, "--beam", "6", "--events", str(events), "--seed", str(seed),
                              "-A", "12", "-Z", "6", "--output", out, *map(str, args)],
                             capture_output=True, text=True, env=self.env)
        self.assertEqual(run.returncode, 0, run.stdout + run.stderr)
        with uproot.open(out) as f:
            a = f["events"].arrays(["weight", "e_px", "e_pz", "res_px", "res_py", "res_pz", "res_E", "nuc_pdg"], library="np")
            a["attempts"] = float(f["hAttempts"].values()[0])
            a["provenance"] = str(f["gcf_generate"].member("fTitle"))
        return a

    def params(self, *args):
        return harness.runner.nuclear_params(*args)

    def assertSameEvents(self, a, b):
        import numpy as np
        for k in ("e_px", "e_pz", "res_E"):
            self.assertTrue(np.array_equal(a[k], b[k]), k)

    @staticmethod
    def residual_mass(a):
        import numpy as np
        return np.sqrt(a["res_E"]**2 - a["res_px"]**2 - a["res_py"]**2 - a["res_pz"]**2)

    def test_contacts_equal_to_the_table_change_nothing(self):
        import numpy as np
        c = self.params("2N", [12, 6])["contacts"]
        a, b = self.generate(), self.generate("--cpp0", c["cpp0"][0], "--cnn0", c["cnn0"][0], "--cpn0", c["cpn0"][0], "--cpn1", c["cpn1"][0])
        self.assertSameEvents(a, b)
        self.assertTrue(np.array_equal(a["weight"], b["weight"]))

    def test_scaling_all_contacts_scales_the_weights(self):
        import numpy as np
        c = self.params("2N", [12, 6])["contacts"]
        a = self.generate()
        b = self.generate("--cpp0", 2 * c["cpp0"][0], "--cnn0", 2 * c["cnn0"][0], "--cpn0", 2 * c["cpn0"][0], "--cpn1", 2 * c["cpn1"][0])
        self.assertSameEvents(a, b)
        self.assertTrue(np.allclose(b["weight"], 2 * a["weight"]))

    def test_one_contact_only_moves_the_weights_of_its_pairs(self):
        import numpy as np
        c = self.params("2N", [12, 6])["contacts"]
        a, b = self.generate(), self.generate("--cpn1", 2 * c["cpn1"][0])
        self.assertSameEvents(a, b)
        ratio = b["weight"] / a["weight"]
        pn = np.array([sorted(p[:2]) == [2112, 2212] for p in a["nuc_pdg"]])
        self.assertTrue(np.allclose(ratio[~pn], 1))
        self.assertTrue((ratio[pn] > 1).all() and (ratio[pn] <= 2 + 1e-9).all())

    def test_zero_contacts_remove_their_pairs(self):
        a = self.generate("--cpp0", 0, "--cnn0", 0)
        self.assertTrue(all(sorted(p[:2]) == [2112, 2212] for p in a["nuc_pdg"]))

    def test_negative_or_all_zero_contacts_are_refused(self):
        for args in (["--cpp0", "-1"], ["--cpp0", "0", "--cnn0", "0", "--cpn0", "0", "--cpn1", "0"]):
            run = subprocess.run([self.binary, "--model", "2N", "--beam", "6", "--check", *args],
                                 capture_output=True, text=True, env=self.env)
            self.assertEqual(run.returncode, 1, args)

    def test_contacts_do_not_apply_to_the_deuteron(self):
        run = subprocess.run([self.binary, "--model", "2N", "--beam", "6", "-A", "2", "-Z", "1", "--check", "--cpn1", "50"],
                             capture_output=True, text=True, env=self.env)
        self.assertEqual(run.returncode, 1)
        self.assertIn("deuteron", run.stderr)

    def test_estar_shifts_the_residual_mass(self):
        import numpy as np
        a, b = self.generate(), self.generate("--estar", 0.02)
        # The events differ, so the mix of the three A-2 residuals does slightly.
        self.assertAlmostEqual(float(np.mean(self.residual_mass(b)) - np.mean(self.residual_mass(a))), 0.02, delta=2e-4)

    def test_3n_estar(self):
        """Without --estar a 3N run is unchanged; with it the residual is heavier by E*."""
        import numpy as np
        a, b = self.generate(model="3N", events=500), self.generate("--estar", 0, model="3N", events=500)
        self.assertSameEvents(a, b)
        self.assertTrue(np.array_equal(a["weight"], b["weight"]))
        c = self.generate("--estar", 0.02, model="3N", events=500)
        self.assertAlmostEqual(float(np.mean(self.residual_mass(c)) - np.mean(self.residual_mass(a))), 0.02, places=5)

    def test_estar_width_draws_per_event(self):
        import numpy as np
        b = self.generate("--estar", 0.02, "--sigma-e", 0.005)
        # Each residual is the A-2 ground state of its pair type plus its own E*.
        ground = self.params("2N", [12, 6])["mRes"]
        kinds = {(2212, 2212): "pp", (2112, 2112): "nn"}
        states = np.array([ground[kinds.get(tuple(int(x) for x in p[:2]), "pn")] for p in b["nuc_pdg"]])
        excitation = self.residual_mass(b) - states
        self.assertTrue((excitation > 0).all())
        self.assertAlmostEqual(float(np.mean(excitation)), 0.02, delta=0.002)
        self.assertAlmostEqual(float(np.std(excitation)), 0.005, delta=0.0015)

    def test_randomize_is_the_same_for_every_job(self):
        """One draw per run: jobs with other event seeds get the same nuclear inputs."""
        a = self.generate("--randomize", 77, seed=1, events=50)
        b = self.generate("--randomize", 77, seed=2, events=50)
        c = self.generate("--randomize", 78, seed=1, events=50)
        inputs = lambda e: [l for l in e["provenance"].splitlines() if l.startswith("nuclear inputs")][0]
        self.assertEqual(inputs(a), inputs(b))
        self.assertNotEqual(inputs(a), inputs(c))

    def test_an_explicit_value_overrides_the_randomized_one(self):
        a = self.generate("--randomize", 77, "--sigma-cm", 0.111, "--cpn1", 50, events=50)
        line = [l for l in a["provenance"].splitlines() if l.startswith("nuclear inputs")][0]
        self.assertIn("sigma_cm=0.111", line)
        self.assertIn("cpn1=50", line)

    def test_light_cone_gives_other_events(self):
        a, b = self.generate(events=1000), self.generate("--light-cone", events=1000)
        self.assertEqual(len(b["weight"]), 1000)
        self.assertNotEqual(a["attempts"], b["attempts"])

    def test_light_cone_excludes_options_it_would_ignore(self):
        for extra in (["--prel-is"], ["--onset", "local"]):
            run = subprocess.run([self.binary, "--model", "2N", "--beam", "6", "--check", "--light-cone", *extra],
                                 capture_output=True, text=True, env=self.env)
            self.assertEqual(run.returncode, 1, extra)

    def test_form_factors_change_the_weights_only(self):
        import numpy as np
        a, b = self.generate(), self.generate("--ff", "ye")
        self.assertSameEvents(a, b)
        self.assertFalse(np.array_equal(a["weight"], b["weight"]))

    def test_params_of_the_setups_the_form_offers(self):
        p = self.params("2N", [12, 6])
        self.assertEqual(p["contacts"]["cpn1"][0], 95.256)
        self.assertEqual(p["sigma_cm"], [0.15, 0.02])
        self.assertIn("k-space", p["contact_source"])
        self.assertEqual(self.params("2N", [12, 6], "N2LO")["contacts"]["cpn1"][0], 63.858)
        self.assertIn("EG2", self.params("2N", [56, 26])["contact_source"])
        self.assertEqual(self.params("3N", [12, 6])["sigma_cm"][0], 0.15)
        with self.assertRaises(harness.runner.RunError):
            self.params("2N", [12, 6], "N3LO")


class TestChannelScales(unittest.TestCase):
    def test_config_round_trip_and_command_line(self):
        from gcf import config, runner
        cfg = config.validate({"run": {"stages": ["generate", "fsi"]}, "generate": {"model": "2N"},
                               "fsi": {"scale": {"pion.absorption": 0.5, "pion.cex": 0}}})
        self.assertEqual(cfg["fsi"]["scale"], {"pion.absorption": 0.5})
        self.assertEqual(cfg["fsi"]["off"], ["pion.cex"])
        import tomllib
        again = tomllib.loads(config.dumps(cfg))
        self.assertEqual(again["fsi"]["scale"], {"pion.absorption": 0.5})
        old, runner.binary = runner.binary, lambda n: n
        try:
            args = runner.fsi_args(cfg["fsi"], "i.root", "o.root", 1)
        finally:
            runner.binary = old
        self.assertIn("pion.absorption=0.5", args)
        self.assertIn("pion.cex=0", args)

    def test_bad_scales_are_refused(self):
        from gcf import config
        for scale in ({"pion.absorption": -1}, {"bogus": 1}, {"nucleon.compound": 0.5}, {"nucleon.cex": 2}):
            with self.assertRaises(config.ConfigError):
                config.validate({"run": {"stages": ["generate", "fsi"]}, "generate": {"model": "2N"},
                                 "fsi": {"model": "hN", "scale": scale}})

    def test_a_scale_changes_the_cascade(self):
        """Pion absorption at half its rate leaves more charged pions."""
        harness.need_binaries()
        gui = Gui()
        self.addCleanup(gui.close)
        for name, fsi in (("nominal", {}), ("half", {"scale": {"pion.absorption": 0.1}})):
            cfg = small_run(name, events=40000, jobs=8)
            cfg["fsi"].update(fsi)
            gui.api("api/runs", cfg)
        a, b = gui.wait("nominal"), gui.wait("half")
        self.assertEqual((a["state"], b["state"]), ("done", "done"), a.get("error") or b.get("error"))
        import awkward as ak
        import numpy as np
        import uproot
        def pions(d):
            with uproot.open(d["outputs"]["fsi"]["path"]) as f:
                return float(np.mean(np.asarray(ak.sum(abs(f["events"]["part_pdg"].array()) == 211, axis=1))))
        self.assertGreater(pions(b), pions(a))


class TestResilience(GuiCase):
    """A run goes on whatever happens to what started it, and loses little if it can't."""

    def status(self, name):
        return json.loads((self.gui.runs_dir / name / "status.json").read_text())

    def wait_for(self, condition, what, timeout=60):
        deadline = time.time() + timeout
        while time.time() < deadline:
            if condition():
                return
            time.sleep(0.2)
        self.fail(f"timed out waiting for {what}")

    def stage_pids(self, name):
        return [s["pid"] for job in self.status(name)["jobs"] for s in job if s.get("pid")]

    def test_jobs_outlive_a_killed_supervisor_and_get_a_new_one(self):
        """They used to die within seconds, of SIGPIPE, with the supervisor."""
        self.gui.api("api/runs", small_run("orphan", model="3N", events=60000, jobs=2, seed=9))
        path = self.gui.runs_dir / "orphan/status.json"
        self.wait_for(lambda: path.exists() and self.status("orphan")["stages"]["fsi"]["done"] > 0, "the cascade to start")
        first = self.status("orphan")
        pids = self.stage_pids("orphan")
        self.assertTrue(pids)
        os.kill(first["pid"], signal.SIGKILL)
        time.sleep(1.0)
        self.assertTrue(any(harness.runner.proc_start(p) is not None for p in pids), "the jobs died with the supervisor")
        self.wait_for(lambda: self.status("orphan")["pid"] != first["pid"], "a new supervisor", 20)
        d = self.gui.wait("orphan")
        self.assertEqual(d["state"], "done", d.get("error"))
        self.assertTrue(d["outputs"]["fsi"]["exists"])
        self.assertEqual(self.gui.api("api/runs/orphan/results", {"stage": "fsi"})["events"], 60000)

    def test_interrupted_runs_continue_when_the_server_comes_back(self):
        """As after a reboot: supervisor and jobs are gone, and no server was running."""
        import uproot
        gui = Gui()
        self.addCleanup(gui.close)
        cfg = small_run("reboot", model="2N", events=40000, jobs=2, seed=12)
        cfg["run"]["chunks"] = 8
        gui.api("api/runs", cfg)
        path = gui.runs_dir / "reboot/status.json"
        self.wait_for(lambda: path.exists() and json.loads(path.read_text())["stages"]["fsi"]["finished"] >= 2, "some finished jobs")
        gui.stop()
        st = json.loads(path.read_text())
        os.kill(st["pid"], signal.SIGKILL)
        for record in (gui.runs_dir / "reboot/logs").glob("*.pid"):
            try:
                os.killpg(int(record.read_text()), signal.SIGKILL)
            except OSError:
                pass
        time.sleep(0.5)
        gui.start()
        d = gui.wait("reboot")
        self.assertEqual(d["state"], "done", d.get("error"))
        with uproot.open(d["outputs"]["fsi"]["path"]) as f:
            self.assertEqual(f["events"].num_entries, 40000)
        self.assertIn("was interrupted; resuming", gui.log.read_text())

    def test_a_run_that_keeps_dying_waits_for_a_person(self):
        from gcf import gui as gui_module
        self.gui.api("api/runs", small_run("dying", events=500, jobs=1, seed=44, fsi=None))
        self.gui.wait("dying")
        run_dir = self.gui.runs_dir / "dying"
        st = self.status("dying")
        st.update(state="running", pid=999999, pid_start=1)
        (run_dir / ".auto_resumes").write_text(str(gui_module.MAX_AUTO_RESUMES))
        (run_dir / "status.json").write_text(json.dumps(st))
        time.sleep(gui_module.SCHEDULE_SECONDS * 2.5)
        self.assertEqual(self.gui.api("api/runs/dying")["state"], "stopped")
        self.assertTrue(self.gui.api("api/runs/dying")["resumable"])

    def test_a_job_the_last_status_missed_is_not_started_twice(self):
        """The supervisor writes its status once a second. A job started just before
        it died isn't in there, and a second copy would write to the same file."""
        cfg = small_run("missed", model="3N", events=60000, jobs=2, seed=19)
        self.gui.api("api/runs", cfg)
        path = self.gui.runs_dir / "missed/status.json"
        self.wait_for(lambda: path.exists() and self.status("missed")["stages"]["fsi"]["done"] > 0, "the cascade to start")
        st = self.status("missed")
        os.kill(st["pid"], signal.SIGKILL)
        records = {p.name: int(p.read_text()) for p in (self.gui.runs_dir / "missed/logs").glob("fsi_*.pid")}
        self.assertTrue(records)
        # Forget them, as a status written a moment earlier would have.
        for job in st["jobs"]:
            for step in job:
                step["pid"] = step["pid_start"] = None
        path.write_text(json.dumps(st))
        self.wait_for(lambda: self.status("missed")["pid"] != st["pid"], "a new supervisor", 20)
        d = self.gui.wait("missed")
        self.assertEqual(d["state"], "done", d.get("error"))
        after = {p.name: int(p.read_text()) for p in (self.gui.runs_dir / "missed/logs").glob("fsi_*.pid")}
        self.assertEqual(records, after, "a job was started a second time")

    def test_runs_outlive_the_server_and_it_finds_them_again(self):
        gui = Gui()
        self.addCleanup(gui.close)
        gui.api("api/runs", small_run("alone", model="3N", events=40000, jobs=2, seed=4))
        time.sleep(1.5)
        token = gui.token
        os.killpg(gui.proc.pid, signal.SIGHUP)     # as when the terminal closes
        gui.proc.wait(timeout=10)
        time.sleep(1)
        st = json.loads((gui.runs_dir / "alone/status.json").read_text())
        self.assertTrue(harness.runner.proc_is(st["pid"], st["pid_start"]))
        gui.start()
        self.assertEqual(gui.token, token, "the address stays valid across a restart")
        d = gui.wait("alone")
        self.assertEqual(d["state"], "done", d.get("error"))

    def test_resume_reruns_only_what_did_not_finish_and_gives_the_same_events(self):
        import numpy as np
        import uproot
        whole = small_run("whole", model="2N", events=40000, jobs=2, seed=77)
        whole["run"]["chunks"] = 8
        self.gui.api("api/runs", whole)
        reference = self.gui.wait("whole")
        self.assertEqual(reference["state"], "done", reference.get("error"))
        cut = small_run("cut", model="2N", events=40000, jobs=2, seed=77)
        cut["run"]["chunks"] = 8
        self.gui.api("api/runs", cut)
        self.wait_for(lambda: (self.gui.runs_dir / "cut/status.json").exists()
                      and self.status("cut")["stages"]["fsi"]["finished"] >= 2, "some finished jobs")
        self.gui.api("api/runs/cut/cancel", {})
        d = self.gui.wait("cut")
        self.assertEqual(d["state"], "cancelled")
        self.assertTrue(d["resumable"])
        finished = [i for i, job in enumerate(self.status("cut")["jobs"]) if all(s["state"] == "done" for s in job)]
        self.assertTrue(finished)
        path = self.gui.runs_dir / f"cut/chunks/fsi_{finished[0]:03d}.root"
        kept = path.stat().st_mtime_ns
        self.gui.api("api/runs/cut/resume", {})
        d = self.gui.wait("cut")
        self.assertEqual(d["state"], "done", d.get("error"))
        with open(self.gui.runs_dir / "cut/config.toml") as f:
            self.assertIn("chunks = 8", f.read())
        # Same seeds and jobs: the same events as the run that was never interrupted.
        for stage in ("generate", "fsi"):
            a = uproot.open(reference["outputs"][stage]["path"])["events"]["weight"].array(library="np")
            b = uproot.open(d["outputs"][stage]["path"])["events"]["weight"].array(library="np")
            self.assertTrue(np.array_equal(a, b), stage)
        self.assertEqual(self.gui.api("api/runs/cut/results", {"stage": "fsi"})["events"], 40000)
        self.assertGreater(kept, 0)

    def test_a_failed_job_does_not_take_the_others_with_it(self):
        cfg = small_run("onefail", model="2N", events=8000, jobs=1, seed=5, fsi=None)
        cfg["run"]["chunks"] = 4
        self.gui.api("api/runs", cfg)
        d = self.gui.wait("onefail")
        self.assertEqual(d["state"], "done", d.get("error"))
        # A run where one job can't write its output.
        cfg = small_run("onefail2", model="2N", events=80000, jobs=1, seed=5, fsi=None)
        cfg["run"]["chunks"] = 4
        run_dir = self.gui.runs_dir / "onefail2"
        (run_dir / "chunks").mkdir(parents=True)
        (run_dir / "chunks/pwia_001.root").mkdir()   # a directory where a file must go,
        (run_dir / "chunks/pwia_001.root/keep").write_text("")   # and not empty, or ROOT replaces it
        (run_dir / "submitted.toml").write_text("")
        import tomllib  # noqa: F401
        from gcf import config
        full = config.validate(dict(cfg, run=dict(cfg["run"], output_dir=str(self.gui.runs_dir))))
        (run_dir / "config.toml").write_text(config.dumps(full))
        out = subprocess.run([sys.executable, str(REPO / "gcf"), "run", "-q", "--resume", str(run_dir / "config.toml")],
                             capture_output=True, text=True, env=harness.runner.environment())
        self.assertEqual(out.returncode, 1, out.stderr)
        st = self.status("onefail2")
        states = [job[0]["state"] for job in st["jobs"]]
        self.assertEqual(states.count("failed"), 1, states)
        self.assertEqual(states.count("done"), 3, states)
        self.assertIn("1 of 4 jobs failed; 3 finished and are kept", st["error"])
        # Fixed: only that job runs again.
        shutil.rmtree(run_dir / "chunks/pwia_001.root")
        self.gui.api("api/runs/onefail2/resume", {})
        d = self.gui.wait("onefail2")
        self.assertEqual(d["state"], "done", d.get("error"))

    def test_finished_jobs_can_be_merged_with_the_right_yield(self):
        cfg = small_run("part", model="2N", events=40000, jobs=2, seed=31)
        cfg["run"]["chunks"] = 8
        self.gui.api("api/runs", cfg)
        self.wait_for(lambda: (self.gui.runs_dir / "part/status.json").exists()
                      and self.status("part")["stages"]["fsi"]["finished"] >= 2, "some finished jobs")
        self.gui.api("api/runs/part/cancel", {})
        d = self.gui.wait("part")
        self.assertTrue(d["mergeable"])
        out = self.gui.api("api/runs/part/merge", {})
        d = self.gui.wait("part")
        self.assertEqual(d["state"], "partial")
        self.assertTrue(d["outputs"]["fsi"]["exists"])
        self.assertEqual(out["partial"]["fsi"]["total"], 40000)
        s = self.gui.api("api/runs/part/results", {"stage": "fsi"})
        self.assertEqual(s["events"], out["partial"]["fsi"]["events"])
        # Its attempts are those of the merged jobs, so the yield isn't diluted.
        self.assertAlmostEqual(s["yield"], 9.4, delta=2.5)
        self.assertTrue(d["resumable"])
        self.gui.api("api/runs/part/resume", {})
        d = self.gui.wait("part")
        self.assertEqual(d["state"], "done", d.get("error"))
        self.assertEqual(self.gui.api("api/runs/part/results", {"stage": "fsi"})["events"], 40000)

    def test_a_run_has_one_supervisor(self):
        cfg_path = None
        self.gui.api("api/runs", small_run("single", model="3N", events=30000, jobs=1, seed=8))
        self.wait_for(lambda: (self.gui.runs_dir / "single/config.toml").exists(), "the run's config")
        cfg_path = self.gui.runs_dir / "single/config.toml"
        out = subprocess.run([sys.executable, str(REPO / "gcf"), "run", "-q", "--resume", str(cfg_path)],
                             capture_output=True, text=True, env=harness.runner.environment())
        self.assertEqual(out.returncode, 2)
        self.assertIn("being run by another process", out.stderr)
        self.gui.wait("single")

    def test_long_runs_get_checkpoints(self):
        from gcf import runner
        self.assertEqual(runner.checkpoint_jobs({"jobs": 8}, 100000), 8)
        self.assertEqual(runner.checkpoint_jobs({"jobs": 8}, 100000000), 100)
        self.assertEqual(runner.checkpoint_jobs({"jobs": 8, "chunks": 20}, 100000000), 20)


class TestQueue(unittest.TestCase):
    def test_runs_wait_for_cpus_in_order_and_can_be_removed(self):
        harness.need_binaries()
        gui = Gui(args=["--max-cpus", "2"])
        self.addCleanup(gui.close)
        gui.api("api/runs", small_run("q1", model="3N", events=30000, jobs=2, seed=1, fsi="hN"))
        gui.api("api/runs", small_run("q2", events=500, jobs=1, seed=2, fsi=None))
        gui.api("api/runs", small_run("q3", events=500, jobs=1, seed=3, fsi=None))
        states = {r["name"]: r["state"] for r in gui.api("api/runs")["runs"]}
        self.assertEqual(states["q2"], "queued")
        self.assertEqual(states["q3"], "queued")
        machine = gui.api("api/runs")["machine"]
        self.assertEqual((machine["max_cpus"], machine["queued"]), (2, 2))
        gui.api("api/runs/q3/cancel", {})
        self.assertEqual(gui.api("api/runs/q3")["state"], "cancelled")
        self.assertEqual(gui.wait("q1")["state"], "done")
        self.assertEqual(gui.wait("q2")["state"], "done")
        self.assertEqual(gui.api("api/runs/q3")["state"], "cancelled")

    def test_too_many_jobs_for_the_limit_are_refused(self):
        harness.need_binaries()
        gui = Gui(args=["--max-cpus", "2"])
        self.addCleanup(gui.close)
        status, out = gui.request("api/runs", small_run("big", jobs=3, events=100, fsi=None))
        self.assertEqual(status, 400)
        self.assertIn("--max-cpus", out["error"])

    def test_the_queue_survives_a_server_restart(self):
        harness.need_binaries()
        gui = Gui(args=["--max-cpus", "1"])
        self.addCleanup(gui.close)
        gui.api("api/runs", small_run("first", model="3N", events=30000, jobs=1, seed=1))
        gui.api("api/runs", small_run("second", events=500, jobs=1, seed=2, fsi=None))
        gui.stop()
        gui.start(["--max-cpus", "1"])
        self.assertEqual(gui.wait("first")["state"], "done")
        self.assertEqual(gui.wait("second")["state"], "done")


class TestEstimate(GuiCase):
    def test_check_reports_time_disk_and_statistics(self):
        r = self.gui.api("api/check", small_run("est", events=200000, jobs=4))
        self.assertTrue(r["ok"], r["problems"])
        e = r["estimate"]
        self.assertGreater(e["seconds"], 0)
        # 2N with FSI: a few hundred bytes per event in each file.
        self.assertTrue(100 < e["bytes_per_event"]["generate"] < 1000, e)
        self.assertTrue(200 < e["bytes_per_event"]["fsi"] < 1500, e)
        self.assertEqual(e["peak_bytes"], 2 * e["bytes"])
        self.assertIn("job", r["estimate_text"])
        self.assertIn("N_eff", r["estimate_text"])

    def test_estimate_is_close_to_a_real_run(self):
        e = self.gui.api("api/check", small_run("est2", events=30000, jobs=2))["estimate"]
        self.gui.api("api/runs", small_run("est2", events=30000, jobs=2))
        d = self.gui.wait("est2")
        size = d["outputs"]["generate"]["size"] + d["outputs"]["fsi"]["size"]
        self.assertLess(abs(size - e["bytes"]) / size, 0.5)
        took = d["finished"] - d["started"]
        self.assertLess(e["seconds"], 5 * took + 5)

    def test_impossible_disk_use_is_refused_before_it_starts(self):
        status, out = self.gui.request("api/runs", small_run("huge", events=10**15, jobs=1))
        self.assertEqual(status, 400)
        self.assertIn("disk space", out["error"])
        self.assertFalse((self.gui.runs_dir / "huge").exists())


class TestFiles(GuiCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.gui.api("api/runs", small_run("files", events=3000, jobs=2, seed=2))
        cls.details = cls.gui.wait("files")

    def fetch(self, path, headers=None):
        import urllib.request
        req = urllib.request.Request(f"http://localhost:{self.gui.port}/{path}",
                                     headers=dict(headers or {}, **{"X-GCF-Token": self.gui.token}))
        with urllib.request.urlopen(req, timeout=60) as r:
            return r.status, dict(r.headers), r.read()

    def test_event_files_download_whole_and_in_part(self):
        path = Path(self.details["outputs"]["fsi"]["path"])
        status, headers, data = self.fetch(f"api/runs/files/{path.name}")
        self.assertEqual(status, 200)
        self.assertEqual(data, path.read_bytes())
        self.assertIn("attachment", headers["Content-Disposition"])
        status, headers, data = self.fetch(f"api/runs/files/{path.name}", {"Range": "bytes=100-199"})
        self.assertEqual(status, 206)
        self.assertEqual(data, path.read_bytes()[100:200])
        self.assertEqual(headers["Content-Range"], f"bytes 100-199/{path.stat().st_size}")

    def test_job_logs_can_be_read_but_nothing_outside_the_run(self):
        self.assertEqual(self.gui.request("api/runs/files/logs/generate_000.log")[0], 200)
        for bad in ("api/runs/files/logs/..%2F..%2F..%2Fserver.log", "api/runs/files/..%2Fserver.log",
                    "api/runs/..%2F..%2Fserver.log/x", "api/runs/files/logs/generate_000.exit"):
            self.assertIn(self.gui.request(bad)[0], (400, 404), bad)

    def test_unauthorized_requests_get_nothing(self):
        for path in ("api/runs", "api/runs/files/config.toml", "static/runs.js", "api/runs/files/logs/generate_000.log"):
            self.assertEqual(self.gui.request(path, token="wrong")[0], 403, path)

    def test_examples_and_saved_configs_load(self):
        from gcf import cli
        for name, text in cli.EXAMPLES.items():
            cfg = self.gui.api("api/parse-config", {"text": text})
            self.assertIn("run", cfg, name)
        saved = (self.gui.runs_dir / "files/config.toml").read_text()
        self.assertEqual(self.gui.api("api/parse-config", {"text": saved})["generate"]["events"], 3000)
        status, out = self.gui.request("api/parse-config", {"text": "this is not toml ["})
        self.assertEqual(status, 400)
        self.assertIn("TOML", out["error"])
        status, out = self.gui.request("api/parse-config", {"text": "[generate]\nbogus = 1\n"})
        self.assertEqual(status, 400)
        self.assertIn("bogus", out["error"])


class TestCommands(unittest.TestCase):
    def setUp(self):
        self.home = tempfile.mkdtemp(prefix="gcf-cli-")
        self.addCleanup(shutil.rmtree, self.home, ignore_errors=True)
        self.env = dict(os.environ, XDG_CONFIG_HOME=self.home)

    def gcf(self, *args):
        return subprocess.run([sys.executable, str(REPO / "gcf"), *map(str, args)], capture_output=True,
                              text=True, env=self.env, timeout=120)

    def test_background_server_lifecycle(self):
        harness.need_binaries()
        port = harness.free_port()
        runs = Path(self.home) / "runs"
        out = self.gcf("gui", "--background", "--port", port, "--runs-dir", runs)
        self.assertEqual(out.returncode, 0, out.stderr)
        self.addCleanup(self.gcf, "gui", "--stop", "--port", port)
        token = re.search(r"token=(\S+)", out.stdout).group(1)
        self.assertIn("log out whenever you like", out.stdout)
        # The token file is private to the user.
        mode = (Path(self.home) / "gcf/gui-token").stat().st_mode & 0o777
        self.assertEqual(mode, 0o600)
        again = self.gcf("gui", "--background", "--port", port, "--runs-dir", runs)
        self.assertIn("Already running", again.stdout)
        self.assertIn(token, again.stdout)
        status = self.gcf("gui", "--status", "--port", port)
        self.assertIn(token, status.stdout)
        self.assertIn("is running", status.stdout)
        clash = self.gcf("gui", "--port", port, "--runs-dir", runs)
        self.assertEqual(clash.returncode, 2)
        self.assertIn("already running", clash.stderr)
        self.assertEqual(self.gcf("gui", "--stop", "--port", port).returncode, 0)
        deadline = time.time() + 10
        while time.time() < deadline and "No gcf gui" not in self.gcf("gui", "--status", "--port", port).stdout:
            time.sleep(0.2)
        self.assertIn("No gcf gui", self.gcf("gui", "--status", "--port", port).stdout)
        # A new token invalidates the old address.
        out = self.gcf("gui", "--background", "--port", port, "--runs-dir", runs, "--new-token")
        self.assertNotIn(token, out.stdout)

    def test_detached_run_and_status(self):
        harness.need_binaries()
        from gcf import config
        runs = Path(self.home) / "runs"
        path = Path(self.home) / "c.toml"
        cfg = small_run("detached", events=2000, jobs=2)
        cfg["run"]["output_dir"] = str(runs)
        path.write_text(config.dumps(config.validate(cfg)))
        out = self.gcf("run", "--detach", path)
        self.assertEqual(out.returncode, 0, out.stderr)
        self.assertIn("on its own", out.stdout)
        deadline = time.time() + 120
        while time.time() < deadline:
            status = self.gcf("status", "--runs-dir", runs)
            if re.search(r"detached\s+done", status.stdout):
                break
            time.sleep(0.5)
        self.assertRegex(status.stdout, r"detached\s+done")

    def test_doctor(self):
        out = self.gcf("doctor")
        self.assertIn("gcf_generate", out.stdout)
        self.assertIn("runs survive logout", out.stdout)

    def test_runs_directory_default(self):
        from gcf import runner
        self.assertEqual(runner.default_runs_dir(), REPO / "runs")
        old = os.environ.get("GCF_RUNS_DIR")
        os.environ["GCF_RUNS_DIR"] = "/some/where"
        try:
            self.assertEqual(runner.default_runs_dir(), Path("/some/where"))
        finally:
            del os.environ["GCF_RUNS_DIR"]
            if old:
                os.environ["GCF_RUNS_DIR"] = old


class TestValidation(unittest.TestCase):
    """What the config refuses as you type, and that the binary agrees."""

    def problem(self, generate=None, fsi=None, stages=("generate", "fsi")):
        from gcf import config
        raw = {"run": {"stages": list(stages)}}
        if "generate" in stages:
            raw["generate"] = {"model": "2N", **(generate or {})}
        if "fsi" in stages:
            raw["fsi"] = dict(fsi or {})
        try:
            config.validate(raw)
        except config.ConfigError as e:
            return str(e)
        return None

    def test_values_outside_their_range_are_refused(self):
        for gen, key in (({"estar": -1.0}, "estar"), ({"estar": 5.0}, "estar"), ({"sigma_cm": -0.1}, "sigma_cm"),
                         ({"sigma_cm": 50.0}, "sigma_cm"), ({"beam": 0.0}, "beam"), ({"beam": 1e6}, "beam"),
                         ({"events": 0}, "events"), ({"cpn1": -1.0}, "cpn1"), ({"sigma_e": 0.0}, "sigma_e"),
                         ({"kcut": 2.0}, "kcut"), ({"q2": [-1.0, 5.0]}, "q2")):
            message = self.problem(gen)
            self.assertIsNotNone(message, gen)
            self.assertIn(f"[generate] {key}", message)
        self.assertIn("[generate] theta_e", self.problem({"model": "3N", "theta_e": [10.0, 200.0]}))
        self.assertIn("nucleon_mfp_scale", self.problem(fsi={"nucleon_mfp_scale": 0.0}))

    def test_values_in_range_pass(self):
        for gen in ({"estar": 0.02}, {"estar": 0.0}, {"sigma_cm": 0.0}, {"sigma_cm": 0.15}, {"beam": 437.0},
                    {"cpp0": 0.0}, {"kcut": 0.3}, {"q2": [1.0, 10.0]}):
            self.assertIsNone(self.problem(gen), gen)

    def test_settings_that_would_silently_cancel_each_other(self):
        """Each is valid alone; together one of them would be ignored (#123)."""
        self.assertIn("[generate] kcut", self.problem({"kcut": 0.4, "onset": "local"}))
        self.assertIn("[generate] sigma_cm", self.problem({"model": "3N", "sigma_cm": 0.2, "cm": False}))
        self.assertIn("[generate] kF_fg", self.problem({"model": "MF", "kF_fg": 0.2, "fermi_gas": "local"}))
        self.assertIn("[generate] prel_is", self.problem({"light_cone": True, "prel_is": True}))
        self.assertIn("[generate] onset", self.problem({"light_cone": True, "onset": "local"}))

    def test_what_means_nothing_for_a_nucleus_is_refused(self):
        d = {"target": [2, 1]}
        for key in ({"sigma_cm": 0.1}, {"estar": 0.02}, {"sigma_e": 0.01}, {"cpn1": 50.0}, {"randomize": True}):
            self.assertIn("deuteron", self.problem({**d, **key}), key)
        he3 = {"model": "3N", "target": [3, 2]}
        self.assertIn("A = 3", self.problem({**he3, "sigma_cm": 0.1}))
        self.assertIn("A = 3", self.problem({**he3, "estar": 0.02}))

    def test_all_zero_contacts_are_refused(self):
        self.assertIn("zero", self.problem({"cpp0": 0.0, "cnn0": 0.0, "cpn0": 0.0, "cpn1": 0.0}))
        self.assertIsNone(self.problem({"cpp0": 0.0, "cnn0": 0.0}))

    def test_mf_is_offered_only_where_it_works(self):
        self.assertIsNone(self.problem({"model": "MF", "target": [12, 6]}))
        self.assertIn("target", self.problem({"model": "MF", "target": [56, 26]}))

    def test_hA_has_no_compound_channel(self):
        """hA reads the flag but never uses it (#177)."""
        self.assertIn("nucleon.compound", self.problem(fsi={"model": "hA", "off": ["nucleon.compound"]}))
        self.assertIsNone(self.problem(fsi={"model": "hN", "off": ["nucleon.compound"]}))

    def test_the_binary_refuses_the_same(self):
        harness.need_binaries()
        env = harness.runner.environment()
        generate = harness.runner.binary("gcf_generate")
        for args in (["--model", "2N", "--kcut", "0.4", "--onset", "local"], ["--model", "3N", "--sigma-cm", "0.2", "--no-cm"],
                     ["--model", "MF", "--kF-fg", "0.2", "--fermi-gas", "local"], ["--model", "2N", "--beam", "1e6"],
                     ["--model", "2N", "--sigma-cm", "50"], ["--model", "2N", "-A", "2", "-Z", "1", "--estar", "0.02"],
                     ["--model", "3N", "-A", "3", "-Z", "2", "--sigma-cm", "0.1"]):
            full = args if "--beam" in args else ["--beam", "6", *args]
            out = subprocess.run([generate, "--check", *full], capture_output=True, text=True, env=env, timeout=120)
            self.assertEqual(out.returncode, 1, (args, out.stdout[-200:], out.stderr[-200:]))
        fsi = harness.runner.binary("gcf_fsi")
        out = subprocess.run([fsi, "--input", "x.root", "--output", "y.root", "--fsi-model", "hA", "--channel", "nucleon.compound=0"],
                             capture_output=True, text=True, env=env, timeout=60)
        self.assertEqual(out.returncode, 1)
        self.assertIn("nucleon.compound", out.stderr)


class TestFsiPreflight(GuiCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.gui.api("api/runs", small_run("base", events=1500, jobs=1, seed=3))
        cls.details = cls.gui.wait("base")

    def fsi_config(self, name, path, **fsi):
        return {"run": {"name": name, "stages": ["fsi"], "seed": 5}, "fsi": {"input": str(path), "model": "hA", **fsi}}

    def test_a_file_that_had_fsi_is_caught_by_the_check(self):
        cfg = self.fsi_config("again", self.details["outputs"]["fsi"]["path"])
        r = self.gui.api("api/check", cfg)
        self.assertFalse(r["ok"])
        self.assertIn("already been through FSI", r["problems"][0])
        self.assertTrue(r["problems"][0].startswith("[fsi] input"))
        status, out = self.gui.request("api/runs", cfg)
        self.assertEqual(status, 400)
        self.assertFalse((self.gui.runs_dir / "again").exists())

    def test_replacing_earlier_fsi_when_asked(self):
        cfg = self.fsi_config("redo", self.details["outputs"]["fsi"]["path"], overwrite=True)
        self.assertTrue(self.gui.api("api/check", cfg)["ok"])
        self.gui.api("api/runs", cfg)
        d = self.gui.wait("redo")
        self.assertEqual(d["state"], "done", d.get("error"))
        self.assertEqual(self.gui.api("api/runs/redo/results", {"stage": "fsi"})["events"], 1500)

    def test_pwia_file_passes(self):
        self.assertTrue(self.gui.api("api/check", self.fsi_config("ok", self.details["outputs"]["generate"]["path"]))["ok"])

    def test_missing_file_and_missing_branches_are_explained(self):
        r = self.gui.api("api/check", self.fsi_config("missing", "/no/such/file.root"))
        self.assertFalse(r["ok"])
        self.assertTrue(r["problems"][0].startswith("[fsi] input"), r["problems"])
        import uproot
        import numpy as np
        path = Path(self.gui.tmp) / "bare.root"
        with uproot.recreate(path) as f:
            f["events"] = {"weight": np.ones(5)}
        r = self.gui.api("api/check", self.fsi_config("bare", path))
        self.assertFalse(r["ok"])
        self.assertIn("nuc_pdg", r["problems"][0])

    def test_another_tree_name(self):
        import numpy as np
        import uproot
        path = Path(self.gui.tmp) / "other.root"
        with uproot.open(self.details["outputs"]["generate"]["path"]) as f:
            a = f["events"].arrays(["nuc_pdg", "nuc_px", "nuc_py", "nuc_pz", "nuc_E"])
        with uproot.recreate(path) as f:
            f.mktree("mytree", {k: "var * int32" if k == "nuc_pdg" else "var * float64" for k in a.fields} | {"weight": "float64"})
            f["mytree"].extend({k: a[k] for k in a.fields} | {"weight": np.ones(len(a))})
        wrong = self.gui.api("api/check", self.fsi_config("t1", path, target=[12, 6]))
        self.assertFalse(wrong["ok"])
        cfg = self.fsi_config("t2", path, target=[12, 6], tree="mytree")
        self.assertTrue(self.gui.api("api/check", cfg)["ok"])
        self.gui.api("api/runs", cfg)
        d = self.gui.wait("t2")
        self.assertEqual(d["state"], "done", d.get("error"))

    def test_an_interaction_without_contacts_points_at_that_field(self):
        r = self.gui.api("api/check", small_run("n3lo", nn="N3LO", events=100, fsi=None))
        self.assertFalse(r["ok"])
        self.assertTrue(r["problems"][0].startswith("[generate] nn"), r["problems"])


class TestScans(GuiCase):
    def scan(self, name, axes, **generate):
        cfg = small_run(name, events=generate.pop("events", 1500), jobs=1, fsi=None, seed=generate.pop("seed", 21), **generate)
        return cfg, axes

    def wait_scan(self, name, timeout=240):
        deadline = time.time() + timeout
        while time.time() < deadline:
            t = self.gui.api(f"api/scans/{name}")
            if not any(p["state"] in ("running", "merging", "starting", "queued") for p in t["points"]):
                return t
            time.sleep(0.5)
        self.fail("scan still running")

    def test_values_are_read_as_written_or_as_a_range(self):
        from gcf import scan
        self.assertEqual(scan.parse_values("0.1, 0.15 0.2", "float"), [0.1, 0.15, 0.2])
        self.assertEqual(scan.parse_values("0.1:0.3:3", "float"), [0.1, 0.2, 0.3])
        self.assertEqual(scan.parse_values("1:4:4", "int"), [1, 2, 3, 4])
        self.assertEqual(scan.parse_values("cc1,cc2", "str"), ["cc1", "cc2"])
        for bad, kind in (("", "float"), ("a,b", "float"), ("1.5", "int"), ("1,1", "float"), ("1:2:1", "float")):
            with self.assertRaises(ValueError, msg=bad):
                scan.parse_values(bad, kind)

    def test_points_and_names(self):
        from gcf import scan
        raw = small_run("s", events=100, fsi=None)
        pts = scan.expand(raw, [{"path": "generate.sigma_cm", "values": "0.1,0.2"}, {"path": "generate.nn", "values": "AV18,N2LO"}])
        self.assertEqual([p[0] for p in pts], ["s_sigma_cm-0.1_nn-AV18", "s_sigma_cm-0.1_nn-N2LO", "s_sigma_cm-0.2_nn-AV18", "s_sigma_cm-0.2_nn-N2LO"])
        self.assertEqual(pts[3][2]["generate"]["nn"], "N2LO")
        self.assertEqual(pts[3][2]["run"]["scan"], "s")
        self.assertEqual({p[2]["run"]["seed"] for p in pts}, {1234})
        unseeded = small_run("u", events=100, fsi=None, seed=0)
        self.assertEqual(len({p[2]["run"]["seed"] for p in scan.expand(unseeded, [{"path": "generate.sigma_cm", "values": "0.1,0.2"}])}), 1)
        chan = scan.expand(small_run("c", events=100), [{"path": "fsi.scale:pion.absorption", "values": "0.5,1.5"}])
        self.assertEqual(chan[1][2]["fsi"]["scale"], {"pion.absorption": 1.5})
        for bad in ("generate.model", "generate.events", "fsi.input", "run.name", "fsi.scale:nucleon.compound", "bogus.x", "generate.bogus"):
            with self.assertRaises(ValueError, msg=bad):
                scan.option(bad)
        with self.assertRaises(ValueError):
            scan.expand(raw, [{"path": "generate.sigma_cm", "values": "0:1:30"}, {"path": "generate.estar", "values": "0:1:30"}])

    def test_a_contact_scan_on_shared_events_scales_as_it_should(self):
        """The whole chain, with a known answer: pp and nn stay, so the yield is linear in C_pn(S=1)."""
        base = self.gui.api("api/nuclear?model=2N&A=12&Z=6&nn=")["contacts"]["cpn1"][0]
        cfg = small_run("c1", events=5000, jobs=1, fsi=None, seed=33)
        r = self.gui.api("api/scans", {"config": cfg, "axes": [{"path": "generate.cpn1", "values": f"{base}, {2 * base}"}]})
        self.assertEqual(r["points"], 2)
        t = self.wait_scan("c1")
        self.assertEqual([p["state"] for p in t["points"]], ["done", "done"])
        a, b = t["points"]
        ratio = b["yield"] / a["yield"]
        # The pn weights nearly double: 1 + (pn share) x (1 - small S=0 part).
        self.assertTrue(1.5 < ratio < 2.0, ratio)
        self.assertEqual(a["events"], 5000)
        self.assertEqual(a["attempts"], b["attempts"], "the same events")
        self.assertTrue(any("share a seed" in n for n in t["notes"]))
        listed = self.gui.api("api/runs")
        self.assertEqual(next(s for s in listed["scans"] if s["name"] == "c1")["states"], {"done": 2})
        self.assertEqual(next(r for r in listed["runs"] if r["name"] == a["name"])["scan"], "c1")

    def test_a_mistake_in_one_point_starts_none(self):
        cfg = small_run("bad1", events=100, jobs=1, fsi=None)
        status, out = self.gui.request("api/scans", {"config": cfg, "axes": [{"path": "generate.estar", "values": "0.01, -1"}]})
        self.assertEqual(status, 400)
        self.assertIn("estar", out["error"])
        self.assertFalse(list(self.gui.runs_dir.glob("bad1*")))
        status, out = self.gui.request("api/scans", {"config": cfg, "axes": [{"path": "generate.nn", "values": "AV18, N3LO"}]})
        self.assertEqual(status, 400)
        self.assertIn("contacts", out["error"])
        self.assertFalse(list(self.gui.runs_dir.glob("bad1*")))

    def test_preview_and_duplicates(self):
        cfg = small_run("pv", events=100, jobs=1, fsi=None)
        axes = [{"path": "generate.sigma_cm", "values": "0.1:0.2:3"}]
        p = self.gui.api("api/scans/preview", {"config": cfg, "axes": axes})
        self.assertEqual((p["points"], p["events"]), (3, 300))
        self.assertEqual(self.gui.api("api/scans", {"config": cfg, "axes": axes})["points"], 3)
        status, out = self.gui.request("api/scans", {"config": cfg, "axes": axes})
        self.assertEqual(status, 400)
        self.assertIn("already exists", out["error"])
        self.wait_scan("pv")
        self.gui.api("api/scans/pv/delete", {})
        self.assertEqual(self.gui.request("api/scans/pv")[0], 404)
        self.assertFalse(list(self.gui.runs_dir.glob("pv_*")))

    def test_points_respect_the_cpu_limit(self):
        gui = Gui(args=["--max-cpus", "1"])
        self.addCleanup(gui.close)
        cfg = small_run("lim", model="3N", events=3000, jobs=1, fsi=None)
        gui.api("api/scans", {"config": cfg, "axes": [{"path": "generate.sigma_cm", "values": "0.1,0.15,0.2"}]})
        peak = 0
        deadline = time.time() + 120
        while time.time() < deadline:
            states = [p["state"] for p in gui.api("api/scans/lim")["points"]]
            peak = max(peak, sum(s in ("running", "merging", "starting") for s in states))
            if all(s == "done" for s in states):
                break
            time.sleep(0.1)
        self.assertEqual(states, ["done"] * 3)
        self.assertLessEqual(peak, 1)

    def test_scan_from_the_command_line(self):
        import tempfile as tf
        from gcf import config
        with tf.TemporaryDirectory() as tmp:
            cfg = small_run("cl", events=500, jobs=1, fsi=None)
            cfg["run"]["output_dir"] = str(Path(tmp) / "runs")
            path = Path(tmp) / "c.toml"
            path.write_text(config.dumps(config.validate(cfg)))
            out = subprocess.run([sys.executable, str(REPO / "gcf"), "scan", str(path), "--set", "generate.sigma_cm=0.1,0.2", "-q"],
                                 capture_output=True, text=True, timeout=120)
            self.assertEqual(out.returncode, 0, out.stderr)
            self.assertTrue((Path(tmp) / "runs/cl_sigma_cm-0.2/cl_sigma_cm-0.2_pwia.root").is_file())
            self.assertTrue((Path(tmp) / "runs/.scans/cl.json").is_file())


class TestOtherGenerators(GuiCase):
    def config(self, name, program, events=4000, jobs=2, **options):
        return {"run": {"name": name, "stages": ["program"], "jobs": jobs, "seed": 11},
                "program": {"program": program, "events": events, "target": [12, 6], **options}}

    def run_program(self, program, suffix="", **options):
        harness.need_suite(program)
        name = f"g_{program}{suffix}"
        from gcf import schema
        if schema.PROGRAMS[program].beam:
            # Above the J/psi threshold for genElectro.
            options.setdefault("beam", 11.0)
        self.gui.api("api/runs", self.config(name, program, **options))
        d = self.gui.wait(name)
        self.assertEqual(d["state"], "done", d.get("error"))
        self.assertTrue(d["outputs"]["program"]["exists"])
        return d, self.gui.api(f"api/runs/{name}/results", {"stage": "program"})

    def test_every_generator_runs_and_merges(self):
        for program in ("genDIS", "genCCQE", "genElectro", "genIon", "sfSampler", "sfSamplerLC"):
            with self.subTest(program=program):
                d, r = self.run_program(program)
                self.assertGreater(r["events"], 0)
                self.assertEqual(r["attempts"], 4000)
                self.assertIsNotNone(r["n_eff"])
                self.assertEqual(d["stages"]["program"]["done"], 4000)

    def test_genphoto_options_reach_it(self):
        d, r = self.run_program("genPhoto", reaction="phi", fixed_beam=8.0)
        self.assertIn("-R phi -B 8.0", (self.gui.runs_dir / "g_genPhoto/logs/program_000.log").read_text().splitlines()[0])
        self.assertGreater(r["events"], 0)

    def test_gendis_seed_makes_it_repeatable(self):
        import numpy as np
        import uproot
        harness.need_suite("genDIS")
        for name in ("rep1", "rep2"):
            self.gui.api("api/runs", self.config(name, "genDIS", events=5000, beam=6.0))
        a, b = self.gui.wait("rep1"), self.gui.wait("rep2")
        wa = uproot.open(a["outputs"]["program"]["path"])["genT"]["weight"].array(library="np")
        wb = uproot.open(b["outputs"]["program"]["path"])["genT"]["weight"].array(library="np")
        self.assertTrue(np.array_equal(wa, wb))

    def test_options_of_another_generator_are_refused(self):
        status, out = self.gui.request("api/check", self.config("x", "genCCQE", reaction="phi"))
        self.assertEqual(status, 400)
        self.assertIn("reaction", out["error"])
        status, out = self.gui.request("api/check", {"run": {"name": "x", "stages": ["program", "fsi"]}, "program": {"program": "genDIS"}})
        self.assertEqual(status, 400)

    def test_plots_and_comparing_say_why_not(self):
        harness.need_suite("genCCQE")
        d, _ = self.run_program("genCCQE", suffix="_noplots")
        status, out = self.gui.request("api/runs/g_genCCQE_noplots/quicklook", {"stage": "program"})
        self.assertEqual(status, 400)
        self.assertIn("its own", out["error"])
        self.assertEqual(d["state"], "done")

    def test_no_stored_events_are_explained(self):
        """genElectro below the J/psi threshold: every attempt has zero weight."""
        harness.need_suite("genElectro")
        r = self.gui.api("api/check", self.config("thr", "genElectro", beam=6.0))
        self.assertTrue(r["ok"], r["problems"])
        self.assertIn("stored no events", r["estimate_text"])
        d, s = self.run_program("genElectro", suffix="_thr", beam=6.0)
        self.assertEqual(s["events"], 0)
        self.assertTrue(any(level == "bad" and "No attempt" in text for level, text in s["notes"]), s["notes"])

    def test_check_estimates_a_generator(self):
        harness.need_suite("genDIS")
        r = self.gui.api("api/check", self.config("e", "genDIS", events=100000, beam=6.0))
        self.assertTrue(r["ok"], r["problems"])
        self.assertGreater(r["estimate"]["bytes"], 0)
        self.assertIn("job", r["estimate_text"])

    def test_a_target_it_has_no_table_for_is_caught(self):
        harness.need_suite("genCCQE")
        cfg = self.config("t", "genCCQE")
        cfg["program"]["target"] = [13, 6]
        r = self.gui.api("api/check", cfg)
        self.assertFalse(r["ok"])
        self.assertIn("[program] target", r["problems"][0])


class TestBrowser(GuiCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.driver = harness.browser()
        cls.addClassCleanup(cls.driver.quit)

    def open(self):
        from selenium.webdriver.common.by import By
        from selenium.webdriver.support.ui import WebDriverWait
        self.driver.get(self.gui.url)
        self.driver.execute_script("localStorage.clear()")   # each test starts from the default form
        self.driver.get(self.gui.url)
        WebDriverWait(self.driver, 20).until(lambda d: d.find_elements(By.CSS_SELECTOR, "#grid-run .field"))

    def set_value(self, field_id, value):
        from selenium.webdriver.common.by import By
        field = self.driver.find_element(By.ID, field_id)
        field.clear()
        field.send_keys(str(value))

    def test_start_message_stays_visible(self):
        """'Started run' used to be cleared by the same call that reset the form."""
        from selenium.webdriver.common.by import By
        from selenium.webdriver.support.ui import Select, WebDriverWait
        self.open()
        Select(self.driver.find_element(By.ID, "f-generate-model")).select_by_value("2N")
        self.set_value("f-generate-events", 500)
        name = self.driver.find_element(By.ID, "f-run-name").get_attribute("value")
        time.sleep(0.5)
        self.driver.find_element(By.ID, "btn-start").click()
        WebDriverWait(self.driver, 60).until(lambda d: d.find_element(By.ID, "d-title").text == name)
        time.sleep(0.8)        # past the debounced validation of the reset form
        self.assertIn(f"Submitted run {name}", self.driver.find_element(By.ID, "form-msg").text)
        self.gui.wait(name)

    def test_outputs_show_numbers_and_warnings(self):
        from selenium.webdriver.common.by import By
        from selenium.webdriver.support.ui import WebDriverWait
        self.gui.api("api/runs", small_run("nums", model="3N", events=300, jobs=1, fsi=None))
        self.gui.wait("nums")
        self.open()
        self.driver.execute_script("showDetails('nums')")
        WebDriverWait(self.driver, 60).until(lambda d: d.find_elements(By.CSS_SELECTOR, "#d-body .numbers"))
        text = self.driver.find_element(By.ID, "d-body").text
        self.assertIn("N_eff", text)
        self.assertIn("Σw / attempts", text)
        self.assertIn("a few events dominate", text)
        self.assertIn("placeholder (#73)", text)

    def test_compare_dialog(self):
        from selenium.webdriver.common.by import By
        from selenium.webdriver.support.ui import WebDriverWait
        for name, seed in (("dlg_a", 5), ("dlg_b", 6)):
            self.gui.api("api/runs", small_run(name, events=3000, jobs=2, seed=seed))
        self.gui.wait("dlg_a")
        self.gui.wait("dlg_b")
        self.open()
        WebDriverWait(self.driver, 20).until(lambda d: d.find_elements(By.CSS_SELECTOR, "tr.run"))
        self.driver.execute_script("openCompare('dlg_a')")
        self.driver.execute_script("""const b = document.getElementById('cmp-b-run'); b.value = 'dlg_b';
                                      b.dispatchEvent(new Event('change'));""")
        self.driver.execute_script("[...document.querySelectorAll('#compare-dialog button')].find(x => x.textContent == 'Compare').click()")
        WebDriverWait(self.driver, 200).until(lambda d: d.find_elements(By.CSS_SELECTOR, "#compare-out table"))
        text = self.driver.find_element(By.ID, "compare-out").text
        self.assertIn("Weighted p", text)
        self.assertIn("dlg_b (with FSI)", text)

    def test_physics_tab_shows_table_values_and_sets_options(self):
        from selenium.webdriver.common.by import By
        from selenium.webdriver.support.ui import Select, WebDriverWait
        self.open()
        Select(self.driver.find_element(By.ID, "f-generate-model")).select_by_value("2N")
        self.driver.execute_script("document.getElementById('tab-physics').click()")
        WebDriverWait(self.driver, 30).until(lambda d: "contacts" in d.find_element(By.ID, "physics-head").text)
        head = self.driver.find_element(By.ID, "physics-head").text
        self.assertIn("12C", head)
        self.assertIn("k-space", head)
        field = self.driver.find_element(By.ID, "f-generate-cpn1")
        self.assertEqual(field.get_attribute("placeholder"), "95.256")
        self.assertIn("95.256 ± 10.62", self.driver.find_element(By.CSS_SELECTOR, '.field[data-key="cpn1"] .tableval').text)
        self.assertEqual(self.driver.find_element(By.ID, "f-generate-sigma_cm").get_attribute("placeholder"), "0.15")
        # Another target or interaction updates the values.
        self.driver.execute_script("""const n = document.getElementById('f-generate-nn'); n.value = 'N2LO';
                                      n.dispatchEvent(new Event('change', {bubbles: true}));""")
        WebDriverWait(self.driver, 30).until(lambda d: d.find_element(By.ID, "f-generate-cpn1").get_attribute("placeholder") == "63.858")
        self.driver.execute_script("""const n = document.getElementById('f-generate-nn'); n.value = '';
                                      n.dispatchEvent(new Event('change', {bubbles: true}));""")
        WebDriverWait(self.driver, 30).until(lambda d: d.find_element(By.ID, "f-generate-cpn1").get_attribute("placeholder") == "95.256")
        # A value set on the tab reaches the config, and the tab counts it.
        self.set_value("f-generate-cpn1", 50)
        WebDriverWait(self.driver, 10).until(lambda d: "1 changed" in d.find_element(By.ID, "tab-count-physics").text)
        cfg = self.driver.execute_script("return collect()")
        self.assertEqual(cfg["generate"]["cpn1"], 50)
        # The fields of other models are not shown, and a 3N setup says why.
        self.driver.execute_script("""const m = document.getElementById('f-generate-model'); m.value = '3N';
                                      m.dispatchEvent(new Event('change', {bubbles: true}));""")
        self.assertNotIn("hidden", self.driver.find_element(By.ID, "physics-note").get_attribute("class"))
        self.assertIn("hidden", self.driver.find_element(By.CSS_SELECTOR, '.field[data-key="cpn1"]').get_attribute("class"))
        self.assertNotIn("hidden", self.driver.find_element(By.CSS_SELECTOR, '.field[data-key="estar"]').get_attribute("class"))
        self.driver.execute_script("resetForm(); showFormTab('setup')")

    def test_channel_scale_reaches_the_config(self):
        from selenium.webdriver.common.by import By
        self.open()
        self.driver.execute_script("resetForm()")
        self.assertTrue(self.driver.find_element(By.ID, "scale-pion.absorption").is_displayed())
        self.set_value("scale-pion.absorption", 0.5)
        self.assertEqual(self.driver.execute_script("return collect()")["fsi"]["scale"], {"pion.absorption": 0.5})
        self.driver.execute_script("resetForm()")
        self.assertNotIn("scale", self.driver.execute_script("return collect()")["fsi"])

    def test_first_visit_guide_and_small_test_run(self):
        from selenium.webdriver.common.by import By
        from selenium.webdriver.support.ui import WebDriverWait
        gui = Gui()
        self.addCleanup(gui.close)
        self.driver.get(gui.url)
        self.driver.execute_script("localStorage.clear()")
        self.driver.get(gui.url)
        WebDriverWait(self.driver, 20).until(lambda d: d.find_elements(By.CSS_SELECTOR, ".guide"))
        self.assertIn("GETTING STARTED", self.driver.find_element(By.ID, "runs").text.upper())
        self.driver.execute_script("fillTestRun()")
        name = self.driver.find_element(By.ID, "f-run-name").get_attribute("value")
        self.assertTrue(name.startswith("first-test-"))
        time.sleep(0.5)
        self.driver.execute_script("document.getElementById('btn-start').click()")
        WebDriverWait(self.driver, 60).until(lambda d: d.find_element(By.ID, "d-title").text == name)
        d = gui.wait(name)
        self.assertEqual(d["state"], "done", d.get("error"))
        WebDriverWait(self.driver, 120).until(lambda x: x.find_elements(By.CSS_SELECTOR, "#d-body img.plot"))
        self.driver.execute_script("localStorage.clear()")

    def test_stopped_run_offers_resume_and_add_fsi(self):
        from selenium.webdriver.common.by import By
        from selenium.webdriver.support.ui import WebDriverWait
        cfg = small_run("brk", events=40000, jobs=2, seed=3)
        cfg["run"]["chunks"] = 8
        self.gui.api("api/runs", cfg)
        deadline = time.time() + 60
        path = self.gui.runs_dir / "brk/status.json"
        while time.time() < deadline:
            if path.exists() and json.loads(path.read_text())["stages"]["fsi"]["finished"] >= 2:
                break
            time.sleep(0.2)
        self.gui.api("api/runs/brk/cancel", {})
        self.gui.wait("brk")
        self.open()
        self.driver.execute_script("showDetails('brk')")
        WebDriverWait(self.driver, 20).until(lambda x: "Resume" in x.find_element(By.ID, "d-body").text)
        text = self.driver.find_element(By.ID, "d-body").text
        self.assertRegex(text, r"Resume \(\d+ of 8 jobs left\)")
        self.assertIn("Use the finished jobs", text)
        self.driver.execute_script("[...document.querySelectorAll('#d-body button')].find(b => b.textContent.startsWith('Resume')).click()")
        WebDriverWait(self.driver, 20).until(lambda x: self.gui.api("api/runs/brk")["state"] != "cancelled")
        d = self.gui.wait("brk")
        self.assertEqual(d["state"], "done", d.get("error"))
        self.driver.execute_script("showDetails('brk')")
        WebDriverWait(self.driver, 60).until(lambda x: "Add FSI" in x.find_element(By.ID, "d-body").text)
        self.driver.execute_script("[...document.querySelectorAll('#d-body button')].find(b => b.textContent.startsWith('Add FSI')).click()")
        self.assertEqual(self.driver.execute_script("return collect()")["run"]["stages"], ["fsi"])
        self.assertTrue(self.driver.execute_script("return collect()")["fsi"]["input"].endswith("brk_pwia.root"))
        self.driver.execute_script("resetForm()")

    def test_page_reconnects_after_a_server_restart(self):
        from selenium.webdriver.common.by import By
        from selenium.webdriver.support.ui import WebDriverWait
        gui = Gui()
        self.addCleanup(gui.close)
        self.driver.get(gui.url)
        WebDriverWait(self.driver, 20).until(lambda d: d.find_elements(By.CSS_SELECTOR, "#grid-run .field"))
        gui.stop()
        WebDriverWait(self.driver, 15).until(lambda d: "Cannot reach" in d.find_element(By.ID, "banner").text)
        self.assertIn("Your runs carry on", self.driver.find_element(By.ID, "banner").text)
        gui.start()
        WebDriverWait(self.driver, 15).until(lambda d: "hidden" in d.find_element(By.ID, "banner").get_attribute("class"))

    def test_examples_fill_the_form(self):
        from selenium.webdriver.common.by import By
        from selenium.webdriver.support.ui import WebDriverWait
        self.open()
        self.driver.execute_script("""const s = document.getElementById('example-select'); s.value = '2N';
                                      s.dispatchEvent(new Event('change'));""")
        WebDriverWait(self.driver, 15).until(lambda d: d.find_element(By.ID, "f-generate-model").get_attribute("value") == "2N")
        self.assertEqual(self.driver.find_element(By.ID, "f-generate-beam").get_attribute("value"), "4.5")
        self.driver.execute_script("resetForm()")

    def test_scan_tab_from_form_to_table_and_chart(self):
        from selenium.webdriver.common.by import By
        from selenium.webdriver.support.ui import Select, WebDriverWait
        self.open()
        Select(self.driver.find_element(By.ID, "f-generate-model")).select_by_value("2N")
        self.driver.execute_script("document.querySelector(\"input[name=stages][value='generate']\").click()")
        self.set_value("f-generate-events", 2000)
        self.set_value("f-run-name", "uiscan")
        self.driver.execute_script("showFormTab('scan')")
        Select(self.driver.find_element(By.CSS_SELECTOR, ".scan-path")).select_by_value("generate.sigma_cm")
        self.driver.find_element(By.CSS_SELECTOR, ".scan-values").send_keys("0.1, 0.15, 0.2")
        WebDriverWait(self.driver, 10).until(lambda d: "3 runs" in d.find_element(By.ID, "scan-preview").text)
        self.assertIn("6,000 events", self.driver.find_element(By.ID, "scan-preview").text)
        self.driver.find_element(By.ID, "btn-scan").click()
        WebDriverWait(self.driver, 200).until(lambda d: len(d.find_elements(By.CSS_SELECTOR, ".scan-chart circle.scan-dot")) == 3)
        body = self.driver.find_element(By.ID, "d-body").text
        self.assertIn("Ratio to first", body)
        self.assertIn("3 done", body)
        self.assertIn("uiscan", self.driver.find_element(By.ID, "scans").text)
        # Its points are listed with the scan, not among the runs, unless asked for.
        self.assertFalse(self.driver.find_elements(By.CSS_SELECTOR, 'tr.run[data-run^="uiscan_"]'))
        self.driver.execute_script("document.getElementById('show-points').click()")
        self.assertEqual(len(self.driver.find_elements(By.CSS_SELECTOR, 'tr.run[data-run^="uiscan_"]')), 3)
        self.driver.execute_script("document.getElementById('show-points').click()")
        self.driver.save_screenshot(str(Path(self.gui.tmp) / "scan.png"))
        self.driver.execute_script("resetForm(); showFormTab('setup')")

    def test_other_generators_in_the_form(self):
        from selenium.webdriver.common.by import By
        self.open()
        self.driver.execute_script("document.querySelector(\"input[name=stages][value='program']\").click()")
        self.assertTrue(self.driver.find_element(By.ID, "fs-program").is_displayed())
        self.assertFalse(self.driver.find_element(By.ID, "fs-generate").is_displayed())
        self.assertIn("W²", self.driver.find_element(By.ID, "program-note").text)
        self.assertTrue(self.driver.find_element(By.ID, "tab-physics").get_attribute("disabled"))
        self.assertTrue(self.driver.find_element(By.ID, "tab-scan").get_attribute("disabled"))
        cfg = self.driver.execute_script("return collect()")
        self.assertEqual(cfg["run"]["stages"], ["program"])
        self.assertEqual(cfg["program"]["program"], "genDIS")
        self.driver.execute_script("""const p = document.getElementById('f-program-program'); p.value = 'genPhoto';
                                      p.dispatchEvent(new Event('change', {bubbles: true}));""")
        self.assertNotIn("hidden", self.driver.find_element(By.CSS_SELECTOR, '.field[data-section="program"][data-key="reaction"]').get_attribute("class"))
        self.assertIn("hidden", self.driver.find_element(By.CSS_SELECTOR, '.field[data-section="program"][data-key="beam"]').get_attribute("class"))
        self.driver.execute_script("resetForm()")

    def test_internal_run_options_are_not_in_the_form(self):
        from selenium.webdriver.common.by import By
        self.open()
        self.assertFalse(self.driver.find_elements(By.CSS_SELECTOR, '.field[data-key="scan"]'))
        self.assertFalse(self.driver.find_elements(By.CSS_SELECTOR, '.field[data-key="output_dir"]'))

    def test_a_bad_value_marks_its_field(self):
        from selenium.webdriver.common.by import By
        from selenium.webdriver.support.ui import WebDriverWait
        self.open()
        self.driver.execute_script("""const m = document.getElementById('f-generate-model'); m.value = '2N';
                                      m.dispatchEvent(new Event('change', {bubbles: true})); showFormTab('physics');""")
        self.assertEqual(self.driver.find_element(By.ID, "f-generate-estar").get_attribute("min"), "0")
        self.set_value("f-generate-estar", -1)
        WebDriverWait(self.driver, 10).until(lambda d: "invalid" in d.find_element(By.CSS_SELECTOR, '.field[data-key="estar"]').get_attribute("class"))
        self.assertIn("at least 0", self.driver.find_element(By.CSS_SELECTOR, '.field[data-key="estar"] .ferr').text)
        self.driver.execute_script("resetForm(); showFormTab('setup')")

    def test_run_row_keeps_keyboard_focus(self):
        """The table is rebuilt every few seconds; that must not drop the focus."""
        self.gui.api("api/runs", small_run("focus", events=200, jobs=1))
        self.gui.wait("focus")
        self.open()
        from selenium.webdriver.common.by import By
        from selenium.webdriver.support.ui import WebDriverWait
        WebDriverWait(self.driver, 20).until(lambda d: d.find_elements(By.CSS_SELECTOR, "tr.run"))
        self.driver.execute_script("document.querySelector('tr.run').focus()")
        time.sleep(3.5)        # the table refreshes every 3 s
        self.assertEqual(self.driver.execute_script("return document.activeElement.tagName"), "TR")


class TestCommandLine(unittest.TestCase):
    def test_quicklook_without_packages_is_a_clean_error(self):
        """No traceback when quick-look's Python packages are missing."""
        code = ("import sys; sys.modules['awkward'] = None; "
                f"sys.path.insert(0, {str(REPO / 'tools')!r}); "
                "from gcf.cli import main; sys.exit(main(['quicklook', 'no_such_file.root']))")
        out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
        self.assertEqual(out.returncode, 2, out.stderr)
        self.assertIn("quicklook needs", out.stderr)
        self.assertNotIn("Traceback", out.stderr)


class TestProvenance(unittest.TestCase):
    def test_version_stamp_matches_the_tree(self):
        """WriteVersion.cmake records the commit of the tree it runs on, and only
        rewrites the file when the version changes."""
        if not shutil.which("cmake"):
            self.skipTest("no cmake")
        head = subprocess.run(["git", "-C", str(REPO), "describe", "--always", "--tags"],
                              capture_output=True, text=True).stdout.strip()
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "gcf_version.h"
            command = ["cmake", f"-DSRC_ROOT={REPO / 'src'}", f"-DOUT={out}", "-P",
                       str(REPO / "src/programs/gcf_stages/WriteVersion.cmake")]
            subprocess.run(command, check=True)
            self.assertIn(head, out.read_text())
            stamp = out.stat().st_mtime_ns
            subprocess.run(command, check=True)
            self.assertEqual(out.stat().st_mtime_ns, stamp)


if __name__ == "__main__":
    unittest.main()
