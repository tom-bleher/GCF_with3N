#!/usr/bin/env python3
"""Tests of the web interface, against a real `gcf gui` and the real stages.

    python3 -m unittest discover -s tools/gcf/tests -v

Browser tests are skipped without selenium and Firefox.
"""

import shutil
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
        self.assertIn(f"Started run {name}", self.driver.find_element(By.ID, "form-msg").text)
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
