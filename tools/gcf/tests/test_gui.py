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
