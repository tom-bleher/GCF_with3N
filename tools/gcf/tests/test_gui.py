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
