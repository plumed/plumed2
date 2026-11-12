"""Check that Vim generation skips target executables that cannot run."""

from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

PACKAGE = Path(__file__).resolve().parents[1]
ROOT = (
    PACKAGE / "plumed_source"
    if (PACKAGE / "plumed_source").is_dir()
    else PACKAGE.parent
)


@unittest.skipIf(sys.platform == "win32", "the Makefile build is not used on Windows")
@unittest.skipUnless(shutil.which("make"), "make is unavailable")
class TestVimMake(unittest.TestCase):
    def test_cross_compilation_does_not_run_generator(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            vim = root / "vim"
            vim.mkdir()
            for name in ("Makefile", "vimsyntax.sh"):
                shutil.copy2(ROOT / "vim" / name, vim / name)
            shutil.copytree(ROOT / "vim" / "templates", vim / "templates")
            (root / "sourceme.sh").write_text("")
            binary = root / "src" / "lib" / "plumed"
            binary.parent.mkdir(parents=True)
            binary.write_text("#!/bin/sh\n: > invoked\nexit 19\n")
            binary.chmod(0o755)

            skipped = subprocess.run(
                ["make", "program_can_run=no"],
                cwd=vim,
                capture_output=True,
                text=True,
            )
            self.assertEqual(skipped.returncode, 0, skipped.stdout + skipped.stderr)
            self.assertIn("cannot run on this machine", skipped.stdout)
            self.assertFalse((vim / "invoked").exists())
            self.assertFalse((vim / ".syntax.stamp").exists())

            # Prove the same fixture reaches the failing executable when enabled.
            enabled = subprocess.run(
                ["make", "program_can_run=yes"],
                cwd=vim,
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(enabled.returncode, 0)
            self.assertTrue((vim / "invoked").exists())
            self.assertFalse((vim / ".syntax.stamp").exists())


if __name__ == "__main__":
    unittest.main()
