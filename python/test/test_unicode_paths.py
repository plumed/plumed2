# File names with non-ASCII characters, passed to the kernel as UTF-8.
# On Windows the kernel converts them to UTF-16 (the narrow functions of the
# C library would read them in the ANSI code page of python.exe).
# The names are written with escapes, to keep this file ASCII.
import array
import os
import tempfile
import unittest

import plumed

NAME = "é中"  # e acute, and a CJK character


def run(directory, log):
    """A distance printed to a file whose name and directory are not ASCII."""
    colvar = os.path.join(directory, "colvar-" + NAME).replace("\\", "/")
    include = os.path.join(directory, "include-" + NAME + ".dat").replace("\\", "/")
    with open(include, "w", encoding="utf-8") as f:
        f.write("d: DISTANCE ATOMS=1,2\n")
    with plumed.Plumed() as p:
        p.cmd("setNatoms", 2)
        p.cmd("setLogFile", os.path.join(directory, log))
        p.cmd("init")
        p.cmd("readInputLine", "INCLUDE FILE=" + include)
        p.cmd("readInputLine", "PRINT ARG=d FILE=" + colvar)
        box = array.array("d", [10, 0, 0, 0, 10, 0, 0, 0, 10])
        virial = array.array("d", [0] * 9)
        forces = array.array("d", [0] * 6)
        positions = array.array("d", [0, 0, 0, 1, 2, 2])
        p.cmd("setStep", 0)
        p.cmd("setBox", box)
        p.cmd("setPositions", positions)
        p.cmd("setForces", forces)
        p.cmd("setVirial", virial)
        p.cmd("calc")
    return colvar


class TestUnicodePaths(unittest.TestCase):
    def test_input_output_and_backup(self):
        with tempfile.TemporaryDirectory(prefix="plumed-" + NAME + "-") as directory:
            log = "log-" + NAME
            colvar = run(directory, log)
            self.assertTrue(os.path.isfile(os.path.join(directory, log)))
            with open(colvar, encoding="utf-8") as f:
                lines = [line for line in f if not line.startswith("#")]
            self.assertEqual(len(lines), 1)
            self.assertAlmostEqual(float(lines[0].split()[1]), 3.0)
            # a second run backs up the existing files (renamed to bck.0.*)
            run(directory, log)
            names = os.listdir(directory)
            self.assertIn("bck.0.colvar-" + NAME, names)
            self.assertIn("bck.0.log-" + NAME, names)


if __name__ == "__main__":
    unittest.main()
