# LOAD, from the python interface, an action compiled at runtime (with plumed mklib).
# The python interface loads the kernel with RTLD_LOCAL, so that its symbols are not
# visible to the libraries loaded later: this works only if the compiled library
# records its dependency on the kernel, which linkers using --as-needed (the default
# e.g. on Ubuntu) drop when the kernel precedes the objects in the link command.
import array
import os
from pathlib import Path
import re
import shlex
import shutil
import sys
import tempfile
import unittest

import plumed


@unittest.skipIf(sys.platform == "win32", "LOAD is not supported on Windows")
class Test(unittest.TestCase):
    def test_load_compiled_action(self):
        root = Path(plumed._guessplumedroot())
        options = root / "src" / "config" / "compile_options.sh"
        template = root / "src" / "colvar" / "Template.cpp"
        if not options.is_file() or not template.is_file():
            self.skipTest(
                "compile_options.sh or Template.cpp not found in " + str(root)
            )
        text = options.read_text()
        if "__PLUMED_HAS_DLOPEN" not in text:
            self.skipTest("PLUMED compiled without dlopen")
        compile_line = re.search(r"^compile='(.*?) -c ", text, re.M)
        self.assertIsNotNone(
            compile_line, "cannot read the compiler from " + str(options)
        )
        compiler = shlex.split(compile_line.group(1))[0]
        if shutil.which(compiler) is None:
            self.skipTest(compiler + " not found")

        # out of the source tree, headers are included as "module/header.h"
        code = template.read_text().replace('"Colvar.h"', '"colvar/Colvar.h"')
        code = code.replace("TEMPLATE", "LOADTEST")

        current = os.getcwd()
        with tempfile.TemporaryDirectory() as temporary:
            os.chdir(temporary)
            try:
                Path("LoadTest.cpp").write_text(code)
                p = plumed.Plumed()
                p.cmd("setNatoms", 2)
                p.cmd("setLogFile", "test.log")
                p.cmd("init")
                p.cmd("readInputLine", "LOAD FILE=LoadTest.cpp")
                p.cmd("readInputLine", "d: LOADTEST ATOMS=1,2")
                p.cmd("readInputLine", "RESTRAINT ARG=d AT=0 KAPPA=1")

                box = array.array("d", [10, 0, 0, 0, 10, 0, 0, 0, 10])
                virial = array.array("d", [0] * 9)
                forces = array.array("d", [0] * 6)
                positions = array.array("d", [0, 0, 0, 1, 2, 3])
                p.cmd("setStep", 0)
                p.cmd("setBox", box)
                p.cmd("setPositions", positions)
                p.cmd("setForces", forces)
                p.cmd("setVirial", virial)
                p.cmd("calc")
                bias = array.array("d", [0])
                p.cmd("getBias", bias)
                p.finalize()
            finally:
                os.chdir(current)
        # same as DISTANCE: the distance squared over two
        self.assertAlmostEqual(bias[0], 7.0)
        self.assertAlmostEqual(forces[0], 1.0)


if __name__ == "__main__":
    unittest.main()
