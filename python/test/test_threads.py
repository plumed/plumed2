# Plumed objects used at the same time by several threads: cmd releases the
# GIL while PLUMED runs, and its arguments must stay valid meanwhile (the key of
# a command with a string value was freed before the call)
import os
import tempfile
import threading
import unittest

import plumed

try:
    import numpy as np

    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False


@unittest.skipUnless(HAS_NUMPY, "requires numpy")
class TestThreads(unittest.TestCase):
    def test_concurrent_objects(self):
        errors = []
        with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as logdir:

            def work(name):
                try:
                    for i in range(50):
                        p = plumed.Plumed()
                        p.cmd("setMDEngine", "python")
                        p.cmd("setNatoms", 2)
                        p.cmd("setLogFile", os.path.join(logdir, f"{name}.log"))
                        p.cmd("init")
                        p.cmd("readInputLine", "d: DISTANCE ATOMS=1,2")
                        colvar = os.path.join(logdir, f"{name}.colvar")
                        p.cmd("readInputLine", f"PRINT ARG=d FILE={colvar}")
                        for step in range(25):
                            x = 1.0 + 0.01 * step
                            p.cmd("setStep", step)
                            p.cmd("setBox", np.diag(10.0 * np.ones(3)))
                            p.cmd("setMasses", np.ones(2))
                            p.cmd("setCharges", np.zeros(2))
                            p.cmd("setPositions", np.array([[0.0, 0, 0], [x, 0, 0]]))
                            p.cmd("setForces", np.zeros((2, 3)))
                            p.cmd("setVirial", np.zeros((3, 3)))
                            p.cmd("calc")
                        p.finalize()
                except Exception as e:
                    errors.append(f"{name}: {e}")

            threads = [threading.Thread(target=work, args=(f"t{i}",)) for i in range(4)]
            for t in threads:
                t.start()
            for t in threads:
                t.join()
        self.assertEqual(errors, [])


if __name__ == "__main__":
    unittest.main()
