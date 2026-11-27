import os
import unittest

import calculation_context
from utilities_for_test import cd, create_plumed_var, preparePlumed, setUpTraj

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


class TestCalculationContext(unittest.TestCase):
    def test_cv_and_function_context(self):
        with cd(THIS_DIR):
            os.environ["PLUMED_MAXBACKUP"] = "0"
            calculation_context.counts.clear()
            traj, _, natoms, box, virial, masses, forces, charges = setUpTraj(
                "traj.xyz"
            )
            p = preparePlumed(natoms)
            cv = create_plumed_var(p, "cv", "PYCVINTERFACE IMPORT=calculation_context")
            fun = create_plumed_var(
                p,
                "fun",
                "PYFUNCTION ARG=cv IMPORT=calculation_context INIT=functionInit",
            )
            p.cmd("readInputLine", "PRINT ARG=cv,fun FILE=context.out")
            for step, exchange, auxiliary, expected in (
                (0, False, False, 100),
                (1, False, True, 110),
                (1, True, True, 111),
                (2, False, False, 200),
            ):
                p.cmd("setStep", step)
                p.cmd("setExchangeStep", int(exchange))
                p.cmd("setBox", box)
                p.cmd("setMasses", masses)
                p.cmd("setCharges", charges)
                p.cmd("setPositions", traj[0])
                p.cmd("setForces", forces)
                p.cmd("setVirial", virial)
                p.cmd("prepareCalc")
                p.cmd("performCalcNoForces" if auxiliary else "performCalc")
                self.assertEqual(float(cv[0]), expected)
                self.assertEqual(float(fun[0]), expected)
            p.finalize()
