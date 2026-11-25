import unittest
import numpy as np

from utilities_for_test import *

import os

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


class TestPyCVCalls(unittest.TestCase):
    def test_INIT(self):
        with cd(THIS_DIR):
            os.environ["PLUMED_MAXBACKUP"] = "0"
            traj, _, num_atoms, box, virial, masses, forces, charges = setUpTraj(
                "traj.xyz"
            )
            plmd = preparePlumed(num_atoms)

            cvPy = create_plumed_var(
                plmd, "cvPy", "PYCVINTERFACE ATOMS=4 IMPORT=justInit INIT=myInit"
            )

            plmd.cmd("readInputLine", "PRINT FILE=colvar.out ARG=*")

            step = 0

            plmd.cmd("setStep", step)
            plmd.cmd("setBox", box)
            plmd.cmd("setMasses", masses)
            plmd.cmd("setCharges", charges)
            plmd.cmd("setPositions", traj[step])
            plmd.cmd("setForces", forces)
            plmd.cmd("setVirial", virial)
            plmd.cmd("calc")
            np.testing.assert_almost_equal(cvPy, 10.0, decimal=4)

    def test_INITDICT(self):
        with cd(THIS_DIR):
            os.environ["PLUMED_MAXBACKUP"] = "0"
            traj, _, num_atoms, box, virial, masses, forces, charges = setUpTraj(
                "traj.xyz"
            )
            plmd = preparePlumed(num_atoms)

            cvPy = create_plumed_var(
                plmd, "cvPy", "PYCVINTERFACE ATOMS=4 IMPORT=justInitDict"
            )

            plmd.cmd("readInputLine", "PRINT FILE=colvar.out ARG=*")

            step = 0

            plmd.cmd("setStep", step)
            plmd.cmd("setBox", box)
            plmd.cmd("setMasses", masses)
            plmd.cmd("setCharges", charges)
            plmd.cmd("setPositions", traj[step])
            plmd.cmd("setForces", forces)
            plmd.cmd("setVirial", virial)
            plmd.cmd("calc")
            np.testing.assert_almost_equal(cvPy, 10.0, decimal=4)

    def test_UPDATE(self):
        with cd(THIS_DIR):
            os.environ["PLUMED_MAXBACKUP"] = "0"
            traj, _, num_atoms, box, virial, masses, forces, charges = setUpTraj(
                "traj.xyz"
            )
            plmd = preparePlumed(num_atoms)

            cvPy = create_plumed_var(
                plmd, "cvPy", "PYCVINTERFACE ATOMS=4 IMPORT=justUpdate UPDATE=myUpdate"
            )

            plmd.cmd("readInputLine", "PRINT FILE=colvar.out ARG=*")

            step = 0

            plmd.cmd("setStep", step)
            plmd.cmd("setBox", box)
            plmd.cmd("setMasses", masses)
            plmd.cmd("setCharges", charges)
            plmd.cmd("setPositions", traj[step])
            plmd.cmd("setForces", forces)
            plmd.cmd("setVirial", virial)
            plmd.cmd("calc")
            np.testing.assert_almost_equal(cvPy, 10.0, decimal=4)

    def test_PREPARE(self):
        with cd(THIS_DIR):
            os.environ["PLUMED_MAXBACKUP"] = "0"
            traj, _, num_atoms, box, virial, masses, forces, charges = setUpTraj(
                "traj.xyz"
            )
            plmd = preparePlumed(num_atoms)
            # atoms=4 but the module choses 1
            cvPy = create_plumed_var(
                plmd,
                "cvPy",
                "PYCVINTERFACE ATOMS=4 IMPORT=justPrepare PREPARE=plumedPrepare",
            )
            plmd.cmd("readInputLine", "PRINT FILE=colvar.out ARG=*")
            # Open an output file

            step = 0

            plmd.cmd("setStep", step)
            plmd.cmd("setBox", box)
            plmd.cmd("setMasses", masses)
            plmd.cmd("setCharges", charges)
            plmd.cmd("setPositions", traj[step])
            plmd.cmd("setForces", forces)
            plmd.cmd("setVirial", virial)
            plmd.cmd("calc")

            np.testing.assert_almost_equal(cvPy, 5.0, decimal=4)


if __name__ == "__main__":
    # Output to four decimal places only
    np.set_printoptions(precision=4)
    unittest.main()
