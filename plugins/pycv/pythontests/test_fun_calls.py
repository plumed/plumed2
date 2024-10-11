import unittest
import numpy as np

from utilities_for_test import *

import os

THIS_DIR = os.path.dirname(os.path.abspath(__file__))

# CALCULATE is tested "by default"


class TestFunctionCalls(unittest.TestCase):
    def test_INIT(self):
        with cd(THIS_DIR):
            os.environ["PLUMED_MAXBACKUP"] = "0"
            traj, num_frames, num_atoms, box, virial, masses, forces, charges = (
                setUpTraj("traj.xyz")
            )
            plmd = preparePlumed(num_atoms)
            d1 = create_plumed_var(plmd, "d1", "DISTANCE ATOMS=1,2")
            d2 = create_plumed_var(plmd, "d2", "DISTANCE ATOMS=1,3")

            fPy = create_plumed_var(
                plmd,
                "fPY",
                "PYFUNCTION IMPORT=pyfuncINITfunc INIT=initForF CALCULATE=function ARG=d1,d2",
            )

            plmd.cmd("readInputLine", "PRINT FILE=colvar.out ARG=*")

            for step in range(num_frames):
                plmd.cmd("setStep", step)
                plmd.cmd("setBox", box)
                plmd.cmd("setMasses", masses)
                plmd.cmd("setCharges", charges)
                plmd.cmd("setPositions", traj[step])
                plmd.cmd("setForces", forces)
                plmd.cmd("setVirial", virial)
                plmd.cmd("calc")

                np.testing.assert_almost_equal(fPy, d1 * d2, decimal=4)

    def test_INITDICT(self):
        with cd(THIS_DIR):
            os.environ["PLUMED_MAXBACKUP"] = "0"
            traj, num_frames, num_atoms, box, virial, masses, forces, charges = (
                setUpTraj("traj.xyz")
            )
            plmd = preparePlumed(num_atoms)

            plmd.cmd("readInputLine", "dc: DISTANCE ATOMS=1,2 COMPONENTS")
            shape = np.zeros(1, dtype=np.int_)
            plmd.cmd("getDataRank " + "dc.x", shape)
            plmd.cmd("getDataRank " + "dc.y", shape)
            plmd.cmd("getDataRank " + "dc.z", shape)
            dcx = np.zeros((1))
            dcy = np.zeros((1))
            dcz = np.zeros((1))
            plmd.cmd("setMemoryForData dc.x", dcx)
            plmd.cmd("setMemoryForData dc.y", dcy)
            plmd.cmd("setMemoryForData dc.z", dcz)

            d = create_plumed_var(plmd, "d", "DISTANCE ATOMS=1,2")

            fPy = create_plumed_var(
                plmd,
                "fPY",
                "PYFUNCTION IMPORT=pyfuncINITdict INIT=initForF CALCULATE=function ARG=dc.x,dc.y,dc.z,d",
            )

            plmd.cmd("readInputLine", "PRINT FILE=colvar.out ARG=*")

            for step in range(num_frames):
                plmd.cmd("setStep", step)
                plmd.cmd("setBox", box)
                plmd.cmd("setMasses", masses)
                plmd.cmd("setCharges", charges)
                plmd.cmd("setPositions", traj[step])
                plmd.cmd("setForces", forces)
                plmd.cmd("setVirial", virial)
                plmd.cmd("calc")

                np.testing.assert_almost_equal(
                    fPy, np.abs(np.sqrt(dcx**2 + dcy**2 + dcz**2) - d[0]), decimal=4
                )


if __name__ == "__main__":
    # Output to four decimal places only
    np.set_printoptions(precision=4)
    unittest.main()
