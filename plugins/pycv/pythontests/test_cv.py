import unittest
import numpy as np

from utilities_for_test import *

import os

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


class TestPyCV(unittest.TestCase):
    def test_nat(self):
        with cd(THIS_DIR):
            os.environ["PLUMED_MAXBACKUP"] = "0"
            traj, num_frames, num_atoms, box, virial, masses, forces, charges = (
                setUpTraj("traj.xyz")
            )
            plmd = preparePlumed(num_atoms)
            cvPy = create_plumed_var(plmd, "cvPy", "PYCVINTERFACE IMPORT=atoms_number")
            plmd.cmd("readInputLine", "PRINT FILE=colvar.out ARG=*")
            # Open an output file
            with open("logfile", "w+") as of:
                # Now analyze the trajectory
                for step in range(0, num_frames):
                    of.write("RUNNING ANALYSIS FOR STEP " + str(step) + "\n")
                    plmd.cmd("setStep", step)
                    plmd.cmd("setBox", box)
                    plmd.cmd("setMasses", masses)
                    plmd.cmd("setCharges", charges)
                    plmd.cmd("setPositions", traj[step])
                    plmd.cmd("setForces", forces)
                    plmd.cmd("setVirial", virial)
                    plmd.cmd("calc")

                    self.assertEqual(cvPy, 2)

    def test_atomPositions(self):
        with cd(THIS_DIR):
            os.environ["PLUMED_MAXBACKUP"] = "0"
            traj, num_frames, num_atoms, box, virial, masses, forces, charges = (
                setUpTraj("traj.xyz")
            )
            plmd = preparePlumed(num_atoms)

            cvPy = create_plumed_var(
                plmd,
                "cvPy",
                "PYCVINTERFACE ATOMS=1,4 IMPORT=pydistancegetAtPos CALCULATE=pydist",
            )
            cvCPP = create_plumed_var(plmd, "cvCPP", "DISTANCE ATOMS=1,4 NOPBC")

            plmd.cmd("readInputLine", "PRINT FILE=colvar.out ARG=*")
            # Open an output file
            with open("logfile", "w+") as of:
                # Now analyze the trajectory
                for step in range(0, num_frames):
                    of.write("RUNNING ANALYSIS FOR STEP " + str(step) + "\n")
                    plmd.cmd("setStep", step)
                    plmd.cmd("setBox", box)
                    plmd.cmd("setMasses", masses)
                    plmd.cmd("setCharges", charges)
                    plmd.cmd("setPositions", traj[step])
                    plmd.cmd("setForces", forces)
                    plmd.cmd("setVirial", virial)
                    plmd.cmd("calc")

                    np.testing.assert_almost_equal(cvPy, cvCPP, decimal=4)

    def test_atomPositionsPBC(self):
        with cd(THIS_DIR):
            os.environ["PLUMED_MAXBACKUP"] = "0"
            traj, num_frames, num_atoms, box, virial, masses, forces, charges = (
                setUpTraj("traj.xyz")
            )
            plmd = preparePlumed(num_atoms)

            cvPy = create_plumed_var(
                plmd,
                "cvPy",
                "PYCVINTERFACE ATOMS=1,4 IMPORT=pydistancePBCs CALCULATE=pydist",
            )
            cvCPP = create_plumed_var(plmd, "cvCPP", "DISTANCE ATOMS=1,4")

            plmd.cmd("readInputLine", "PRINT FILE=colvar.out ARG=*")
            # Open an output file
            with open("logfile", "w+") as of:
                # Now analyze the trajectory
                for step in range(0, num_frames):
                    of.write("RUNNING ANALYSIS FOR STEP " + str(step) + "\n")
                    plmd.cmd("setStep", step)
                    plmd.cmd("setBox", box)
                    plmd.cmd("setMasses", masses)
                    plmd.cmd("setCharges", charges)
                    plmd.cmd("setPositions", traj[step])
                    plmd.cmd("setForces", forces)
                    plmd.cmd("setVirial", virial)
                    plmd.cmd("calc")

                    np.testing.assert_almost_equal(cvPy, cvCPP, decimal=4)

    def test_newFrameNewAtomSTR(self):
        with cd(THIS_DIR):
            os.environ["PLUMED_MAXBACKUP"] = "0"
            traj, num_frames, num_atoms, box, virial, masses, forces, charges = (
                setUpTraj("trajnewFrameNewAtom.xyz")
            )
            plmd = preparePlumed(num_atoms)

            cvPy = create_plumed_var(
                plmd,
                "cvPy",
                "PYCVINTERFACE ATOMS=@mdatoms IMPORT=pycvPerFrameSTR CALCULATE=pydist PREPARE=changeAtom",
            )

            plmd.cmd("readInputLine", "PRINT FILE=colvar.out ARG=*")
            # Open an output file
            with open("logfile", "w+") as of:
                # Now analyze the trajectory
                for step in range(0, num_frames):
                    of.write("RUNNING ANALYSIS FOR STEP " + str(step) + "\n")
                    plmd.cmd("setStep", step)
                    plmd.cmd("setBox", box)
                    plmd.cmd("setMasses", masses)
                    plmd.cmd("setCharges", charges)
                    plmd.cmd("setPositions", traj[step])
                    plmd.cmd("setForces", forces)
                    plmd.cmd("setVirial", virial)
                    plmd.cmd("calc")

                    np.testing.assert_almost_equal(cvPy, step + 1.0, decimal=4)

    def test_newFrameNewAtom(self):
        with cd(THIS_DIR):
            os.environ["PLUMED_MAXBACKUP"] = "0"
            traj, num_frames, num_atoms, box, virial, masses, forces, charges = (
                setUpTraj("trajnewFrameNewAtom.xyz")
            )
            plmd = preparePlumed(num_atoms)

            cvPy = create_plumed_var(
                plmd,
                "cvPy",
                "PYCVINTERFACE ATOMS=@mdatoms IMPORT=pycvPerFrame CALCULATE=pydist PREPARE=changeAtom",
            )

            plmd.cmd("readInputLine", "PRINT FILE=colvar.out ARG=*")
            # Open an output file
            with open("logfile", "w+") as of:
                # Now analyze the trajectory
                for step in range(0, num_frames):
                    of.write("RUNNING ANALYSIS FOR STEP " + str(step) + "\n")
                    plmd.cmd("setStep", step)
                    plmd.cmd("setBox", box)
                    plmd.cmd("setMasses", masses)
                    plmd.cmd("setCharges", charges)
                    plmd.cmd("setPositions", traj[step])
                    plmd.cmd("setForces", forces)
                    plmd.cmd("setVirial", virial)
                    plmd.cmd("calc")

                    np.testing.assert_almost_equal(cvPy, step + 1.0, decimal=4)

    def test_loadAmodule_and_persistData(self):
        """This test loads a module that is a directory and stores some data within plumed"""
        with cd(THIS_DIR):
            os.environ["PLUMED_MAXBACKUP"] = "0"
            traj, num_frames, num_atoms, box, virial, masses, forces, charges = (
                setUpTraj("trajnewFrameNewAtom.xyz")
            )
            plmd = preparePlumed(num_atoms)

            cvPy = create_plumed_var(
                plmd,
                "cvPy",
                "PYCVINTERFACE ATOMS=@mdatoms IMPORT=pycvPersistentData CALCULATE=pydist INIT=pyinit",
            )

            plmd.cmd("readInputLine", "PRINT FILE=colvar.out ARG=*")

            # Now analyze the trajectory
            for step in range(0, num_frames):
                plmd.cmd("setStep", step)
                plmd.cmd("setBox", box)
                plmd.cmd("setMasses", masses)
                plmd.cmd("setCharges", charges)
                plmd.cmd("setPositions", traj[step])
                plmd.cmd("setForces", forces)
                plmd.cmd("setVirial", virial)
                plmd.cmd("calc")
                # this cv sums the number of the step till now:
                # to future me: "//" is integer (floor) division
                np.testing.assert_almost_equal(
                    cvPy, ((step) * (step + 1)) // 2, decimal=4
                )


if __name__ == "__main__":
    # Output to four decimal places only
    np.set_printoptions(precision=4)
    unittest.main()
