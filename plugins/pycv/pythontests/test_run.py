import unittest
import numpy as np
from plumed import Plumed

import os
from contextlib import contextmanager

THIS_DIR = os.path.dirname(os.path.abspath(__file__))


@contextmanager
def cd(newdir):
    prevdir = os.getcwd()
    os.chdir(newdir)
    try:
        yield
    finally:
        os.chdir(prevdir)


def read_xyz(filename: str):
    xyz = open(filename)
    n_atoms = int(xyz.readline())
    _ = xyz.readline()
    trajectory = []
    while True:
        atom_type = np.zeros(n_atoms).astype(str)
        coordinates = np.zeros([n_atoms, 3])
        for i in range(0, n_atoms):
            line = xyz.readline()
            atom, x, y, z = line.split()
            atom_type[i] = atom
            coordinates[i, :] = np.array([x, y, z], dtype=np.float64)
        trajectory.append(coordinates)
        nextline = xyz.readline()
        if nextline == "":
            break
        c_atoms = int(nextline)
        if c_atoms != n_atoms:
            break
        _ = xyz.readline()
    xyz.close()
    return trajectory


def create_plumed_var(plmd: Plumed, name: str, command: str):
    plmd.cmd("readInputLine", name + ": " + command)
    shape = np.zeros(1, dtype=np.int_)
    plmd.cmd("getDataRank " + name, shape)
    data = np.zeros((1))
    plmd.cmd("setMemoryForData " + name, data)
    return data


class TestPyCV(unittest.TestCase):
    def setUpTraj(self):
        self.traj = read_xyz("traj.xyz")
        self.num_frames = len(self.traj)
        self.num_atoms = self.traj[0].shape[0]

        # Create arrays for stuff
        self.box = np.diag(12.41642 * np.ones(3, dtype=np.float64))
        self.virial = np.zeros((3, 3), dtype=np.float64)
        self.masses = np.ones(self.num_atoms, dtype=np.float64)
        self.forces = np.random.rand(self.num_atoms, 3)
        self.charges = np.zeros(self.num_atoms, dtype=np.float64)

    def preparePlumed(self):
        from pycv import getPythonCVInterface

        # Create PLUMED object and read input
        plmd = Plumed()

        # not really needed, used to check https://github.com/plumed/plumed2/issues/916
        plumed_version = np.zeros(1, dtype=np.intc)
        plmd.cmd("getApiVersion", plumed_version)
        plmd.cmd("setMDEngine", "python")
        plmd.cmd("setTimestep", 1.0)
        plmd.cmd("setKbT", 1.0)
        plmd.cmd("setNatoms", self.num_atoms)
        plmd.cmd("setLogFile", "test.log")
        plmd.cmd("init")
        plmd.cmd("readInputLine", f"LOAD FILE={getPythonCVInterface()}")
        return plmd

    def test_nat(self):
        with cd(THIS_DIR):
            os.environ["PLUMED_MAXBACKUP"] = "0"
            self.setUpTraj()
            plmd = self.preparePlumed()
            cvPy = create_plumed_var(plmd, "cvPy", "PYCVINTERFACE IMPORT=atoms_number")
            plmd.cmd("readInputLine", "PRINT FILE=colvar.out ARG=*")
            # Open an output file
            with open("logfile", "w+") as of:
                # Now analyze the trajectory
                for step in range(0, self.num_frames):
                    of.write("RUNNING ANALYSIS FOR STEP " + str(step) + "\n")
                    plmd.cmd("setStep", step)
                    plmd.cmd("setBox", self.box)
                    plmd.cmd("setMasses", self.masses)
                    plmd.cmd("setCharges", self.charges)
                    plmd.cmd("setPositions", self.traj[step])
                    plmd.cmd("setForces", self.forces)
                    plmd.cmd("setVirial", self.virial)
                    plmd.cmd("calc")

                    self.assertEqual(cvPy, 2)

    def test_atomPositions(self):
        with cd(THIS_DIR):
            os.environ["PLUMED_MAXBACKUP"] = "0"
            self.setUpTraj()
            plmd = self.preparePlumed()
            
            cvPy=create_plumed_var(plmd, "cvPy",  "PYCVINTERFACE ATOMS=1,4 IMPORT=pydistancegetAtPos CALCULATE=pydist")
            cvCPP=create_plumed_var(plmd, "cvCPP", "DISTANCE ATOMS=1,4 NOPBC")

            plmd.cmd("readInputLine", "PRINT FILE=colvar.out ARG=*")
            # Open an output file
            with open("logfile", "w+") as of:
                # Now analyze the trajectory
                for step in range(0, self.num_frames):
                    of.write("RUNNING ANALYSIS FOR STEP " + str(step) + "\n")
                    plmd.cmd("setStep", step)
                    plmd.cmd("setBox", self.box)
                    plmd.cmd("setMasses", self.masses)
                    plmd.cmd("setCharges", self.charges)
                    plmd.cmd("setPositions", self.traj[step])
                    plmd.cmd("setForces", self.forces)
                    plmd.cmd("setVirial", self.virial)
                    plmd.cmd("calc")

                    np.testing.assert_almost_equal(cvPy, cvCPP,decimal=4)

    def test_atomPositionsPBC(self):
        with cd(THIS_DIR):
            os.environ["PLUMED_MAXBACKUP"] = "0"
            self.setUpTraj()
            plmd = self.preparePlumed()
            
            cvPy=create_plumed_var(plmd, "cvPy",  "PYCVINTERFACE ATOMS=1,4 IMPORT=pydistancePBCs CALCULATE=pydist")
            cvCPP=create_plumed_var(plmd, "cvCPP", "DISTANCE ATOMS=1,4")

            plmd.cmd("readInputLine", "PRINT FILE=colvar.out ARG=*")
            # Open an output file
            with open("logfile", "w+") as of:
                # Now analyze the trajectory
                for step in range(0, self.num_frames):
                    of.write("RUNNING ANALYSIS FOR STEP " + str(step) + "\n")
                    plmd.cmd("setStep", step)
                    plmd.cmd("setBox", self.box)
                    plmd.cmd("setMasses", self.masses)
                    plmd.cmd("setCharges", self.charges)
                    plmd.cmd("setPositions", self.traj[step])
                    plmd.cmd("setForces", self.forces)
                    plmd.cmd("setVirial", self.virial)
                    plmd.cmd("calc")

                    np.testing.assert_almost_equal(cvPy, cvCPP,decimal=4)            


if __name__ == "__main__":
    # Output to four decimal places only
    np.set_printoptions(precision=4)
    unittest.main()
