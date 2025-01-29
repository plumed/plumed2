from plumed import Plumed
from contextlib import contextmanager
import os
import numpy as np

__all__ = ["read_xyz", "setUpTraj", "preparePlumed", "cd", "create_plumed_var"]


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


def setUpTraj(trajfile: str):
    traj = read_xyz(trajfile)
    num_frames = len(traj)
    num_atoms = traj[0].shape[0]

    # Create arrays for stuff
    box = np.diag(12.41642 * np.ones(3, dtype=np.float64))
    virial = np.zeros((3, 3), dtype=np.float64)
    masses = np.ones(num_atoms, dtype=np.float64)
    forces = np.random.rand(num_atoms, 3)
    charges = np.zeros(num_atoms, dtype=np.float64)

    return traj, num_frames, num_atoms, box, virial, masses, forces, charges


def preparePlumed(num_atoms: int):
    from pycv import getPythonCVInterface

    # Create PLUMED object and read input
    plmd = Plumed()

    # not really needed, used to check https://github.com/plumed/plumed2/issues/916
    plumed_version = np.zeros(1, dtype=np.intc)
    plmd.cmd("getApiVersion", plumed_version)
    plmd.cmd("setMDEngine", "python")
    plmd.cmd("setTimestep", 1.0)
    plmd.cmd("setKbT", 1.0)
    plmd.cmd("setNatoms", num_atoms)
    plmd.cmd("setLogFile", "test.log")
    plmd.cmd("init")
    plmd.cmd("readInputLine", f"LOAD FILE={getPythonCVInterface()}")
    return plmd


def create_plumed_var(plmd: Plumed, name: str, command: str):
    plmd.cmd("readInputLine", name + ": " + command)
    shape = np.zeros(1, dtype=np.int_)
    plmd.cmd("getDataRank " + name, shape)
    data = np.zeros((1))
    plmd.cmd("setMemoryForData " + name, data)
    return data
