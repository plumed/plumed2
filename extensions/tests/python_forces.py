#!/usr/bin/env python3
"""Evaluate a fixed VES bias through Python against a selected installed kernel."""

import argparse
from pathlib import Path

import numpy as np
import plumed

parser = argparse.ArgumentParser()
parser.add_argument("--kernel", required=True)
parser.add_argument("--fixture", required=True, type=Path)
parser.add_argument("--output", required=True, type=Path)
parser.add_argument("--load-module", action="store_true")
args = parser.parse_args()
fixture = args.fixture.resolve()
lines = (fixture / "traj.gro").read_text().splitlines()
natoms = int(lines[1])
positions = np.array(
    [
        [float(line[start : start + 8]) for start in (20, 28, 36)]
        for line in lines[2 : 2 + natoms]
    ],
    dtype=np.float64,
)
forces = np.zeros_like(positions)
virial = np.zeros((3, 3), dtype=np.float64)
box = np.diag([float(v) for v in lines[2 + natoms].split()[:3]])
bias = np.zeros(1, dtype=np.float64)
with plumed.Plumed(kernel=str(Path(args.kernel).resolve())) as p:
    p.cmd("setNatoms", natoms)
    p.cmd("setLogFile", "python-ves.log")
    p.cmd("setTimestep", 0.002)
    p.cmd("init")
    if args.load_module:
        p.cmd("readInputLine", "LOAD MODULE=ves")
    p.cmd("readInputLine", "phi: TORSION ATOMS=5,7,9,15 NOPBC")
    p.cmd("readInputLine", "bf1: BF_FOURIER ORDER=6 MINIMUM=-pi MAXIMUM=pi")
    p.cmd(
        "readInputLine",
        "ves1: VES_LINEAR_EXPANSION ARG=phi BASIS_FUNCTIONS=bf1 "
        "TEMP=300 GRID_BINS=100 COEFFS=" + str(fixture / "coeffs-in.data"),
    )
    # Live metadata must expose the actions visible to this exact object.
    builder = plumed.InputBuilder(toplumed=p)
    assert hasattr(builder, "BF_FOURIER")
    p.cmd("setStep", 0)
    p.cmd("setBox", box)
    p.cmd("setPositions", positions)
    p.cmd("setMasses", np.ones(natoms, dtype=np.float64))
    p.cmd("setCharges", np.zeros(natoms, dtype=np.float64))
    p.cmd("setForces", forces)
    p.cmd("setVirial", virial)
    p.cmd("calc")
    p.cmd("getBias", bias)
np.savez(args.output, bias=bias, forces=forces, virial=virial)
print("VES bias:", bias[0], "force norm:", np.linalg.norm(forces))
