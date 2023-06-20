# This is only a partial example. This function does not compute the
# gradient, so it is useless for biasing. See the other regression
# tests for how to auto-grad.

# And, of course, one should not call slow functions (such as print)
# in the CV calculation.

import numpy as np
import plumedCommunications

# import plumedUtilities


def getAtoms01(
    action: plumedCommunications.PythonCVInterface,
) -> "list[plumedCommunications.Vector3D]":
    return [
        action.getPosition(0),
        action.getPosition(1),
    ]


def pymod(action: plumedCommunications.PythonCVInterface):
    atoms = getAtoms01(action)
    d = plumedCommunications.modulo(atoms[0] - atoms[1])
    return d

def pymod2(action: plumedCommunications.PythonCVInterface):
    atoms = getAtoms01(action)
    d = plumedCommunications.modulo2(atoms[0] - atoms[1])
    return d

def pymodOBJ(action: plumedCommunications.PythonCVInterface):
    atoms = getAtoms01(action)
    d = atoms[0] - atoms[1]
    return d.modulo()

def pymod2OBJ(action: plumedCommunications.PythonCVInterface):
    atoms = getAtoms01(action)
    d = atoms[0] - atoms[1]
    return d.modulo2()

def pysum(action: plumedCommunications.PythonCVInterface):
    atoms = getAtoms01(action)
    d = atoms[0] + atoms[1]
    return d[0]


def pysub(action: plumedCommunications.PythonCVInterface):
    atoms = getAtoms01(action)
    d = atoms[0] - atoms[1]
    return d[0]

def pymul(action: plumedCommunications.PythonCVInterface):
    atoms = getAtoms01(action)
    d = atoms[0]*atoms[1][0]
    return d[0]

def pymulfv(action: plumedCommunications.PythonCVInterface):
    atoms = getAtoms01(action)
    d = atoms[1][0]*atoms[0]
    return d[0]

def pyneg(action: plumedCommunications.PythonCVInterface):
    atom = action.getPosition(0)
    d = -atom
    return d[0]