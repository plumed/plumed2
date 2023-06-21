import numpy as np
import plumedCommunications


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
    d = atoms[0] * atoms[1][0]
    return d[0]


def pymulfv(action: plumedCommunications.PythonCVInterface):
    atoms = getAtoms01(action)
    d = atoms[1][0] * atoms[0]
    return d[0]


def pyneg(action: plumedCommunications.PythonCVInterface):
    atom = action.getPosition(0)
    d = -atom
    return d[0]


log = open("pyArray.log", "w")

print("Checking if vector3 is correclty converted to numpy", file=log)


def pyArray(action: plumedCommunications.PythonCVInterface):
    # The number are formatted to ensure reproducibility
    atom = action.getPosition(0)
    print(f"{atom[0]:.3f} {atom[1]:.3f} {atom[2]:.3f}", file=log)
    d = atom.toArray()
    print(
        f"{d[0]:.3f} {d[1]:.3f} {d[2]:.3f} {d.shape} {isinstance(d, np.ndarray)}",
        file=log,
    )
    return 0
