import plumedCommunications as PLMD
from sys import stderr as log
import numpy as np

print("Imported unitTest", file=log)


def myInit(action: PLMD.PythonCVInterface):
    t = np.loadtxt("massCharges.dat")
    action.data["masses"] = t[:, 1]
    action.data["charges"] = t[:, 2]
    print("Calling myInit", file=log)
    return {
        "COMPONENTS": {
            "mass0": PLMD.defaults.COMPONENT_NODEV,
            "mass1": PLMD.defaults.COMPONENT_NODEV,
            "charge0": PLMD.defaults.COMPONENT_NODEV,
            "charge1": PLMD.defaults.COMPONENT_NODEV,
        }
    }


def mypytest(action: PLMD.PythonCVInterface):
    masses = action.masses()
    charges = action.charges()
    action.data["masses"]
    action.data["charges"]
    ret = {}
    absoluteIndexes = action.absoluteIndexes()
    for i in range(action.nat):
        ret[f"mass{i}"] = (
            action.data["masses"][absoluteIndexes[i]] == masses[i]
        )
        ret[f"charge{i}"] = (
            action.data["charges"][absoluteIndexes[i]] == charges[i]
        )
    return ret
