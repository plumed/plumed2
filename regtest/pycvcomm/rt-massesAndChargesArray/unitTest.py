import plumedCommunications as PLMD
from sys import stderr as log
import numpy as np

print("Imported unitTest", file=log)


def myInit(action: PLMD.PythonCVInterface):
    t = np.loadtxt("massCharges.dat")
    action.data["masses"] = t[:, 1]
    action.data["charges"] = t[:, 2]


def mypytest(action: PLMD.PythonCVInterface):
    masses = action.masses()
    charges = action.charges()
    action.data["masses"]
    action.data["charges"]
    ret = {}
        
    for i in range(action.nat):
        ret[f"mass{i}"] = (
            action.data["masses"][action.absoluteIndexes[i].index] == masses[i]
        )
        ret[f"charge{i}"] = (
            action.data["charges"][action.absoluteIndexes[i].index] == charges[i]
        )
    return ret
