# This is only a partial example. This function does not compute the
# gradient, so it is useless for biasing. See the other regression
# tests for how to auto-grad.

# And, of course, one should not call slow functions (such as print)
# in the CV calculation.

import numpy as np
import plumedCommunications

log = open("pycv.log", "w")

print("Imported.", file=log)
plumedInit={"Value": plumedCommunications.defaults.COMPONENT,}

def changeAtom(plmdAction: plumedCommunications.PythonCVInterface):
    print(f"pyCVCALLED")
    toret = {"setAtomRequest": [0, int(plmdAction.getStep()) + 1]}
    # this is just for "fun"
    if plmdAction.getStep() == 3:
        toret["setAtomRequest"][1] = 1
    print(toret)
    return toret


def pydist(action: plumedCommunications.PythonCVInterface):
    at: np.ndarray = action.getPositions()
    d = at[0] - at[1]
    d = np.linalg.norm(d)
    print(f"{at[0]},{at[1]},{d}", file=log)

    return d
