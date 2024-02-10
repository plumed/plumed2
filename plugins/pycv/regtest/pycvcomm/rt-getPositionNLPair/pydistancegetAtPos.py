# This is only a partial example. This function does not compute the
# gradient, so it is useless for biasing. See the other regression
# tests for how to auto-grad.

# And, of course, one should not call slow functions (such as print)
# in the CV calculation.

import numpy as np
import plumedCommunications
from sys import stderr as log

# import plumedUtilities
log = open("pydist.log", "w")

print("Imported my pydist+.", file=log)

plumedInit = {
    "COMPONENTS": {
        "x": plumedCommunications.defaults.COMPONENT,
        "y": plumedCommunications.defaults.COMPONENT,
        "z": plumedCommunications.defaults.COMPONENT,
    }
}


def pydistInPair(action: plumedCommunications.PythonCVInterface):
    # NB: This is not a realistic case of using the neigbour list!!!
    # cvPY: PYCVINTERFACE GROUPA=1,4 IMPORT=pydistancegetAtPos CALCULATE=pydist

    atoms = action.getPositions()
    nl = action.getNeighbourList()
    assert nl.size == 3
    assert len(nl) == nl.size

    x, y, z = [atoms[pair[0]] - atoms[pair[1]] for pair in nl.getClosePairs()]

    zero = np.zeros(atoms.shape)
    # MUST work with tuple or lists
    return {"x": [x[0], zero], "y": (y[1], zero), "z": (z[2], zero)}
