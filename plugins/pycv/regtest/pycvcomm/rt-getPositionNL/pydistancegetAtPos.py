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

plumedInit = {"Value": plumedCommunications.defaults.COMPONENT_NODEV}


def pydist(action: plumedCommunications.PythonCVInterface):
    # NB: This is not a realistic case of using the neigbour list!!!
    # cvPY: PYCVINTERFACE GROUPA=1,4 IMPORT=pydistancegetAtPos CALCULATE=pydist
    # ^ using this line should behave like calling "ATOMS=1,4" (with this function)
    atoms = action.getPositions()
    nl = action.getNeighbourList()
    assert nl.size == 1
    assert len(nl) == nl.size
    NLlist = nl.getClosePairs()[0]

    d = atoms[NLlist[0]] - atoms[NLlist[1]]
    d = np.linalg.norm(d)
    print(f"{atoms[0]},{atoms[1]},{d}", file=log)

    return d
