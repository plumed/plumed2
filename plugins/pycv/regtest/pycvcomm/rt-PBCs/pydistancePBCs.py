# This is only a partial example. This function does not compute the
# gradient, so it is useless for biasing. See the other regression
# tests for how to auto-grad.

# And, of course, one should not call slow functions (such as print)
# in the CV calculation.

import numpy as np
import plumedCommunications

# import plumedUtilities
log = open("pydist.log", "w")

print("Imported my pydist.", file=log)

plumedInit={"Value":plumedCommunications.defaults.COMPONENT_NODEV}

def pydist(action: plumedCommunications.PythonCVInterface):
    at: np.ndarray = action.getPositions()

    # print(at0,file=log)
    d = at[0] - at[1]
    d = action.getPbc().apply([d])
    assert d.shape[1] == 3, "d is not a (*,3) array"
    d = np.linalg.norm(d[0])

    # print(f"{at1},{at2}",file=log)
    return d
