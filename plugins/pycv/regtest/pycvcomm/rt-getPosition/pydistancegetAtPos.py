# This is only a partial example. This function does not compute the
# gradient, so it is useless for biasing. See the other regression
# tests for how to auto-grad.

# And, of course, one should not call slow functions (such as print)
# in the CV calculation.

import numpy as np
import plumedCommunications

# import plumedUtilities
log = open("pydist.log", "w")

print("Imported my pydist+.", file=log)

plumedInit = {"Value": plumedCommunications.defaults.COMPONENT_NODEV}


def pydist(action: plumedCommunications.PythonCVInterface):
    at: list = [
        action.getPosition(0),
        action.getPosition(1),
    ]
    d = at[0] - at[1]
    d = np.linalg.norm(d)
    print(f"{at[0]},{at[1]},{d}", file=log)

    # print(f"{at1},{at2}",file=log)
    return d
