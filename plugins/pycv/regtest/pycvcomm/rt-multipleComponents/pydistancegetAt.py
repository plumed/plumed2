# This is only a partial example. This function does not compute the
# gradient, so it is useless for biasing. See the other regression
# tests for how to auto-grad.

# And, of course, one should not call slow functions (such as print)
# in the CV calculation.

import numpy as np
import plumedCommunications
from sys import stderr as log

# import plumedUtilities
# log = open("pydist.log", "w")

print("Imported my pydist.", file=log)
plumedInit = dict(
    COMPONENTS=dict(
        d01=plumedCommunications.defaults.COMPONENT,
        d02=plumedCommunications.defaults.COMPONENT,
    )
)


def pydist(action: plumedCommunications.PythonCVInterface):
    at: np.ndarray = action.getPositions()
    nat = at.shape[0]
    # print(at0,file=log)
    ret = {}
    for rn, j in [["d01", 1], ["d02", 2]]:
        d = at[0] - at[j]
        val = np.linalg.norm(d)
        grad = np.zeros((nat, 3))
        grad[0] = d / val
        grad[j] = -d / val
        # print(f"{rn} {val} {grad}", file=log)
        ret[rn] = (val, grad)
    return ret
