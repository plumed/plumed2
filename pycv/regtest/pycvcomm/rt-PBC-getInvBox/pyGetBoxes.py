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


def pyInvBox(action: plumedCommunications.PythonCVInterface):
    invBox = action.getPbc().getInvBox()
    print(f"{invBox=}", file=log)
    ret = {}
    for i, name in enumerate(["a", "b", "c"]):
        for j, coord in enumerate(["x", "y", "z"]):
            ret[f"{name}i{coord}"] = (invBox[i, j], np.zeros((1, 3)))
    print(ret, file=log)
    return ret
