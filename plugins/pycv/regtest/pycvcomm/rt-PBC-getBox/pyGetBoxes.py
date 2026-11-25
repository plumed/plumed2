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

plumedInit = {
    "COMPONENTS": {
        "ax": plumedCommunications.defaults.COMPONENT,
        "ay": plumedCommunications.defaults.COMPONENT,
        "az": plumedCommunications.defaults.COMPONENT,
        "bx": plumedCommunications.defaults.COMPONENT,
        "by": plumedCommunications.defaults.COMPONENT,
        "bz": plumedCommunications.defaults.COMPONENT,
        "cx": plumedCommunications.defaults.COMPONENT,
        "cy": plumedCommunications.defaults.COMPONENT,
        "cz": plumedCommunications.defaults.COMPONENT,
    }
}


def pyBox(action: plumedCommunications.PythonCVInterface):
    d = action.getPbc().getBox()
    print(f"{d=}", file=log)
    ret = {}
    grad = {}
    for i, name in enumerate(["a", "b", "c"]):
        for j, coord in enumerate(["x", "y", "z"]):
            ret[f"{name}{coord}"] = (d[i, j], np.zeros((1, 3)))
    print(ret, file=log)
    return ret
