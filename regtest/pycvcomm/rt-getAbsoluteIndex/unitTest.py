# This is only a partial example. This function does not compute the
# gradient, so it is useless for biasing. See the other regression
# tests for how to auto-grad.

# And, of course, one should not call slow functions (such as print)
# in the CV calculation.

import numpy as np
import plumedCommunications as PLMD

# import plumedUtilities
log = open("pydist.log", "w")

print("Imported my pydist+.", file=log)


def mypytest(action: PLMD.PythonCVInterface):
    
    print(f"{action.getAbsoluteIndex(0).index=}, {action.getAbsoluteIndex(0).serial=}", file=log)

    return {"at0":action.getAbsoluteIndex(0).index,"at1":action.getAbsoluteIndex(0).serial}
