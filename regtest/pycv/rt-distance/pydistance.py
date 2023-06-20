# This is only a partial example. This function does not compute the
# gradient, so it is useless for biasing. See the other regression
# tests for how to auto-grad.

# And, of course, one should not call slow functions (such as print)
# in the CV calculation.

import numpy as np
import plumedCommunications
#import plumedUtilities
log=open("pydist.log","w")

print("Imported pydist.",file=log)

def pydist(x):
    return 0

def pydist_(action:plumedCommunications.PythonCVInterface):
    at0=action.getPosition(0)
    at1=action.getPosition(1)
    
    d = at0-at1
    return plumedUtilities.modulo(d)