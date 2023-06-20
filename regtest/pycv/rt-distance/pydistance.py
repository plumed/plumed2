# This is only a partial example. This function does not compute the
# gradient, so it is useless for biasing. See the other regression
# tests for how to auto-grad.

# And, of course, one should not call slow functions (such as print)
# in the CV calculation.

import numpy as np
import plumedCommunications
#import plumedUtilities
log=open("pydist.log","w")

print("Imported my pydist.",file=log)

def pydist_(x):
    #print("call",file=log)
    return 0

def pydist(action:plumedCommunications.PythonCVInterface):
    print("call",file=log)
    at0:plumedCommunications.Vector3D=action.getPosition(0)
    #at0=action.getPosition(0)
    at1:plumedCommunications.Vector3D=action.getPosition(1)
    
    #print(at0,file=log)
    d = plumedCommunications.modulo(at0-at1)
    print(f"{at0},{at1},{d}",file=log)
    
    #print(f"{at1},{at2}",file=log)
    return d
    