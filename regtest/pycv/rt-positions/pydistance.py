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
    at:list[plumedCommunications.Vector3D]=action.getPositions()
    
    #print(at0,file=log)
    d = plumedCommunications.modulo(at[0]-at[1])
    print(f"{at[0]},{at[1]},{d}",file=log)
    
    #print(f"{at1},{at2}",file=log)
    return d

def pyX(action:plumedCommunications.PythonCVInterface):
    #this tests that this should work also with only one atom in the request
    at:list[plumedCommunications.Vector3D]=action.getPositions()
    
    #print(at0,file=log)
    d = at[0][0]

    
    #print(f"{at1},{at2}",file=log)
    return d

def pyCount(action:plumedCommunications.PythonCVInterface):
    #this tests that this should work also with only one atom in the request
    at:list[plumedCommunications.Vector3D]=action.getPositions()
    #print(at0,file=log)
    return len(at)