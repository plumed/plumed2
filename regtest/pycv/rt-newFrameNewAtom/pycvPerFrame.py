# This is only a partial example. This function does not compute the
# gradient, so it is useless for biasing. See the other regression
# tests for how to auto-grad.

# And, of course, one should not call slow functions (such as print)
# in the CV calculation.

import numpy as np
import plumedCommunications
log=open("pycv.log","w")

print("Imported.",file=log)

def changeAtom(plmdAction:plumedCommunications.PythonCVInterface):
    print(f"pyCVCALLED")
    toret={
        "setAtomRequest":[0,int(plmdAction.getStep())+1]
    }
    #this is just for "fun"
    if plmdAction.getStep() == 3:
        toret["setAtomRequest"][1]=1
    print(toret)
    return toret
    
def pydist(action:plumedCommunications.PythonCVInterface):
    at:list[plumedCommunications.Vector3D]=[action.getPosition(0),action.getPosition(1)]
    
    d = plumedCommunications.modulo(at[0]-at[1])
    print(f"{at[0]},{at[1]},{d}",file=log)

    return d
