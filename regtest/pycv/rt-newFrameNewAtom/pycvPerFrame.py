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
    if plmdAction.getStep() == 3:
        toret["setAtomRequest"][1]=1
    print(toret)
    return toret
    

def cv1(x):
    print(x,file=log)           # But don't
    d = x[0,:]-x[1,:]
    # If computing a gradient, return it as a 2nd return value.
    return np.sqrt(np.dot(d,d))
