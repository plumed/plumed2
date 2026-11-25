# This is only a partial example. This function does not compute the
# gradient, so it is useless for biasing. See the other regression
# tests for how to auto-grad.

# And, of course, one should not call slow functions (such as print)
# in the CV calculation.

import numpy as np
import plumedCommunications
from plumedCommunications.defaults import COMPONENT_NODEV
from sys import stderr as log
#log = open("pycv.log", "w")

print("Imported pycvPersistentData", file=log)

def pyinit(plmdAction: plumedCommunications.PythonCVInterface):
    print("Calling pyinit", file=log)
    plmdAction.log("---Calling pyinit---")
    plmdAction.lognl("Logging from Python :)")
    print(f"{plmdAction.data=}", file=log)
    plmdAction.data["pycv"]=0
    print(f"{plmdAction.data=}", file=log)
    return {"Value":COMPONENT_NODEV}

def pydist(plmdAction: plumedCommunications.PythonCVInterface):
    plmdAction.log("Calling pydist: ")
    plmdAction.lognl(plmdAction.getStep())
    print("Calling pydist", file=log)
    print(f"{plmdAction.data=}, {plmdAction.getStep()=}", file=log)
    
    plmdAction.data["pycv"]+=plmdAction.getStep()
    d=plmdAction.data["pycv"]
    return d
