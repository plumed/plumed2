import plumedCommunications as PLMD
import numpy

log = open("pydist.log", "w")

print("Imported my pydist+.", file=log)

initForF={"Value": PLMD.defaults.COMPONENT_NODEV}

def function(action: PLMD.PythonFunction):
    arg = [action.argument(0),
           action.argument(1)
           ]
    return arg[0]*arg[0]*arg[1]