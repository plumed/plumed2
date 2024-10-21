import plumedCommunications as PLMD
import numpy

log = open("pydist.log", "w")

print("Imported my pydist+.", file=log)

def initForF(_: PLMD.PythonFunction): return {"Value": PLMD.defaults.COMPONENT}

def function(action: PLMD.PythonFunction):
    arg = [action.argument(0),
           action.argument(1)
           ]
    return arg[0]*arg[1], [arg[1], arg[0]]
