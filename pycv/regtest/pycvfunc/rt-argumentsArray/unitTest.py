import plumedCommunications as PLMD
import numpy

log = open("pydist.log", "w")

print("Imported my pydist+.", file=log)

initForF = {"Value": PLMD.defaults.COMPONENT_NODEV}


def function(action: PLMD.PythonFunction):
    arg = action.arguments()

    return numpy.abs(numpy.sqrt(arg[0] ** 2 + arg[1] ** 2 + arg[2] ** 2) - arg[3])
