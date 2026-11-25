import plumedCommunications as PLMD
import numpy

log = open("pydist.log", "w")

print("Imported my pydist+.", file=log)

initForF = {
    "COMPONENTS": {
        "x": PLMD.defaults.COMPONENT_NODEV,
        "y": PLMD.defaults.COMPONENT_NODEV,
        "z": PLMD.defaults.COMPONENT_NODEV,
    }
}


def function(action: PLMD.PythonFunction):
    arg = action.arguments()

    return {
        "x": arg[0],
        "y": arg[1],
        "z": arg[2],
    }
