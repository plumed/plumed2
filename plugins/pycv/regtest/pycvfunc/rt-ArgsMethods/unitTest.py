import plumedCommunications as PLMD
import numpy

log = open("pydist.log", "w")

print("Imported my pydist+.", file=log)

initForF = {
    "COMPONENTS": {
        "difference": PLMD.defaults.COMPONENT_NODEV,
        "bringBackInPbc": PLMD.defaults.COMPONENT_NODEV,
    }
}


def function(action: PLMD.PythonFunction):
    arg = action.arguments()

    return {
        "difference": action.difference(0, action.getStep(), arg[0]),
        "bringBackInPbc": action.bringBackInPbc(0, arg[0] * action.getStep()),
    }
