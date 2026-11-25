import plumedCommunications as PLMD
import numpy

# import plumedUtilities
log = open("pydist.log", "w")

print("Imported my pydist+.", file=log)
plumedInit = {
    "COMPONENTS": {
        "dbefore": PLMD.defaults.COMPONENT_NODEV,
        "dafter": PLMD.defaults.COMPONENT_NODEV,
    }
}


def mypytest(action: PLMD.PythonCVInterface):
    atoms = action.getPositions()
    d = atoms[0] - atoms[1]
    ret = {
        "dbefore": numpy.linalg.norm(d),
    }
    action.makeWhole()
    atoms = action.getPositions()
    d = atoms[0] - atoms[1]
    ret["dafter"] = numpy.linalg.norm(d)
    return ret
