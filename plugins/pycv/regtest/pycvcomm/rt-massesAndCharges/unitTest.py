import plumedCommunications as PLMD

# import plumedUtilities
log = open("pydist.log", "w")

print("Imported my pydist+.", file=log)

plumedInit = {
    "COMPONENTS": {
        "mass0": PLMD.defaults.COMPONENT_NODEV,
        "mass1": PLMD.defaults.COMPONENT_NODEV,
        "charge0": PLMD.defaults.COMPONENT_NODEV,
        "charge1": PLMD.defaults.COMPONENT_NODEV,
    }
}


def mypytest(action: PLMD.PythonCVInterface):
    ret = {
        "mass0": action.mass(0),
        "charge0": action.charge(0),
        "mass1": action.mass(1),
        "charge1": action.charge(1),
    }

    return ret
