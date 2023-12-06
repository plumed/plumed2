import plumedCommunications as PLMD
import numpy

log = open("pydist.log", "w")

print("Imported my pydist+.", file=log)


def init(action: PLMD.PythonCVInterface):
    return {"COMPONENTS": {"val":{"period": ["0", 1.3]}}}

initForF={"Value": PLMD.defaults.COMPONENT_NODEV}


def mypytest(action: PLMD.PythonCVInterface):
    ret = [action.getStep()]

    return ret

def function(action: PLMD.PythonFunction):
    arg = [action.argument(0)]
    return arg[0]*arg[0]