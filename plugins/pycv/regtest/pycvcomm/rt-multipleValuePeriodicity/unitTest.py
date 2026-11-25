import plumedCommunications as PLMD
import numpy

log = open("pydist.log", "w")

print("Imported my pydist+.", file=log)

init = dict(
    COMPONENTS=dict(
        nonPeriodic=dict(period=None),
        Periodic={"period": ["0", "1.3"]},
        PeriodicPI={"period": ["0", "pi"]},
    ),
    ATOMS="1,2",
)


def mypytest(action: PLMD.PythonCVInterface):
    ret = {
        "nonPeriodic": action.getStep(),
        "Periodic": action.getStep(),
        "PeriodicPI": action.getStep(),
    }

    return ret
