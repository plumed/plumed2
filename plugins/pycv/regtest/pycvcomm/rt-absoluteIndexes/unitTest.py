import plumedCommunications as PLMD
from plumedCommunications.defaults import COMPONENT_NODEV as nodevnoperiod
# import plumedUtilities
log = open("pydist.log", "w")

print("Imported my pydist+.", file=log)
plumedInit = dict(
    COMPONENTS=dict(
        absoluteIndex0index=nodevnoperiod,
        absoluteIndex0serial=nodevnoperiod,
        absoluteIndex1index=nodevnoperiod,
        absoluteIndex1serial=nodevnoperiod,
    )
)


def mypytest(action: PLMD.PythonCVInterface):
    #this is the costly way of working:
    ret = {
        "absoluteIndex0index": action.absoluteIndexes()[0],
        "absoluteIndex0serial":action.absoluteIndexes()[0]+1,
    }
    #this should be faster
    indexes = action.absoluteIndexes()
    ret["absoluteIndex1index"] = indexes[1]
    ret["absoluteIndex1serial"] = indexes[1]+1
    return ret
