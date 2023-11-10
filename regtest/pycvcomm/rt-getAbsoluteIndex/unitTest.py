# This is only a partial example. This function does not compute the
# gradient, so it is useless for biasing. See the other regression
# tests for how to auto-grad.

# And, of course, one should not call slow functions (such as print)
# in the CV calculation.

import plumedCommunications as PLMD
# import plumedUtilities
log = open("pydist.log", "w")

print("Imported my pydist+.", file=log)

def mypytest(action: PLMD.PythonCVInterface):
    ret={"absoluteIndex0index":action.absoluteIndexes[0].index,
        "absoluteIndex0serial":action.absoluteIndexes[0].serial}
    indexes = action.absoluteIndexes
    ret["absoluteIndex1index"]=indexes[1].index
    ret["absoluteIndex1serial"]=indexes[1].serial
    #the following lines are guarantee to fail :)
    #action.absoluteIndexes[0].index=0
    #indexes[1].index=0
    return ret
