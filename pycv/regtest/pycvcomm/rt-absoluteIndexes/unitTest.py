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
