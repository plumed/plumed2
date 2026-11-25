import plumedCommunications as PLMD

from sys import stderr as log
# log = open("pydist.log", "w")

print("Imported my pydist+.", file=log)


def plumedInitio(action: PLMD.PythonCVInterface):
    # return {"Value": {"period": [None,0,0]}}
    # return {"Value": {"period": None}}
    return {"Value": {"period": ["0",0.3]}}
plumedInit={"Value": {"period": None},"NOPBC":False,
            "ATOMS":"1,2"}

def plumedCalculate(action: PLMD.PythonCVInterface):
    ret = [action.nat]
    print ("Hello from plumedCalculate")
    return ret
