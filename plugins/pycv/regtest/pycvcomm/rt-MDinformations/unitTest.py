import plumedCommunications as PLMD

# import plumedUtilities
log = open("pytest.log", "w")

def myPrint(*args,**kwargs): print(*args,**kwargs, file=log)

def myInit(action: PLMD.PythonCVInterface):
    myPrint(f"action label: {action.label}")
    return{"Value": PLMD.defaults.COMPONENT_NODEV,}

def myPrepare(action: PLMD.PythonCVInterface):
    myPrint(f"@step {action.getStep()}")
    myPrint(f"{action.getTime()=}")
    myPrint(f"{action.getTimeStep()=}")
    myPrint(f"{action.isRestart()=}")
    myPrint(f"{action.isExchangeStep()=}")
    return {}

    
def mypytest(action: PLMD.PythonCVInterface):
    return 0.0
