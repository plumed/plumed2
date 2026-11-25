import plumedCommunications as PLMD

# import plumedUtilities
log = open("pytest.log", "w")

def myPrint(*args,**kwargs): print(*args,**kwargs, file=log)

def plumedInit(action: PLMD.PythonFunction):
    myPrint(f"action label: {action.label}")
    return{"Value": PLMD.defaults.COMPONENT_NODEV,}

def plumedCalculate(action: PLMD.PythonFunction):
    myPrint(f"@step {action.getStep()}")
    myPrint(f"{action.getTime()=}")
    myPrint(f"{action.getTimeStep()=}")
    myPrint(f"{action.isRestart()=}")
    myPrint(f"{action.isExchangeStep()=}")
    return 0.0
