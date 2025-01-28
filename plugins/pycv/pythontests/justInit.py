import numpy as np
import plumedCommunications


def myInit(_: plumedCommunications.PythonCVInterface):
    return{"Value": plumedCommunications.defaults.COMPONENT_NODEV,}

def plumedCalculate(action: plumedCommunications.PythonCVInterface):
    at: np.ndarray = action.getPositions()
    return at[0][0]+at[0][1]
