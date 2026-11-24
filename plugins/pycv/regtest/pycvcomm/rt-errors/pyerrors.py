import numpy as np
import plumedCommunications as PLMD

plumedInit = {"Value": PLMD.defaults.COMPONENT}


def ok(action: PLMD.PythonCVInterface):
    return action.getPosition(1)[0], np.zeros((action.nat, 3))


def raises(action):
    raise ValueError("raised in python")


def position(action: PLMD.PythonCVInterface):
    return action.getPosition(action.nat)[0], np.zeros((action.nat, 3))


def mass(action: PLMD.PythonCVInterface):
    return action.mass(-1), np.zeros((action.nat, 3))


def charge(action: PLMD.PythonCVInterface):
    return action.charge(action.nat), np.zeros((action.nat, 3))


def pbc(action: PLMD.PythonCVInterface):
    action.getPbc().apply(np.zeros((2, 2)))
    return 0.0, np.zeros((action.nat, 3))


def neighbourlist(action: PLMD.PythonCVInterface):
    return float(len(action.getNeighbourList())), np.zeros((action.nat, 3))


def argument(action: PLMD.PythonFunction):
    return action.argument(action.nargs), np.zeros(action.nargs)
