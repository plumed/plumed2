import numpy as np
import plumedCommunications as PLMD

# With a thread-MPI build of PLUMED the replicas are threads of the same process:
# they import this module once, and share its global variables. The state of an
# action is in its own data dictionary.

plumedInit = {"Value": PLMD.defaults.COMPONENT_NODEV}
initFunction = {"Value": PLMD.defaults.COMPONENT_NODEV}


def distance(action: PLMD.PythonCVInterface):
    calls = action.data.get("calls", 0)
    if calls != action.getStep():
        raise RuntimeError(f"call {calls} at step {action.getStep()}")
    action.data["calls"] = calls + 1
    at = action.getPositions()
    return np.linalg.norm(at[0] - at[1])


def square(action: PLMD.PythonFunction):
    return action.argument(0) ** 2
