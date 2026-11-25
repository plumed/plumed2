# This is only a partial example. This function does not compute the
# gradient, so it is useless for biasing. See the other regression
# tests for how to auto-grad.

# And, of course, one should not call slow functions (such as print)
# in the CV calculation.

import plumedCommunications
import pydoc

def plumedInit(_):
    #rather than the other doc example this creates a very ugly html manua;
    pydoc.writedoc(plumedCommunications)
    pydoc.writedoc(plumedCommunications.defaults)
    pydoc.writedoc(plumedCommunications.PythonCVInterface)
    return {"Value":plumedCommunications.defaults.COMPONENT_NODEV, "ATOMS":"1"}

def plumedCalculate(_):
    return 0
