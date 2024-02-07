# This is only a partial example. This function does not compute the
# gradient, so it is useless for biasing. See the other regression
# tests for how to auto-grad.

# And, of course, one should not call slow functions (such as print)
# in the CV calculation.

import plumedCommunications
import pydoc

def plumedInit(_):
    with open('PythonCVInterface.help.txt', 'w') as f:
        h = pydoc.Helper(output=f)
        h(plumedCommunications.PythonCVInterface)
    with open('plumedCommunications.help.txt', 'w') as f:
        h = pydoc.Helper(output=f)
        h(plumedCommunications)
    with open('plumedCommunications.defaults.help.txt', 'w') as f:
        h = pydoc.Helper(output=f)
        h(plumedCommunications.defaults)
    return {"Value":plumedCommunications.defaults.COMPONENT_NODEV, "ATOMS":"1"}

def plumedCalculate(_):
    return 0
