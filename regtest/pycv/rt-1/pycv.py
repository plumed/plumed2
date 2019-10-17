# This is only a partial example. This function does not compute the
# gradient, so it is useless for biasing. See the other regression
# tests for how to auto-grad.

# And, of course, one should not call slow functions (such as print)
# in the CV calculation.

import numpy as np
log=open("pycv.log","w")

print("Imported.",file=log)

def cv1(x):
    print(x,file=log)           # But don't
    d = x[0,:]-x[1,:]
    d2 = np.dot(d,d)
    # If computing a gradient, return it as a 2nd return value.
    return np.sqrt(d2)
