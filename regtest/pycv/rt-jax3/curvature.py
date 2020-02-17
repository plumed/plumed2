# Same calculation (and results) w.r.t. https://doi.org/10.1016/j.cpc.2018.02.017
#
# (first frame returning NANs was removed)
# https://github.com/tonigi/plumed2-automatic-gradients/tree/automatic-gradient-computation/regtest/curvature_codegen


# Import the JAX library
import jax.numpy as np
from jax import grad, jit

# Implementation of the angle function. @jit really improves speed
@jit
def r_f(x):
    r21 = x[0,:]-x[1,:]
    r23 = x[2,:]-x[1,:]
    r13 = x[2,:]-x[0,:]

    cos2theta = np.dot(r21,r23)**2 / (np.dot(r21,r21) * np.dot(r23,r23))
    sin2theta = 1-cos2theta
    
    R2= np.dot(r13,r13)/sin2theta/4.0
    return np.sqrt(R2)

# Use JAX to auto-gradient it
r_g = grad(r_f)

# The CV function actually called
def r(x):
    return r_f(x), r_g(x)


# Differentiation is trivial here
@jit
def inv(x):
    return 1/x, -(x**-2)

