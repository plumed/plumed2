# Import the JAX library
import jax.numpy as np
from jax import grad, jit

# Define the distance function
@jit
def dist(x):
    d = x[0,:]-x[1,:]
    d2 = np.dot(d,d)
    return np.sqrt(d2)

# Use JAX to auto-gradient it
grad_dist = grad(dist)

# The CV function actually called
def cv1(x):
    return dist(x), grad_dist(x)
