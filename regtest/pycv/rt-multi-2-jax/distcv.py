import jax.numpy as np
from jax import jacrev, jit


# Define the distance function
@jit
def cv_f(X):
    return {
     'd12': np.linalg.norm( X[0,:]-X[1,:] ),
     'd13': np.linalg.norm( X[0,:]-X[2,:] )
     }

cv_j=jacrev(cv_f)


# The CV function actually called
def cv(X):
    return cv_f(X), cv_j(X)

