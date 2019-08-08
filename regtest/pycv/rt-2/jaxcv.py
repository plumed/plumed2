
import jax.numpy as np
from jax import grad, jit, vmap


log=open("pycv.log","w")

print("Imported.",file=log)


def dist(x):
    d = x[0,:]-x[1,:]
    d2 = np.dot(d,d)
    return np.sqrt(d2)

grad_dist = grad(dist)

def cv1(x):
    print(x,file=log)
    return dist(x), grad_dist(x)
