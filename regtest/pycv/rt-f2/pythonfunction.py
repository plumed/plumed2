#import jax.numpy as np
from jax import grad,jit

@jit
def cc1f(xx):
    x, y, z = xx
    return x**2+y+0*z+1
cc1g=grad(cc1f)
def cc1(x): return cc1f(x), cc1g(x)

@jit
def cc2f(xx):
    x, y, z = xx
    return 0*x**3+0*y**3+0*z**3
cc2g=grad(cc2f)
def cc2(x): return cc2f(x), cc2g(x)

@jit
def cc3f(xx):
    x, y = xx
    return x**2+y**2
cc3g=grad(cc3f)
def cc3(x): return cc3f(x), cc3g(x)


