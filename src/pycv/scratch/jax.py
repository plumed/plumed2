import jax
import jax.numpy as np

def d2(x):
    return x.dot(x)

d2grad = jax.grad(d2)

