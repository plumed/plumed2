The PYCV module for PLUMED 2
====================================

Enable the implementation of PLUMED collective variables (CVs) and
functions in the Python language.

Advantages of using PYCV over standard development of CVs in C++ are:
 1. functions may be prototyped in  high-level code, using
    extensive mathematical libraries, and no boilerplate;
 2. just-in-time compilation
    occurs transparently: there are no compilation and link delays
    for code changes;
 3. CVs may be automatically differentiated in common cases.



Documentation
------------------------------------

The PYCV module defines the following actions:

 * `PYTHONCV`, to implement single- and multi-component CVs in Python;
 * `PYTHONFUNCTION`, to implement arbitrary functions.

In both cases, a Python interpreter is first started; the Python file
indicated in the `IMPORT=` keyword is then imported; from it, an
user-chosen function (`FUNC=`) is called to perform the computations
at each timestep. Modules can be shared for multiple functions, and
contain one-time initialization.

Transparent auto-differentiation, JIT compilation, and vectorization
are available through Google's [JAX
library](https://github.com/google/jax) (recommended).

Please see the generated documentation and regression tests. 



Prerequisites
------------------------------------

Python 3 is required. Under OSX, [Homebrew](https://brew.sh)'s Python
3 package works well. The Conda version is problematic (slow
compilation, linking failures).

You will also need to install `numpy`, either via `pip3 install
numpy` or your distribution's packages.

Automatic differentiation examples require the JAX library: install
it with `pip3 install jax`.



Installation
------------------------------------

Follow the usual PLUMED 2 configuration procedure, with added
Python-specific flags:

```
./configure --enable-modules=+pycv --enable-python PYTHON_BIN=python3 LDFLAGS="`python3-config --ldflags`"
```

It is also possible to compile the module as a `LOAD`-able dynamic
object.  Once in the `src/pycv` directory, issue `make PYCV.so`
(`PYCV.dylib` under OSX). The compilation step *should* pick
Python-specific flags as long as the correct `python3-config`
executable is in your path.


Quickstart
------------------------------------

Here's a quick example to whet your appetite (regression test `rt-jax2`).

**File plumed.dat**

```
cv1:  PYTHONCV ATOMS=1,4,3 IMPORT=jaxcv FUNCTION=cv1
```

**File jaxcv.py**

```py
# Import the JAX library
import jax.numpy as np
from jax import grad, jit, vmap

# Implementation of the angle function. @jit really improves speed
@jit
def angle(x):
    r1 = x[0,:]-x[1,:]
    r2 = x[2,:]-x[1,:]

    costheta = np.dot(r1,r2) / np.sqrt(np.dot(r1,r1) * np.dot(r2,r2))
    # OR: costheta = np.dot(r1,r2) / np.linalg.norm(r1) / np.linalg.norm(r2)
    theta = np.arccos(costheta)
    return theta

# Use JAX to auto-gradient it
grad_angle = grad(angle)

# The CV function actually called
def cv1(x):
    return angle(x), grad_angle(x)

```





Demos and regression tests
------------------------------------

A number of PLUMED-style regression tests are provided to test and
demonstrate the package features.

Summary of the provided regression tests:

Name   | Feature
-------|------------
`rt-1` | Basic test (distance of two atoms), pure Python, no gradient
`rt-2` | CV and gradient explicitly coded, pure Python
`rt-3` | Multiple PYTHON actions in the same file
`rt-f1`| Function, pure Python
`rt-f2`| Function, with JIT and auto-gradient, with MATHEVAL equivalent
`rt-jax1` | Distance CV with JAX reverse-mode differentiation
`rt-jax2` | Angle CV with JAX reverse-mode differentiation and JIT decorator
`rt-multi1` | Multi-component CV, pure Python
`rt-multi2` | Multi-component CV, with JAX reverse-mode differentiation


To run:

```
cd [path_to_repository]
export PLUMED_PROGRAM_NAME=$PWD/src/lib/plumed
cd regtest/pycv
make
```





Common errors
------------------------------------

* Missing modules at run time: you may need to set the `PYTHONHOME`
  environment variable.

* `uncaught exception [...] The interpreter is already running` This
  occurs if the module is loaded twice. The most common cause is when
  it is both built-in in the current kernel, and loaded as a module.



Limitations
------------------------------------

Behavior with MPI, and JAX's GPU/TPU offloading are untested.



Author
------------------------------------

Toni Giorgino <toni.giorgino@cnr.it>


Copyright
------------------------------------

PYCV is distributed under the LGPL terms: see COPYRIGHT.

PYCV includes the [pybind11](https://github.com/pybind/pybind11)
library, distributed under its own license terms (BSD-3).


