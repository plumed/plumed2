# The PYCV module/plugin for PLUMED 2

The [PYCV module](https://giorginolab.github.io/plumed2-pycv) enables
PLUMED2 Collective Variables (CVs) and arbitrary functions to be
defined and auto-differentiated in the Python language.

Advantages of using PYCV over standard development of CVs in C++ are:
 1. functions may be prototyped in  high-level code, using
    extensive mathematical libraries, and no boilerplate;
 2. just-in-time compilation
    occurs transparently: there are no compilation and link delays
    for code changes;
 3. CVs may be automatically differentiated in common cases.

The code is organized as a standard PLUMED module: development occurs
in a [fork of the original
repository](https://github.com/giorginolab/plumed2-pycv/tree/v2.5.2-pycv/src/pycv). All
code is in the `src/pycv` and `regtest/pycv` directories.

[![DOI](https://joss.theoj.org/papers/10.21105/joss.01773/status.svg)](https://doi.org/10.21105/joss.01773)
[![plumID:19.075](https://www.plumed-nest.org/eggs/19/075/badge.svg)](https://www.plumed-nest.org/eggs/19/075/)


## Documentation

The PYCV module defines the following actions:

 * `PYCVINTERFACE`, to implement single- and multi-component CVs in Python;
 * `PYFUNCTION`, to implement arbitrary functions.

Plumed will start a a Python interpreter.
Then Plumed will import a python module with the `IMPORT=` keyword, this module
must contain at least two objects: a calculate function that will be called at
the `calculate()` step, and an init function or dictionary that will be used at
time of constructing the pycv.
The module that can be imported can be a `*.py` file or a directory that contains
an `__init__.py`, or a module in the your python path.

`PYCVINTERFACE` will accept also functions during the `prepare()` and `update()`
steps.

`PYCVINTERFACE` will accept also the `ATOMS` keyword (see `DISTANCE` cv)
**or** a series of keyword relative to the neighbour list (see `COORDINATION` cv).
If both the flavour of asking atom is used, `PYCVINTERFACE` will raise and error.

The interface to python as anticipated before depends on the following keywords:

| keyword   | description                                  | PYCVINTERFACE | PYFUNCTION |
|-----------|----------------------------------------------|---------------|------------|
| IMPORT    | the module to import                         | ✅            | ✅         |
| INIT      | the function/dict to call/get @ construction | ✅            | ✅         |
| CALCULATE | the function to call at `calculate()` step   | ✅            | ✅         |
| PREPARE   | the function to call at `prepare()` step     | ✅            | ❌         |
| UPDATE    | the function to call at `update()` step      | ✅            | ❌         |

If not specified INIT will default to `"plumedInit` and CALCULATE  to 
`"plumedCalculate"`, on the other hand, if PREPARE and UPDATE are not present, will be ignored.

## Preparation

For compiling the plugin you just need pybind11 and numpy.
I always recomend to create and ad-hoc environment for your projects:
```bash
python3 -m venv pycvenv
source ./pycvenv/bin/activate
pip install -U pip
pip install -r requirements.txt
```
The requirements.txt file is in the home of the plug in

### Standard compilation

If you have a plumed that supports plumed mklib (that will be release in the 2.10 version, but it is avaiable in the master branch) with multiple files you can simply
```bash
./standaloneCompile.sh
```

### Developer compilation

If your plumed version is inferior to 2.10, or you want to contribute to this module,
the procedure is slighly more compex:
```bash
./prepareMakeForDevelop.sh
```
will prepare a Make.inc in this directory that will be included by the Makefile.
Then simply:
```bash
make
```

#### Set up tests

If you are interested in running the test regarding this plugin you can use the same procedure as with the standard plumed, but in the subdir regtest of this plugin.
The only requirement is to copy or to create a symbolic link to the `regtest/scripts` directory in a plumed source. Plumed must be runnable to execute tests
```bash
cd regtest
ln -s path/to/plumed/source/regtest/scripts .
make
```
### About older Plumed versions

If you are using an older plumed version you must know that:
 - On linux the plug-in can be loaded only if `LOAD` support the GLOBAL keyword
 - mklib won't work, so you'll need to use `./prepareMakeForDevelop.sh`

## Quickstart

Here's a quick example to whet your appetite, following the regression test `rt-jax2`.

    export PLUMED_PROGRAM_NAME=$PWD/src/lib/plumed-static 
    make -C regtest/pycv/rt-jax2

The files in `regtest/pycv/rt-jax2` are reproduced here:

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

## EXTRA: JAX

Transparent auto-differentiation, JIT compilation, and vectorization
are available through Google's [JAX
library](https://github.com/google/jax) (recommended).

### Install jax:

Go to the original guide in the [jax documenation](https://jax.readthedocs.io/en/latest/installation.html)

jax has different method of installation, and can be accelerated with various different hardware,
(as stated before, trying to install things in a virtual environment make doing error less costly)
The command for installing should be similar to:
 - example if you have a cuda12 compatible device (a wheel for cuda will be installed alognside jax):
`pip install "jax[cuda12_pip]" -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html`
 - example if you have a cuda12 compatible device, and **cuda already installed on your system**:
`pip install "jax[cuda12_local]" -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html`

## Demos and regression tests

A number of PLUMED-style regression tests are provided to test and
demonstrate the package features.


Name   | Feature
-------|------------
`rt-1` | Basic test (distance of two atoms), pure Python, no gradient
`rt-2` | CV and gradient explicitly coded, pure Python
`rt-3` | Multiple PYTHON actions in the same file
`rt-f1`| Function, pure Python
`rt-f2`| Function, with JIT and auto-gradient, with MATHEVAL equivalent
`rt-jax1` | Distance CV with JAX reverse-mode differentiation
`rt-jax2` | Angle CV with JAX reverse-mode differentiation and JIT decorator
`rt-jax3` | Radius of curvature, as in [doi:10.1016/j.cpc.2018.02.017](http://doi.org/10.1016/j.cpc.2018.02.017)
`rt-multi1` | Multi-component CV, pure Python
`rt-multi2` | Multi-component CV, with JAX reverse-mode differentiation


To run:

```bash
cd [path_to_repository]
export PLUMED_PROGRAM_NAME=$PWD/src/lib/plumed
cd regtest/pycv
make
```

## Common errors

* Missing modules at run time: you may need to set the `PYTHONHOME`
  environment variable.

* `uncaught exception [...] The interpreter is already running` This
  occurs if the module is loaded twice. The most common cause is when
  it is both built-in in the current kernel, and loaded as a module.



## Limitations

- No test have been done with MPI
- JAX's GPU/TPU offloading are unde.



## Authors

Toni Giorgino <toni.giorgino@gmail.com>


## Contributing

Please report bugs and ideas via this repository's *Issues*. 


## Citation

Giorgino T. PYCV: a PLUMED 2 Module Enabling the Rapid Prototyping of
Collective Variables in Python. The Journal of Open Source Software
4(42):1773 

[![DOI](https://joss.theoj.org/papers/10.21105/joss.01773/status.svg)](https://doi.org/10.21105/joss.01773)
[![plumID:19.075](https://www.plumed-nest.org/eggs/19/075/badge.svg)](https://www.plumed-nest.org/eggs/19/075/)


## Copyright

PYCV is distributed under the LGPL terms: see COPYRIGHT.

PYCV includes the [pybind11](https://github.com/pybind/pybind11)
library, distributed under its own license terms (BSD-3).
