# The PYCV module/plugin for PLUMED 2

The PYCV module enables
PLUMED2 Collective Variables (CVs) and arbitrary functions to be
defined and auto-differentiated in the Python language.

Advantages of using PYCV over standard development of CVs in C++ are:
 1. functions may be prototyped in  high-level code, using
    extensive mathematical libraries, and no boilerplate;
 2. just-in-time compilation
    occurs transparently: there are no compilation and link delays
    for code changes (using external tools such as JAX);
 3. CVs may be automatically differentiated in common cases. (using external tools such as JAX)

You can see the original PyCV [here](https://giorginolab.github.io/plumed2-pycv)
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

If you want to contribute to this module,
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
 - On linux the plug-in can be loaded only if `LOAD` supports the `GLOBAL` keyword
 - mklib won't work (supports only single file compilations), so you'll need to use `./prepareMakeForDevelop.sh`

## Getting started

### Initialization: the INIT keyword

Both `PYFUNCTION` and `PYCVINTERFACE` use the keyword INIT to finalize the set up for the action.
If not specified INIT wil default to `"plumedInit"`.
INIT instruct plumed to call a function that returns a dict or read directly a dict.

The function must accept a `plumedCommunications.PythonFunction` or a `plumedCommunications.PythonCVInterface` object that can be used to interact with plumed (this will be explained better in the CALCULATE section)
The init dict must contain at least the "Value" or the "COMPONENTS" keys: 
 - to the "Value" key must be assigned a value-dict
 - to  the "COMPONENTS" key must be assigne a dict whose keys are assigned to a value-dict. The keys of the COMPONTENTS-dict will be used as name of the components of the action.
What I am calling a "value-dict" is a dict with two entries: `{'derivative': bool, 'period': None or [min,max]}` "derivative" simply tells plumed that the component will returns derivatice, the "period" will set up the periodicity for that componentd ("min" an "max" are parsed with lepton, so you can also use strings containing `pi`, `e`, etc...)

To help setting up a value/components plumedCommunications has a submodule defaults that contains 2 possible defaults (see the example): 
  - `plumedCommunications.defaults.COMPONENT = {'derivative': True, 'period': None}`
  - `plumedCommunications.defaults.COMPONENT_NODEV = {'derivative': False, 'period': None}`
 
The init dictionary for `PYCVINTERFACE` can contain nearly all the keywords of the input line:
 - NOPBC -flag, needs True/False explicitly set-
 - ATOMS
 - GROUPA
 - GROUPB
 - PAIR -flag, needs True/False explicitly set-
 - NLIST -flag, needs True/False explicitly set-
 - NL_CUTOFF
 - NL_STRIDE
Note that if the keyword is both in the input line of the plumed file and in the init-dict an error will be raised.

Nothe that if you use ATOMS and the NL keywords (GROUPA, GROUPB, PAIR, NLIST, NL_CUTOFF, NL_STRIDE) an error will be raised.

`PYCVINTERFACE` has a  `.data` attribute that is a dict, that can be used to store long term data, this dict will not be accessed by plumed.

```python
import plumedCommunications as PLMD

def plumedInit(action: PLMD.PythonCVInterface):
    action.data["count"]=0.0
    return {"NOPBC": True, "Value": PLMD.defaults.COMPONENT_NODEV}

def plumedCalculate(action: PLMD.PythonCVInterface):
    #...something is defined...
    if g(something):
        action.data["example"]+=1.0
    return f(something)
```

### Calculate step: the CALCULATE keyword
Both `PYFUNCTION` and `PYCVINTERFACE` use the keyword CALCULATE to call a function in the `calculate()` step
If not specified CALCULATE will default to `"plumedCalculate"`.

Plumed will expect CALCULATE to return either a tuple or a dict:
 - the tuple can have up to 3 components: `[value, derivative, boxDerivative]`: if "derivative" is set to `True` an error will be raised if the tuple contain only one element, or if contains more than one element in the other case.
 - in the case of multiple COMPONETS plumed will expect a dict with a key per defined component. Each key  must contain a tuple with the previous criterions.

Instead of a tuple you can return a single float (or a dict of floats), but plumed will complain with a warning.

```python
import plumedCommunications as PLMD


plumedInit = {
    "COMPONENTS": {
        "first": PLMD.defaults.COMPONENT_NODEV,
        "second": {"derivative": True, "period": ["-pi", "pi"]},
    }
}

def plumedCalculate(action: PLMD.PythonCVInterface):
    #...complex calculation are called here...
    return {
        "first": [resForFirst],
        "second": [resForSecond, derivativeForSecond, boxDerivativeForSecond],
    }
```

### Prepare step: the PREPARE keyword
**Only** `PYCVINTERFACE` use the keyword PREPARE to call a function in the `prepare()` step, before `calculate()`

If not specified PREPARE will be ignored.
Python expect this PREPARE to return a dict.

If that dict contains a key called 'setAtomRequest' with a list of atoms (you can pass a string that will be parsed by plumed like the one in ATOMS) the new list will be requested.

### Update step: the UPDATE keyword
**Only** `PYCVINTERFACE` use the keyword UPDATE to call a function in the `update()` step, after `calculate()`

If not specified UPDATE will be ignored.
Python expect this UPDATE to return a dict.

As now nothing will be done after the python UPDATE call.

You can however act on the `.data` dict or access to all the parameters accessible from python, and for example update a plot or accumulate an histogram or any other thing that may be useful for your analysis/simulation.

## Some Examples

In the following paragraphs I am showing a few examples, for a short description of the python interface look [here](PythonInterface.md)

### Bare minimum
The minimum invocation can be the following:

**plumed.dat**
```
LOAD GLOBAL FILE=path/to/PythonCVInterface.so
cvPY: PYFUNCTION IMPORT=pycv
PRINT FILE=colvar.out ARG=*
```

**pycv.py**
```python
import plumedCommunications as PLMD
plumedInit={"Value":PLMD.defaults.COMPONENT_NODEV}

def plumedCalculate(action: PLMD.PythonFunction):
    action.lognl("Hello, world!")
    return 0.0
```
This simply prints an "Hello, world!" at each step of the simulation/trajectory.

### An example, gradient calculation with jax

Here's a quick example with calculation of an angle between three atoms

**plumed.dat**

```
cv1:  PYTHONCV ATOMS=1,4,3 IMPORT=jaxcv CALCULATE=cv1
```

**jaxcv.py**

```py
# Import the JAX library
import jax.numpy as np
from jax import grad, jit, vmap
import plumedCommunications

plumedInit={"Value": plumedCommunications.defaults.COMPONENT}
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
def cv1(action):
    x=action.getPositions()
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
 - example if you have a cuda12 compatible device (a wheel for cuda will be installed alongside jax):
`pip install "jax[cuda12_pip]" -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html`
 - example if you have a cuda12 compatible device, and **cuda already installed on your system**:
`pip install "jax[cuda12_local]" -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html`


## Limitations

- No test have been done with MPI
- JAX's GPU/TPU offloading are not 100% testes.

## Authors

Original author: Toni Giorgino <toni.giorgino@gmail.com>

Daniele Rapetti <Daniele.Rapetti@sissa.it>
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
