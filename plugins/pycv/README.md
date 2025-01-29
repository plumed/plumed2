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
[![plumID:19.075](https://www.plumed-nest.org/eggs/19/075/badge.svg)](https://www.plumed-nest.org/eggs/19/075/). 

Note that the current syntax is different than the one described in the original paper.


As usual if you are not working in a virtual environment you are doing this at your own risk
`pip install .` should be enough, with a plumed version (>=2.10) avaiable with a pkg-config or in the `PATH`s.


>[!NOTE]
>On linux If you use `make` and `make check` note that the tests will find pycv only if you are working in the same python environment that you had active when compiling/installing plumed

> [!NOTE]
>We have some problem with the macs:
you may need to explicitly declare and export some environmental variables like the `PYTHON_PATH` to make pycv work when plumed is called


## Documentation

>[!NOTE]
>The regtests `regtest/pycvcomm/rt-doc` and `regtest/pycvfunc/rt-doc` can be used both as
 an example and a way to create a (very) basic code documentation for the action/function interface.
>The example in `regtest/pycvcomm/rt-doc` will produce very basic html pages,
 and the one in `regtest/pycvfunc/rt-doc` will dump a few simple textfiles


The PYCV module defines the following actions:

 * `PYCVINTERFACE`, to implement single- and multi-component CVs in Python;
 * `PYFUNCTION`, to implement arbitrary functions.

Plumed will start an embedded Python interpreter.
Then Plumed will import a python module with the `IMPORT=` keyword. This module
must contain at least two objects: a *calculate* function that will be called at
the `calculate()` step, and an *init* function or dictionary that will be used at
time of constructing the PYCV.
The module that can be imported can be a `*.py` file or a directory that contains
an `__init__.py`, or a module in the your Python path.

`PYCVINTERFACE` will accept also functions during the `prepare()` and `update()`
steps.

`PYCVINTERFACE` will accept also the `ATOMS` keyword (see `DISTANCE` CV)
**or** a series of keyword relative to the neighbour list (see `COORDINATION` CV).
If both the flavour of asking atom is used, `PYCVINTERFACE` will raise and error.

The interface to Python, as anticipated above, depends on the following keywords:

| keyword   | description                                  | PYCVINTERFACE | PYFUNCTION |
|-----------|----------------------------------------------|---------------|------------|
| IMPORT    | the module to import                         | ✅            | ✅         |
| INIT      | the function/dict to call/get @ construction | ✅            | ✅         |
| CALCULATE | the function to call at `calculate()` step   | ✅            | ✅         |
| PREPARE   | the function to call at `prepare()` step     | ✅            | ❌         |
| UPDATE    | the function to call at `update()` step      | ✅            | ❌         |

If not specified, `INIT` will default to `plumedInit` and `CALCULATE`  to 
`plumedCalculate`. On the other hand, if `PREPARE` and `UPDATE` are not present, they will be ignored.

## Installation

It should be sufficient to run, ideally in a virtual or conda environment:
```sh
cd plugins/pycv
pip install .
```
You will need to have plumed avaiable in your path or through `pkg-config`

Required dependencies `numpy` and `pybind11` are installed as part of the installation process.

Note that an in-place installation, `pip install -e .`, won't work. 

Now you can get the position of the shared object to load with `python -m pycv`, see [below](#some-examples)

## Regression tests

A suite of regression tests are provided in the `regtest` subdirectory. They can be run e.g. with 
```sh
make -C regtest
```

## Common runtime problems

On some platforms, *embedded* Python interpreters  (such as the one used in PYCV) appear to behave 
differently than the plain ones, raising surprising errors.  For example:

>[!NOTE]
>Some Python configurations (e.g. Conda under Linux) require the
 Python shared library to be found in the LD_LIBRARY_PATH, 
 ignoring the activated environment.  This manifests itself with an
 error like:

```
      libpython3.13.so.1.0: cannot open shared object file: No such file or directory
```

>[!NOTE]
>Similarly, some Python configurations (e.g. MacOS) ignore the current
 environment when searching for packages (e.g. `numpy`). Hence,
 one should set PYTHONPATH manually.  This manifests itself with an
 error like:

```
      No module named numpy
```
### Developing pycv

If you are developing pycv can be compiled and tested from the pycv dir with `make`

```sh
make
make check
```

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

>![IMPORTANT]
> In the following examples you will see often `path/to/PythonCVInterface.so` as the path to load pycv into plumed
> If you are using a standalone plumed you can get the path to that file simply by calling `python -m pycv`
> If you are using plumed within python by importing the `pycv` python module
  and loading the pycv plumed module with the convenient function  `pycv.getPythonCVInterface()`,
  like `plmd.cmd("readInputLine", f"LOAD FILE={pycv.getPythonCVInterface()}")`


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

## JAX 

Transparent auto-differentiation, JIT compilation, neural networks and vectorization
are readily available through Google's [JAX
library](https://github.com/google/jax) (recommended).

### Installation

As described in the  [jax documenation](https://jax.readthedocs.io/en/latest/installation.html), there are several installation routes, and calculations can be accelerated with various hardware. For example:
 - example if you have a cuda12 compatible device (a wheel for cuda will be installed alongside jax):
`pip install "jax[cuda12_pip]" -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html`
 - example if you have a cuda12 compatible device, and **cuda already installed on your system**:
`pip install "jax[cuda12_local]" -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html`



### Example: a CV with automatic gradient calculation 

Here's a numpy-style calculation of an angle between three atoms, with automatic differentiation:

**plumed.dat**

```
cv1:  PYCVINTERFACE ATOMS=1,4,3 IMPORT=jaxcv CALCULATE=cv1
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
grad_angle = jit(grad(angle))

# The CV function actually called
def cv1(action):
    x=action.getPositions()
    return angle(x), grad_angle(x)

```

## Limitations

- No test have been done with MPI
- JAX's GPU/TPU offloading are not 100% tested.


If you are using an older Plumed version you must know that:
 - On linux the plug-in can be loaded only if `LOAD` supports the `GLOBAL` keyword
 - mklib won't work (supports only single file compilations), so you'll need to use `./prepareMakeForDevelop.sh`

## Authors

Original author: Toni Giorgino <toni.giorgino@cnr.it>

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
