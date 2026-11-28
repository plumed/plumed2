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


> [!NOTE]
>We have some problem with the macs:
you may need to explicitly declare and export some environmental variables like the `PYTHONPATH` to make pycv work when plumed is called


## Documentation

>[!NOTE]
>The regtests `regtest/pycv/rt-doc` and `regtest/pycv/rt-func-doc` can be used both as
 an example and a way to create a (very) basic code documentation for the action/function interface.
>The example in `regtest/pycv/rt-doc` will produce very basic html pages,
 and the one in `regtest/pycv/rt-func-doc` will dump a few simple textfiles


The PYCV module defines the following actions:

 * `PYCVINTERFACE`, to implement single- and multi-component CVs in Python;
 * `PYFUNCTION`, to implement arbitrary functions.

Plumed will start an embedded Python interpreter (or use the one of the
Python program, when PLUMED runs within Python).
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

pycv embeds a Python interpreter, so it is built for one Python version: usually
as an extension against an installed PLUMED, once for each Python version that you
want to use, or in the PLUMED kernel, for one Python version (see below).

With conda, install the package `plumed-ext-pycv` next to `plumed`: the variant
that matches the Python of the environment is selected automatically.

Otherwise, from the source tree of the same version of the installed PLUMED,
with a Python that has `pybind11` (>=3.0) and its shared library, and `numpy` at runtime:
```sh
./configure --with-extension-base=/path/to/plumed --enable-extension-modules=pycv \
  PYTHON_BIN=python3.12 --prefix=/path/to/pycv
make
make install
```
(see `extensions/README.md`). The installed module is loaded automatically
(add its `lib/plumed/plugins` directory to `PLUMED_PLUGIN_PATH` if it is not
installed next to PLUMED), or with `LOAD MODULE=pycv` when autoloading is
disabled (`PLUMED_AUTOLOAD=0`). `plumed modules --has pycv` checks that it is available.

When only one Python version is needed (e.g. on a cluster), pycv can instead be
compiled into the kernel, requesting it by name:
```sh
./configure --enable-modules=+pycv PYTHON_BIN=python3.12
```
The kernel then links the shared libpython of that Python (`PYTHON_BIN`, by default
`python3`, which must have `pybind11` >=3.0), and works only with that Python
version; `--enable-modules=all` does not include pycv. When PLUMED runs within a
Python program whose executable contains its own interpreter (as the Python of
Debian, Ubuntu and conda), the kernel must be loaded without `RTLD_DEEPBIND`
(`export PLUMED_LOAD_NODEEPBIND=yes`, otherwise pycv stops with an error):
with `RTLD_DEEPBIND` the kernel would use the libpython it links, a second
interpreter. The extension has no such limitation.

> [!NOTE]
>Until PLUMED 2.10 pycv was a plugin, installed with `pip install plugins/pycv`
 and loaded with `LOAD FILE=$(python -m pycv)`: these `LOAD FILE` lines are not
 needed anymore, and can be removed from the inputs.

## Regression tests

A suite of regression tests is provided in `regtest/pycv`, run with
```sh
make -C regtest/pycv
```
when pycv is available (e.g. after the installation above, with `PLUMED_PLUGIN_PATH`).
The tests in `regtest/pycv/pythontests` run PLUMED from Python (they need the `plumed`
Python package) and are run with `pytest`.

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
>Similarly, some Python configurations (e.g. MacOS, or a virtual environment) ignore the current
 environment when searching for packages (e.g. `numpy`): the embedded interpreter is
 the Python pycv was built for. Hence,
 one should set PYTHONPATH manually.  This manifests itself with an
 error like:

```
      No module named numpy
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
    action.data["count"] = 0.0
    return {"NOPBC": True, "Value": PLMD.defaults.COMPONENT_NODEV}


def plumedCalculate(action: PLMD.PythonCVInterface):
    # ...something is defined...
    if g(something):
        action.data["example"] += 1.0
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
    # ...complex calculation are called here...
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

In the following paragraphs I am showing a few examples, for a short description of the python interface look [below](#the-python-interface)


### Bare minimum
The minimum invocation can be the following:

**plumed.dat**
```
cvPY: PYFUNCTION IMPORT=pycv
PRINT FILE=colvar.out ARG=*
```

**pycv.py**
```python
import plumedCommunications as PLMD

plumedInit = {"Value": PLMD.defaults.COMPONENT_NODEV}


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

plumedInit = {"Value": plumedCommunications.defaults.COMPONENT}


# Implementation of the angle function. @jit really improves speed
@jit
def angle(x):
    r1 = x[0, :] - x[1, :]
    r2 = x[2, :] - x[1, :]

    costheta = np.dot(r1, r2) / np.sqrt(np.dot(r1, r1) * np.dot(r2, r2))
    # OR: costheta = np.dot(r1,r2) / np.linalg.norm(r1) / np.linalg.norm(r2)
    theta = np.arccos(costheta)
    return theta


# Use JAX to auto-gradient it
grad_angle = jit(grad(angle))


# The CV function actually called
def cv1(action):
    x = action.getPositions()
    return angle(x), grad_angle(x)
```

## The Python interface

### Getting the manual
To obtain a <i>**VERY BASIC**</i> reference/manual with all the function and method definitions, run the test rt-doc or call `PYFUNCTION` in with: 

**plumed.dat**
```
PYFUNCTION IMPORT=pycv
```

**pycv.py**
```python
import plumedCommunications
import pydoc


def plumedInit(_):
    pydoc.writedoc(plumedCommunications)
    pydoc.writedoc(plumedCommunications.defaults)
    return {"Value": plumedCommunications.defaults.COMPONENT_NODEV}


def plumedCalculate(_):
    return 0.0
```

### Specifics
In the following section
#### Common interface
Both `PYFUNCTION` and `PYCVINTERFACE` have the following attribute:
 - `label` (readonly) returns the label

And the following functions:
 - `log(s:object)` puts a string in the PLUMED output
 - `lognl(s:object)` puts a string in the PLUMED output (and appends a newline)
 - `getStep()` Returns the current step
 - `getTime()` Return the present time
 - `getTimeStep()` Return the timestep
 - `isExchangeStep()` Check if we are on an exchange step
 - `isCalculationOnly()` In a calculate callback, check whether this is an auxiliary evaluation that must not change persistent state
 - `isRestart()` Return true if we are doing a restart

Python callbacks must guard their own persistent state, including random-number
generation, with `isCalculationOnly()`. This also covers auxiliary evaluations
that are not replica exchanges; `isExchangeStep()` alone is insufficient.
Calculation-only callbacks must still provide values and derivatives required
by their interface: auxiliary shell evaluations require forces.
The calculation-only context is established for calculation callbacks; it is not
a promise that an earlier prepare callback already knows about an unflagged trial.

#### PYFUNCTION: arguments

`PYFUNCTION` accepts arguments in the plumed file with the `ARG` keyword, the arguments are then accessible in python with the functions:
 - `PythonFunction.argument(argID:int)` Get value of the of argID-th argument (as `float`)
 - `PythonFunction.arguments()` Retuns a ndarray with the values of the arguments (as `numpy.ndarray[numpy.float64]`)

and the argument:
- `PythonFunction.nargs` (readonly) Get the number of arguments

`PYFUNCTION` also has a few functions that can be used to interact with the periodicity of the arguments:
 - `PythonFunction.bringBackInPbc(argID:int,x:float)` Takes one value and brings it back into the pbc of argument argID
 - `PythonFunction.difference(argID:int,x:float,y:float)` Takes the difference taking into account pbc for argument argID

#### PYCVINTERFACE: atoms, pbc, neigbourlist, makeWhole
`PYCVINTERFACE` works with atoms:

Has the following attributes:
 - `PythonCVInterface.nat` (readonly) Return the number of atoms

For getting the atomic data you can call the following functions:
 - `PythonCVInterface.getPosition(atomID:int)` Returns an ndarray with the position of the atomID-th atom 
 - `PythonCVInterface.getPositions()` Returns a numpy.array that contains the atomic positions of the atoms 
 - `PythonCVInterface.mass(atomID:int)` Get mass of atomID-th atom
 - `PythonCVInterface.masses()` Returns and ndarray with the masses 
 - `PythonCVInterface.charge(atomID:int)` Get charge of atomID-th atom
 - `PythonCVInterface.charges()` Returns and ndarray with the charges 
 - `PythonCVInterface.absoluteIndexes()` Get the npArray of the absolute indexes (like in AtomNumber::index()).

 And in general you can use some support-function from plumed:
 - `PythonCVInterface.getNeighbourList()` returns an interface to the current Neighborlist
 - `PythonCVInterface.getPbc()` returns an interface to the current pbcs
 - `PythonCVInterface.makeWhole()` Make atoms whole, assuming they are in the proper order

`plumedCommunications.Pbc` and `plumedCommunications.NeighborList` are simple methods container to help with the calculations:

`plumedCommunications.Pbc` has the following functions:
 - `Pbc.apply(d:ndarray)` pply PBC to a set of positions or distance vectors (`d.shape==(n,3)`)
 - `Pbc.getBox()` Get a numpy array of shape (3,3) with the box vectors
 - `Pbc.getInvBox()` Get a numpy array of shape (3,3) with the inverted box vectors

`plumedCommunications.NeighborList` has the following functions and attributes:
- `NeigbourList.getClosePairs()` get the (nl.size,2) nd array with the list of couple indexes
- `NeigbourList.size` and `len(nl:NeigbourList.size)` give the number of couples

## Limitations

- There is one Python interpreter per process, shared by all the pycv actions
  of all the Plumed objects, in all the threads: e.g. the thread-MPI ranks (as in
  GROMACS without real MPI, or `plumed --tmpi-threads`), or Plumed objects used by
  different threads. A module is imported once, and its global variables are
  shared by all the actions that import it, also of different replicas: keep the
  state of an action in its `data` dictionary (PYCVINTERFACE). The calls to Python
  hold the GIL, so the Python code of different threads does not run at the same
  time. With MPI processes every process has its own interpreter.
- When pycv starts the interpreter (PLUMED is not running within Python), it does
  so in a thread of its own, without the Python signal handlers (e.g. ctrl-c stops
  the program as without pycv), and finalizes it when the last pycv action is
  destroyed (or at exit): `threading.main_thread()` is that thread, so Python code
  that must run in the main thread (e.g. `signal.signal`) cannot be used.
- pycv is compiled for one Python version, and it stops with an error if the
  Python interpreter is another version (e.g. when PLUMED runs within Python).
- Once a pycv action has been created, the library of pycv (and so the PLUMED
  kernel and libpython) is never unloaded, as any Python extension: when PLUMED
  runs within Python, the Python code can keep `plumedCommunications`.
- `plumedCommunications` is provided by pycv, not installed as a Python package:
  when PLUMED runs within Python, the Python script can import it (or a module
  that imports it) only after PLUMED has created a pycv action.
- JAX's GPU/TPU offloading are not 100% tested.


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
