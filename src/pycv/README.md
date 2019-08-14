The PYCV module for PLUMED 2
====================================

Enable the implementation of PLUMED collective variables (CVs) and
functions in the Python language.

Advantages of using PYCV over standard development of CVs in C++ are:
 1. code may be prototyped in a high-level language, with access to
    extensive mathematical libraries, and almost no boilerplate;
 2. there are no compilation and link delays (just-in-time compilation
    occurs transparently);
 3. the required gradients of CVs can be
    computed autmatically.



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
library](https://github.com/google/jax).

Please see the generated documentation and regression tests. 



Prerequisites
------------------------------------

Python 3 is required. Under OSX, [Homebrew](https://brew.sh)'s Python
3 package is known to work well. The version distributed with Conda is
problematic (slow compilation, linking failures).

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


Testing
------------------------------------

PLUMED-style regression tests demonstrate and check the package
features. To run them:

```
cd [path_to_repository]
export PLUMED_PROGRAM_NAME=$PWD/src/lib/plumed
cd regtest/pycv
make
```

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




Common errors
------------------------------------

* Missing modules at run time: you may need to set the `PYTHONHOME`
  environment variable.

* `uncaught exception [...] The interpreter is already running` This
  occurs if the module is loaded twice. The most common cause is when
  it is both built-in in the current kernel, and loaded as a module.



Limitations
------------------------------------

Interaction with MPI is undefined.



Author
------------------------------------

Toni Giorgino <toni.giorgino@cnr.it>


Copyright
------------------------------------

PYCV is distributed under the LGPL terms: see COPYRIGHT.

PYCV includes the [pybind11](https://github.com/pybind/pybind11)
library, distributed under its own license terms (BSD-3).


