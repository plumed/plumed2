The PYCV module for PLUMED 2
====================================

Enables the implementation of PLUMED collective variables (`PYTHONCV`
action) or functions (`PYTHONFUNCTION` action) in the Python language.

Transparent auto-differentiation and JIT compilation are provided
through Google's [JAX library](https://github.com/google/jax).


Install
------------------------------------

Follow the usual PLUMED 2 configuration procedure with Python-specific
flags:

```
./configure --enable-modules=+pycv --enable-python PYTHON_BIN=python3 LDFLAGS="`python3-config --ldflags`"
```

It is also possible to compile the module as a stand-alone dynamic
object.  Once in the `src/pycv` directory, try `make PYCV.so` (`make
PYCV.dylib` under OSX). The compilation step *should* pick
Python-specific compilation and linker flags as long as the
`python3-config` executable is in your path.  Use Plumed's `LOAD`
action to load the generated object.


Prerequisites
------------------------------------

Python 3 is required. The version distributed with Conda is known to
create problems (slow compilation, linking failures), particularly
under OSX. There, [Homebrew](https://brew.sh)'s Python 3 package is
known to work. You will also need to install `numpy`, either via `pip3
install numpy` or your distribution's packages.

Automatic differentiation examples require the JAX library: install
it with `pip3 install jax`.

At run time, you may need to set the `PYTHONHOME` environment
variable.



Documentation
------------------------------------
Please see the generated documentation and regression tests. 


Authors
------------------------------------
Toni Giorgino <toni.giorgino@cnr.it>


Copyright
------------------------------------
Distributed under the LGPL terms: see COPYRIGHT.

The software includes the
[pybind11](https://github.com/pybind/pybind11) library, distributed
under its own license terms (BSD-3).


