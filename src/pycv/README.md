The PYCV (Python-CV) module for PLUMED 2
====================================

Enables the implementation of PLUMED collective variables (`PYTHONCV`
action) or functions (`PYTHONFUNCTION` action) in the Python language.

Transparent auto-differentiation and JIT compilation are available
through Google's [JAX library](https://github.com/google/jax).


Install
------------------------------------

Follow the usual PLUMED 2 configuration procedure with Python-specific
flags:

```
./configure --enable-modules=+pycv --enable-python PYTHON_BIN=python3 LDFLAGS="`python3-config --ldflags`"
```

Python 3 is required. Conda is known to create problems (slow
compilation, linking failures) particularly under OSX. There,
[Homebrew](https://brew.sh) Python 3 is known to work.

Automatic differentiation examples require the JAX library: install
with `pip3 install jaxlib`.

At run time, you may need to set the `PYTHONHOME` 
environment libraries.


It may be useful to compile the variable as a stand-alone dynamic
object.  Once in the `src/pycv` directory, try `make PYCV.so` (or
`make PYCV.dylib` on OSX). The compilation step *should* pick
Python-specific compilation and linker flags.  Use Plumed's \ref LOAD
action to load the generated object.




Documentation
------------------------------------
Please see the generated documentation and regression tests. 


Authors
------------------------------------
Toni Giorgino <toni.giorgino@cnr.it>


Copyright
------------------------------------
See COPYRIGHT. This software relies on the
[pybind11](https://github.com/pybind/pybind11) library, which is
distributed under its own license terms (BSD-3).


