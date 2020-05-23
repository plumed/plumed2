The PYCV module for PLUMED 2
====================================

Note: this is a fork of the official [PLUMED
repository](https://github.com/plumed/plumed2). Please see [PLUMED's
homepage](https://www.plumed.org/) for information on the PLUMED
library.


The [PYCV module](https://giorginolab.github.io/plumed2-pycv)
enables PLUMED2 Collective Variables (CVs) and arbitrary functions to
be defined and auto-differentiated in the Python language.

Advantages of using PYCV over standard development of CVs in C++ are:
 1. functions may be prototyped in  high-level code, using
    extensive mathematical libraries, and no boilerplate;
 2. just-in-time compilation
    occurs transparently: there are no compilation and link delays
    for code changes;
 3. CVs may be automatically differentiated in common cases.

Please see the project's
[homepage](https://giorginolab.github.io/plumed2-pycv/),
[paper](https://doi.org/10.21105/joss.01773), and [regression
tests](https://github.com/giorginolab/plumed2-pycv/tree/v2.6-pycv-devel/regtest/pycv)
for detailed instructions.


Last tested with
------------------------------

Software | Version
---------|---------
PLUMED | 2.6.0
Python | 3.7.7 (Homebrew)
Jax | 0.1.68
Jaxlib | 0.1.47
Numpy | 1.17.0





Citation
------------------------------------

Giorgino T. PYCV: a PLUMED 2 Module Enabling the Rapid Prototyping of
Collective Variables in Python. The Journal of Open Source Software
4(42):1773

[![DOI](https://joss.theoj.org/papers/10.21105/joss.01773/status.svg)](https://doi.org/10.21105/joss.01773) 
[![plumID:19.075](https://www.plumed-nest.org/eggs/19/075/badge.svg)](https://www.plumed-nest.org/eggs/19/075/)


Copyright
------------------------------------

PYCV is distributed under the LGPL terms: see COPYRIGHT.

PYCV includes the [pybind11](https://github.com/pybind/pybind11)
library, distributed under its own license terms (BSD-3).

