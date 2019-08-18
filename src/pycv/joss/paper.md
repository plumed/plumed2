---
title: 'PYCV: a PLUMED 2 module enabling Rapid Development of Python-based Collective Variables'
tags:
  - Python
  - PLUMED
  - molecular dynamics
  - collective variables
  - biased sampling
  - gradient
authors:
  - name: Toni Giorgino
    orcid: 0000-0001-6449-0596
    affiliation: 1
affiliations:
 - name: Institute of Biophysics (IBF-CNR), National Research Council of Italy
   index: 1
date: 15 August 2019
bibliography: paper.bib
---

# Summary

Collective variables (CVs) are functions of the coordinates of
particles in a molecular system. The choice of CV is
crucial to capture relevant degrees of freedom of the model being
simulated [@barducci_metadynamics_2011]. This is especially important
when employing *biased sampling* techniques such as umbrella sampling
or metadynamics [@torrie_nonphysical_1977, @laio_escaping_2002], which
apply generalized forces to CVs to enhance the sampling of events
otherwise not observable by direct simulation.   CVs may be
simple geometrical observables (distances, angles, torsions, etc.),
but are often more complex functions designed to capture structural
determinants, such as tertiary and quaternary structure of proteins,
experimental observables, crystal symmetries, etc. [@branduardi_b_2007;
@bonomi_integrative_2017; @pipolo_navigating_2017].

Iterative development of CVs therefore occupies a large fraction of
the effort in the exploration of molecular systems. On the one hand,
the task has been largely facilitated by biasing libraries such as
PLUMED [@tribello_plumed_2014], which provides pre-defined functions and a
*lingua franca* to express CV combinations, atom groups and biasing
schemes, and its active community [@the_plumed_consortium_promoting_2019].
However, users willing to explore CVs beyond the pre-defined
ones have to implement them in C++, together with the corresponding
(often cumbersome) derivatives [@giorgino_how_2018]. Compiled code is
unwieldy for iterative analysis, because it is relatively low-level,
error-prone, and inconvenient in exploratory stages.

This paper introduces **PYCV**, a module for the PLUMED 2 library
which enables users to define CVs and arbitrary functions in the
Python language.  CV implementations may thus be modified and tested
independently of the main code, with essentially no "test latency".
Of note, coordinates are processed as `numpy` arrays, making it
convenient to leverage the vast set of linear algebra and numerical
algorithms provided by `numpy`, `scipy`, and many other open-source
modules. Furthermore, just-in-time compilation and reverse-mode
automatic differentiation are easily accessible using Google's JAX
library.


# Usage

The **PYCV** module registers itself with PLUMED to provide the
following actions:

 * `PYTHONCV`, to implement single- and multi-component CVs;
 * `PYTHONFUNCTION`, to implement arbitrary functions.

The actions are documented in the respective inline manuals (e.g.,
`plumed manual --action PYTHONCV`).  In both cases, an interpreter is
first started; the Python module indicated in the `IMPORT=` keyword is
then loaded; from it, an user-chosen function (`FUNCTION=`) is called
to perform the computations at each timestep. Modules can contain
multiple functions and one-time initialization.



# Example

A self-explanatory example is provided for illustration below. It is
equivalent to the *radius of curvature* example shown in
[@giorgino_how_2018]. Further examples are available in the manual and
in `regtest/pycv`.


## Biasing script

The actions are declared in the PLUMED input file (say,
`plumed.dat`). Here, one declares a CV labelled `rc`, to be computed by
the Python function `curvature.r()`. It will receive a 3-by-3 array
with the coordinates of atoms 1, 4 and 3 (orderly, as rows).  The CV
value will be `PRINT`ed, and the atoms subject to a constant generalized
force pushing to increase the curvature.

```
# Start plumed.dat -----------------------------------------
rc: PYTHONCV ATOMS=1,4,3 IMPORT=curvature FUNCTION=r
    PRINT ARG=rc FILE=colvar.out
    RESTRAINT ARG=rc AT=0 KAPPA=0 SLOPE=1
# End plumed.dat -------------------------------------------
```


## Function definition

The actual function `r` is defined in the `curvature.py` file. It
computes the radius of the circle passing through three given atom
coordinates (the three rows of the input argument, with 0-based
indexing). Note how matrix operations enable a readable translation of
the sine formula.

The function is expected to return two objects, i.e. the value of the
CV at the given coordinates (a scalar), and its gradient with respect
to each of the 9 coordinates (a 3-by-3 array); here the gradient
function is obtained automatically. Although not directly comparable,
the equivalent C++ implementation required approximately 160 lines of
non-trivial code.


```py
# Start curvature.py ----------------------------------------
# Same calculation (and results) as doi:10.1016/j.cpc.2018.02.017

# Import the JAX library
import jax.numpy as np
from jax import grad, jit

# Implementation of the angle function. @jit really improves speed
@jit
def r_f(x):
    r21 = x[0,:]-x[1,:]
    r23 = x[2,:]-x[1,:]
    r13 = x[2,:]-x[0,:]

    cos2theta = np.dot(r21,r23)**2 / (np.dot(r21,r21) * np.dot(r23,r23))
    sin2theta = 1-cos2theta
    
    R2= np.dot(r13,r13)/sin2theta/4.0
    return np.sqrt(R2)

# Use JAX to auto-gradient it
r_g = grad(r_f)

# The CV function actually called
def r(x):
    return r_f(x), r_g(x)

# End curvature.py ---------------------------------------------
```



# Conclusion

**PYCV** enables Python-based prototyping of CVs in PLUMED 2. This
programming model may be an advantage over standard C++-based development in that
(a) functions may be prototyped in high-level language, using extensive
mathematical libraries, without boilerplate; (b) just-in-time
compilation occurs transparently: code changes incur in no compilation
and link delays; and (c) CV code may be automatically differentiated in
common cases.




# Acknowledgements

I would like to thank PLUMED's lead authors (Prof. Giovanni Bussi,
Prof. Carlo Camilloni, Prof. Gareth Tribello and Prof. Massimiliano
Bonomi) for inputs and discussions.

# References