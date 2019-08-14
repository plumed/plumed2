---
title: 'PYCV: Python-based Rapid Development of Collective Variables for PLUMED 2'
tags:
  - Python
  - PLUMED
  - molecular dynamics
  - collective variables
  - biased sampling
  - gradient
  - JAX
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
particles of a molecular system. They may be simple geometrical
observables (distances, angles, torsions, etc.), but are often complex
functions designed to capture structural determinants, such as
tertiary and quaternary structure of proteins, crystal arrangements,
etc. [@branduardi_b_2007; @bonomi_integrative_2017;
@pipolo_navigating_2017].

The system-specific choice of CV is crucial to capture relevant
degrees of freedom of the model being simulated
[@barducci_metadynamics_2011]. This is especially important when
employing *biased sampling* techniques such as umbrella sampling
[@roux_calculation_1995] or metadynamics [@laio_escaping_2002], which
apply generalized forces to CVs to enhance the sampling of events
otherwise not observable by direct simulation.

Iterative development of CVs therefore occupies a large fraction of
the effort in the exploration of molecular systems. On the one hand,
this task has been largely facilitated by biasing libraries such as
PLUMED [@tribello_plumed_2014], which provide pre-defined functions
and a *lingua franca* to express CV combinations, atom groups and
biasing schemes. However, users willing to explore functions beyond
the pre-defined ones have to implement them, as well as the
corresponding derivatives, in C++ [@giorgino_how_2018]. Compiled code
is however unwieldy for iterative analysis, because it is relatively
low-level, error-prone, and generally inconvenient in the development
phase.

Here, we present **PYCV**, a module for the PLUMED 2 library which
enables users to define CVs and arbitrary functions in Python.  CV
implementations may thus be modified and tested independently of the
main (compiled) code, with essentially no "coding impedance".  Of
note, values are exchanged as `numpy` arrays, making it convenient to
access the vast array of numerical algorithms provided by `numpy`,
`scipy`, and countless other modules. Furthermore, Google's JAX
library may serve as a drop-in replacement for `numpy` to obtain
just-in-time compilation and reverse-mode automatic differentiation
(CV functions are required to return their gradient with respect to
coordinates, which is often cumbersome to implement manually).


# Usage

The **PYCV** module registers itself with PLUMED to provide the
following actions:

 * `PYTHONCV`, to implement single- and multi-component CVs;
 * `PYTHONFUNCTION`, to implement arbitrary functions.

The actions are documented in the respective inline documentation
(accessible e.g. with `plumed manual --action PYTHONCV`).  In both
cases, an interpreter is first started; the Python module indicated in
the `IMPORT=` keyword is then imported; from it, an user-chosen
function (`FUNC=`) is called to perform the computations at each
timestep. Modules can be shared for multiple functions, and contain
one-time initialization.



# Example

Here, a self-explanatory example is provided for illustration
below. Further examples are available in the manual and in
`regtest/pycv`.


## Biasing script

The actions are defined in the PLUMED input file (say,
`plumed.dat`). Here, we define a CV labelled `cv1`, to be computed by
the Python function `jaxcv.cv()`. The CV value
will be printed and subject to a constant generalized force in the
positive direction.

```
# Start plumed.dat -----------------------------------------
cv1:  PYTHONCV ATOMS=1,4,3 IMPORT=jaxcv FUNCTION=cv
      PRINT ARG=cv1 FILE=colvar.out
      RESTRAINT ARG=cv1 AT=0 KAPPA=0 SLOPE=1
# End plumed.dat -------------------------------------------
```


## CV function definition

The actual function is defined in the `jaxcv.py` file. It will compute
the angle (in radians) formed at the second atom by the other two
(rows 1, 0 and 2 of the input array respectively, with 0-based
indexing). Note how matrix operations make for a readable translation
of the cosine formula.


```py
# Start jaxcv.py -------------------------------------------
# Import the JAX library and numpy replacement
import jax.numpy as np
from jax import grad, jit

# Compute the angle between segments r12 and r32. @jit for speed
@jit	      	    	       
def angle(x):
    r12 = x[0,:]-x[1,:]
    r32 = x[2,:]-x[1,:]
    costheta = ( np.dot(r12,r32) /
    	         np.sqrt(np.dot(r12,r12) * np.dot(r32,r32)) )
    theta = np.arccos(costheta)
    return theta

# Use JAX to obtain the corresponding reverse-mode gradient function
grad_angle = grad(angle)

# The CV function actually called
def cv(X):
    return angle(X), grad_angle(X)

# End jaxcv.py ---------------------------------------------
```

The function will receive a 3-by-3 array with the coordinates of atoms
1, 4 and 3 (orderly, as rows). It is expected to return two values,
i.e. the value of the CV (a scalar), and its gradient with respect to
each of the 9 coordinates, computed automatically.


# Conclusion




# Acknowledgements

We acknowledge PLUMED's lead authors (Prof. Giovanni Bussi,
Prof. Carlo Camilloni, Prof. Gareth Tribello and Prof. Massimiliano
Bonomi) for inputs and discussions.

# References