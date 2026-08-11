Chemical systems contain an enormous number atoms, which, in most cases makes it simply impossible for
us to understand anything by monitoring the atom positions directly.  Consequently,
we introduce Collective variables (CVs) that describe the chemical processes we are
interested in and monitor these simpler quantities instead.  These CVs are used in many of the methods
implemented in PLUMED - there values can be monitored using [PRINT](PRINT.md), Functions of them can be calculated using the methods in the 
[function](module_function.md) module, they can be analyzed using the tools in the [dimred](module_dimred.md), [gridtools](module_gridtools.md)
or [generic](module_generic.md) moeuls or they can be biased using the tools in the [bias](module_bias.md) module.  Before doing any of these things, 
however, we first have to tell PLUMED how to calculate them.

The simplest collective variables that are implemented in PLUMED take in a
set of atomic positions and output one or multiple scalar CV values.  Many of the variables that operate like this are provided in this 
module.  

Please be aware that many other modules contain implementations other collective variables.  In other words, the colvar module does not 
contain implementations of all the collectivar variables that are available in PLUMED. 

## Reactive soft-Voronoi collective variables

Proton-transfer and acid-base reactions can change molecular identities, so a
fixed list of bonds or permanent ions is often not a suitable reaction
coordinate.  The reactive soft-Voronoi family separates one shared smooth
assignment from three distinct physical reductions:

- [VORONOI_COORDINATION](VORONOI_COORDINATION.md) measures signed or squared
  coordination-defect activity with explicit references, selections, and
  coefficients;
- [VORONOI_DISTANCE](VORONOI_DISTANCE.md) combines defects with explicit
  within-group or cross-group center distances;
- [VORONOI_POSITION](VORONOI_POSITION.md) resolves selected defects along a
  Cartesian direction relative to a declared fixed origin.

The Actions are part of the default `colvar` module and require no external
library.  For development they can also be compiled as one runtime plugin with
`plumed mklib ReactiveVoronoi.cpp` and loaded with
`LOAD FILE=./ReactiveVoronoi.so`, avoiding a full PLUMED rebuild.  Only examples
that use OPES require the separately enabled `opes` module.

Start with the VORONOI_COORDINATION page.  It contains the common mathematics,
installation instructions, chemistry-to-keyword workflow, exact versus NLIST
guidance, OpenMP/MPI scaling and GPU-host CPU allocation guidance, derivative
and biasing cautions, validation checklist, and troubleshooting.  The distance
and position pages add detailed worked inputs for water autoionization, single
reactive O/N sites, solvated glycine, nitrogen reduction, and air/oil-water
interfaces.
