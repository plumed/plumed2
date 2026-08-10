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
coordinate.  [VORONOI_COORDINATION](VORONOI_COORDINATION.md) assigns
transferable atoms smoothly to user-selected centers and measures their
coordination defects.  [VORONOI_DISTANCE](VORONOI_DISTANCE.md) converts these
defects into a defect-pair separation, while
[VORONOI_POSITION](VORONOI_POSITION.md) resolves them along a Cartesian
direction relative to a fixed origin.

These Actions are part of the default `colvar` module.  Their Manual
pages explain the mathematical definition, keyword mapping, neighbor-list
approximation, limitations, and worked examples for water autoionization,
solvated glycine, interfaces, and catalytic proton transfer.
