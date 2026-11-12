This sizeshape module contains methods for calculating a 1D linear projection ([SIZESHAPE_POSITION_LINEAR_PROJ](SIZESHAPE_POSITION_LINEAR_PROJ.md)) and the Mahalanobis distance ([SIZESHAPE_POSITION_MAHA_DIST](SIZESHAPE_POSITION_MAHA_DIST.md)) in size-and-shape space for a given reference configurational distribution.

The reference structure and the precision matrix that define this distribution are typically obtained by clustering an ensemble of configurations
with a size-and-shape space Gaussian mixture model (shapeGMM, see the second paper cited below).
The coefficients of the linear projection can be obtained, for example, by linear discriminant analysis on the aligned positions of two
states, which provides a reaction coordinate for conformational transitions that can be biased (see the first paper cited below).

## Installation 

This module is not installed by default. Add `--enable-modules=sizeshape` to your './configure' command when building PLUMED to enable these features.

