This module implements PINES (Permutationally Invariant Networks for Enhanced Sampling), a permutationally invariant featurization of the atomic configuration that is discussed in the paper cited below.
The [PINES](PINES.md) action computes a vector of switching functions of selected interatomic distances, sorted within blocks of atom-group pairs in the spirit of the [PIV](PIV.md) collective variable.
The elements of this vector are intended to be used as input features of a neural network that learns the collective variables (for example with [PYTORCH_MODEL](PYTORCH_MODEL.md)), in systems such as solutions where molecules of the same species are interchangeable.

## Installation

This module is not installed by default. Add `--enable-modules=pines` to your './configure' command when building PLUMED to enable these features.
