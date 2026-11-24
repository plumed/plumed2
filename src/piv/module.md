This module implements the Permutation Invariant Vector (PIV) collective variable ([PIV](PIV.md)) that is discussed in the first paper cited below.
The PIV of a configuration is built by computing all the distances between atoms of selected types, transforming them with switching functions and sorting them within each block of atom-type pairs.
The resulting vector is invariant under the permutation of identical atoms, which makes it suitable for describing processes such as nucleation, phase transitions or chemical reactions in which atoms exchange their roles.
The PIV-distance between the current configuration and one or more reference configurations can be used directly as a CV or combined with [FUNCPATHMSD](FUNCPATHMSD.md) to define a path CV in PIV space.

## Installation 

This module is not installed by default. Add `--enable-modules=piv` to your './configure' command when building PLUMED to enable these features.
