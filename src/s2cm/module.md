This module implements the S2 contact model CV ([S2CM](S2CM.md)) that is discussed in the paper cited below.
The CV estimates the NMR order parameter $S^2$ of a backbone NH group or of a side-chain methyl group from the number of heavy atoms that surround it. Summing it over many residues gives a proxy of the conformational entropy of a protein that can be used as a collective variable, as done in the paper cited below.

## Installation 

This module is not installed by default. Add `--enable-modules=s2cm` to your './configure' command when building PLUMED to enable these features.

