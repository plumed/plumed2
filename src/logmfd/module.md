The LOGMFD module contains the LogMFD/LogPD method for enhanced sampling in a CV space and for on-the-fly free energy reconstruction along the CVs. 
This module implements the multiple-replica algorithm (see second reference below) as well as the single-replica algorithm (see first reference) the former invoking the Crooks-Jarzynski non-equilibrium work relation. In addition, TAMD/d-AFED, which is discussed in the fourth of the papers cited below, can also be implemented by this module.

## Installation 

This module is not installed by default. Add `--enable-modules=logmfd` to your './configure' command when building PLUMED to enable these features.

