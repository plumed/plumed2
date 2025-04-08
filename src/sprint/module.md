This module contains an implementation of the SPRINT method that was introduced by Pietrucci and Andreoni in the paper that is cited below.
The implementation of this method that has been incldued in PLUMED is a shortcut method that uses features from the [adjmat](module_adjmat.md)
and [matrixtools](module_matrixtools.md) modules.  By looking through these shortcuts you can see how the SPRINT CVs are calculated and how new 
SPRINT-like CVs based on graph theory might be introduced by using different definitions for the adjacency matrix.  
