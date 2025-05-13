This Experiment Directed Simulation module contains methods for adaptively determining linear bias parameters such that each biased CV samples a new target mean value. 
This module implements the stochastic gradient descent algorithm in the original EDS paper that is cited below as well as additional minimization algorithms for Coarse-Grained Directed Simulation 
that are discussed in the second paper cited below.

The third paper cited below is a recent review on the method and its applications. 

Notice that a similar method is available as [MAXENT](MAXENT.md), although with different features and using a different optimization algorithm.

A tutorial using EDS specifically for biasing coordination number can be found on <a href="http://thewhitelab.org/Blog/tutorial/2017/05/10/lammps-coordination-number-tutorial/">Andrew White's webpage</a>.

## Installation 

This module is not installed by default. Add `--enable-modules=eds` to your './configure' command when building PLUMED to enable these features.

