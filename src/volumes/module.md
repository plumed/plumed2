The actions in this module take a vector of $n$ atomic positions in input.  These actions then return an $n$-dimensional vector. The elements
of this vector are one if the correpsonding input atom is in a particular part of the box and zero otherwise.  These actions were implemented so that 
you can calculate functions with the following functional form:

$$
s = \sum_i w(\{X\}_i) g[f(\{X\}_i)]
$$

In this expression, $g$ is a function with one argument and $g$ is a function of a set of atomic positions, $\{X\}_i$.  
Meanwhile, \f$w\f$ is a weight that is also a function of the set of atomic positions. This weight varies between zero and one.  You can 
thus use expressions like the one above to determine the average value of a quanitty in a particular part of the simulation box.
