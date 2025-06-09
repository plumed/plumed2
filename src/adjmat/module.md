An adjacency matrix is an $N \times N$ matrix in which the $i$th, $j$th element tells you whether or not the $i$th
and $j$th atoms/molecules from a set of $N$ atoms/molecules are adjacent or not.  There are various ways of defining
whether a pair of atoms/molecules are adjacent or not.  For example we can say two atoms are adjacent if the distance between
them is less than some cutoff.  Alternatively, if we have a have a pair of molecules, we might state they are adjacent if their
centers of mass are within a certain cutoff and if the two molecules have the same orientation.  Two electronegative atoms
might be said to be adjacent if there is a hydrogen bond between them.  For these reasons then PLUMED contains all the
multiple methods for calculating adjacency matrices.

Once you have calculated an adjacency matrix you can use it to calculate symmetry functions by exploting functionalility in the 
[symfunc](module_symfunc.md) module.  Alternatively, you can find the connected components in the graph representation of the matrix
by using the functionality in the [clusters](module_clusters.md) module.  You can even work with the properties of this graph representation 
of the matrix directly by exploiting the methods described in the [sprint](module_sprint.md) module.

__You will get a colossal speedup by specifying the D_MAX keyword in all switching functions that act on distances.
D_MAX tells PLUMED that the switching function is strictly zero if the distance is greater than this value.  As a result
PLUMED knows that it does not need to calculate these zero terms in what are essentially sums with a very large number of terms.
In fact when D_MAX is set PLUMED uses linked lists when calculating these coordination numbers, which is what
gives you such a dramatic increase in performance.__

