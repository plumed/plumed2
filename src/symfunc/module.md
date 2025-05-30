The actions in this module can be used to calculate symmetry functions that measure the degree of order in the coordination sphere around
the atoms and molecules that make up your system.  The methods within this class are typically used to investigate nucleation of solids 
from the liquid phase.

Notice that each symmetry function you calculate can be allocated to a particular position in space; namely, the position of the central atom whose 
coordination sphere has been used in the calculation of the symmetry function.  This fact has been used in the papers cited below to:

- calculate average values for the symmetry function in particular parts of the simulation cell by using this module together with the [volumes](module_volumes.md) module.
- calculate average values for the symmetry function in the largest cluster by using this module together with the [clusters](module_clusters.md) module.
- calculate phase field models that provide a smooth representation of the average value of the CV as a function of position by using this module together with the [gridtools](module_gridtools.md) module.
- calculate the gradient of the phase field representation of the CV as a function of position by using the module together with the [gridtools](module_gridtools.md) module.
- calculate the position of the interface between the solid and liquid parts of the simulation cell by using this module together with the together with the [gridtools](module_gridtools.md) and [contour](module_contour.md) modules.

Notice that from version 2.10 onwards many of the symmetry functions in this class (particularly the Steinhardt Parameters) have been implemented as shortcuts. The
implementations that PLUMED provides from these methods are thus very flexible and allow you to calculate the many subtle variationts of this technique for determining the 
degree of order in the system that have been used in the literature.  

__You will get a colossal speedup by specifying the D_MAX keyword in all switching functions that act on distances.
D_MAX tells PLUMED that the switching function is strictly zero if the distance is greater than this value.  As a result
PLUMED knows that it does not need to calculate these zero terms in what are essentially sums with a very large number of terms.
In fact when D_MAX is set PLUMED uses linked lists when calculating these coordination numbers, which is what
gives you such a dramatic increase in performance.__ 
