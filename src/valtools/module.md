Each line in a PLUMED input file defines and action. The actions in the input file communicate through the passing of values.
In early versions of PLUMED these values were scalars.  However, from PLUMED 2.10 onwards you can pass values that are scalars, 
vectors, matrices or functions on grids.  

The tools in this module allow you to convert values from one type to another.  For example, [CONCATENATE](CONCATENATE.md) allows you 
to combine multiple scalars and vectors into a single larger vector, [FLATTEN](FLATTEN.md) allows you to convert a matrix to a vector 
and [VSTACK](VSTACK.md) allows you to combine multiple input vectors into a matrix. 

