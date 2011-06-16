/*
  We here compile a version of Plumed.c which is statically bound to the plumed library,
  for which it is not possible to redefine the location of plumed kernel at runtime.
  Since Plumed.c can be compiled with both C and C++, we opt here for C++
  so as to further check compatibilities
*/
#define __PLUMED_STATIC_KERNEL
#include "Plumed.c"
