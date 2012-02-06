#include "Vector.h"
#include "PlumedException.h"
#include <cmath>

namespace PLMD{

/// Small auxiliary class.
/// I use it to test a few things that I am scary of and could introduce bugs.
/// It checks at startup that Vector satifies some requirement so as to allow
/// accessing a vector of tensors as a 3 times longer array of doubles.
static class VectorChecks{
public:
  VectorChecks(){
    if( sizeof(VectorGeneric<2>)==2*sizeof(double)
     && sizeof(VectorGeneric<3>)==3*sizeof(double)
     && sizeof(VectorGeneric<4>)==4*sizeof(double)) return;
    plumed_merror("sizeof(VectorGeneric<x>)!=x*sizeof(double). PLUMED cannot work properly in these conditions.");
  }
} checks;

}



