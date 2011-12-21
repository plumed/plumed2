#include "Vector.h"
#include <cmath>
#include <iostream>
#include <cstdlib>

namespace PLMD{

double Vector::modulo()const{
  return sqrt(modulo2());
}

/// Small auxiliary class.
/// I use it to test a few things that I am scary of and could introduce bugs.
/// It checks at startup that Vector satifies some requirement so as to allow
/// accessing a vector of tensors as a 3 times longer array of doubles.
static class VectorChecks{
public:
  VectorChecks(){
    if(sizeof(Vector)==3*sizeof(double)) return;
    std::cerr<<"Severe error: sizeof(Vector)!=3*sizeof(double)\n";
    std::cerr<<"PLUMED cannot work properly in these conditions\n";
    std::abort();
  }
} checks;

}



