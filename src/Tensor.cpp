#include "Tensor.h"
#include <iostream>
#include <cstdlib>

namespace PLMD{

/// Small auxiliary class.
/// I use it to test a few things that I am scary of and could introduce bugs.
/// It checks at startup that Tensor satifies some requirement so as to allow
/// accessing a vector of tensors as a 9 times longer array of doubles.
static class TensorChecks{
public:
  TensorChecks(){
    if(sizeof(Tensor)==9*sizeof(double)) return;
    std::cerr<<"Severe error: sizeof(Tensor)!=9*sizeof(double)\n";
    std::cerr<<"PLUMED cannot work properly in these conditions\n";
    std::abort();
  }
} checks;

}

