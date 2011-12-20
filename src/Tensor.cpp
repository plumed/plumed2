#include "Tensor.h"
#include <iostream>
#include <cstdlib>

namespace PLMD{

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

