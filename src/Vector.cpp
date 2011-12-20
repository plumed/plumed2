#include "Vector.h"
#include <cmath>
#include <iostream>
#include <cstdlib>

namespace PLMD{

double Vector::modulo()const{
  return sqrt(modulo2());
}

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



