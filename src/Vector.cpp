#include "Vector.h"
#include <cmath>

using namespace PLMD;

double Vector::modulo()const{
  return sqrt(modulo2());
}

