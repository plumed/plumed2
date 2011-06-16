#ifndef __PLUMED_Pbc_h
#define __PLUMED_Pbc_h
#include <cassert>
#include "Vector.h"
#include "Tensor.h"

namespace PLMD{

class Pbc{
  enum {unset,orthorombic,xy,xz,yz,generic} type;
  Tensor box;
  Tensor invBox;
public:
  Pbc();
  Vector distance(const Vector&,const Vector&)const;
  void setBox(const Tensor&);
};

}

#endif
