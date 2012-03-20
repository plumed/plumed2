#ifndef __PLUMED_Pbc_h
#define __PLUMED_Pbc_h

#include "Vector.h"
#include "Tensor.h"

namespace PLMD{

class Pbc{
  enum {unset,orthorombic,xy,xz,yz,generic} type;
  Tensor box;
  Tensor invBox;
public:
  Pbc();
  double distance( const bool pbc, const Vector& v1, const Vector& v2 ) const;
  Vector distance(const Vector&,const Vector&)const;
  void setBox(const Tensor&);
  const Tensor& getBox()const;
  const Tensor& getInvBox()const;
  Vector realToScaled(const Vector&)const;
  Vector scaledToReal(const Vector&)const;
  bool isOrthorombic()const;
};

}

#endif
