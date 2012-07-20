/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
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
