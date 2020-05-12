/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2020 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

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
#ifndef __PLUMED_multicolvar_CatomPack_h
#define __PLUMED_multicolvar_CatomPack_h

#include <vector>
#include "tools/Exception.h"
#include "tools/Tensor.h"
#include "tools/Vector.h"

namespace PLMD {
namespace multicolvar {

class CatomPack {
private:
  std::vector<unsigned> indices;
  std::vector<Tensor> derivs;
public:
  void resize( const unsigned& );
  void setIndex( const unsigned&, const unsigned& );
  void setDerivative( const unsigned&, const Tensor& );
  unsigned getNumberOfAtomsWithDerivatives() const ;
  unsigned getIndex( const unsigned& ) const ;
  double getDerivative( const unsigned&, const unsigned&, const Vector& ) const ;
};

inline
void CatomPack::setIndex( const unsigned& jind, const unsigned& ind ) {
  plumed_dbg_assert( jind<indices.size() );
  indices[jind]=ind;
}

inline
void CatomPack::setDerivative( const unsigned& jind, const Tensor& der ) {
  plumed_dbg_assert( jind<indices.size() );
  derivs[jind]=der;
}

inline
unsigned CatomPack::getNumberOfAtomsWithDerivatives() const {
  return indices.size();
}

inline
unsigned CatomPack::getIndex( const unsigned& jind ) const {
  plumed_dbg_assert( jind<indices.size() );
  return indices[jind];
}

inline
double CatomPack::getDerivative( const unsigned& iatom, const unsigned& jcomp, const Vector& df ) const {
  plumed_dbg_assert( iatom<indices.size() );
  return df[jcomp]*derivs[iatom](jcomp,0) + df[jcomp]*derivs[iatom](jcomp,1) + df[jcomp]*derivs[iatom](jcomp,2);
}

}
}

#endif
