/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2025 The plumed team
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
#ifndef __PLUMED_colvar_ColvarInput_h
#define __PLUMED_colvar_ColvarInput_h

#include <vector>

#include "tools/Pbc.h"
#include "tools/View.h"
#include "tools/View2D.h"
#include "tools/Vector.h"
#include "tools/Tensor.h"
#include "tools/ColvarOutput.h"

namespace PLMD {

class Colvar;

namespace colvar {

template <typename T=double>
struct ColvarInput {
  unsigned mode;
  const Pbc& pbc;
  View2D<const T,helpers::dynamic_extent,3> pos;
  View<const T> mass;
  View<const T> charges;
  ColvarInput( unsigned m,
               unsigned natoms,
               const T* p,
               const T* w,
               const T* q,
               const Pbc& box ) :
    mode(m),
    pbc(box),
    pos(p,natoms),
    mass(w,natoms),
    charges(q,natoms) {
  }
  static ColvarInput createColvarInput( unsigned m, const std::vector<Vector>& p, const Colvar* colv );
  static void setBoxDerivativesNoPbc( const ColvarInput<T>& inpt, ColvarOutput<T>& out );
  /// same as setBoxDerivativesNoPbc but with no extra memory allocations
  static void setBoxDerivativesNoPbc_inplace( const ColvarInput<T>& inpt, ColvarOutput<T>& out );
};

template <typename T>
ColvarInput<T>  ColvarInput<T>::createColvarInput( unsigned m,
    const std::vector<Vector>& p,
    const Colvar* colv ) {
  //this won't work with floating point
  return ColvarInput( m,
                      p.size(),
                      &p[0][0],
                      colv->getMasses().data(),
                      colv->getCharges(true).data(),
                      colv->getPbc() );
}

template <typename T>
void  ColvarInput<T>::setBoxDerivativesNoPbc( const ColvarInput<T>& inpt, ColvarOutput<T>& out ) {
  //Both version passes the tests, we should discuss wicht onw might be better

  unsigned nat=inpt.pos.size();
  for(unsigned i=0; i<out.ncomponents; ++i) {
    Tensor v;
    v.zero();
    for(unsigned j=0; j<nat; j++) {
      v-=Tensor(Vector(inpt.pos[j][0],inpt.pos[j][1],inpt.pos[j][2]),
                out.derivs.getAtomDerivatives(i,j));
    }
    out.virial.set( i, v );
  }
}

template <typename T>
void  ColvarInput<T>::setBoxDerivativesNoPbc_inplace( const ColvarInput<T>& inpt, ColvarOutput<T>& out ) {
  //now with no extra allocated memory:
  unsigned nat=inpt.pos.size();
  for(unsigned i=0; i<out.ncomponents; ++i) {
    auto v = out.virial.getView(i);
    LoopUnroller<9>::_zero(v.data());
    for(unsigned j=0; j<nat; j++) {
      const auto deriv =  out.derivs.getView(i,j);
      v[0] -= inpt.pos[j][0]*deriv[0];
      v[1] -= inpt.pos[j][0]*deriv[1];
      v[2] -= inpt.pos[j][0]*deriv[2];

      v[3] -= inpt.pos[j][1]*deriv[0];
      v[4] -= inpt.pos[j][1]*deriv[1];
      v[5] -= inpt.pos[j][1]*deriv[2];

      v[6] -= inpt.pos[j][2]*deriv[0];
      v[7] -= inpt.pos[j][2]*deriv[1];
      v[8] -= inpt.pos[j][2]*deriv[2];
    }
  }
}
} // namespace colvar
} //namespace PLMD
#endif // __PLUMED_colvar_ColvarInput_h
