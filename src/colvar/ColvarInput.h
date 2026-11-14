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
  // the virial is accumulated in local variables: building a Tensor for each atom
  // and subtracting it made the loop stall on the store-to-load forwarding of the temporary.
  // The sums are in double, as before, also when T is float
  const unsigned nat=inpt.pos.size();
  for(unsigned i=0; i<out.ncomponents; ++i) {
    double v00=0, v01=0, v02=0, v10=0, v11=0, v12=0, v20=0, v21=0, v22=0;
    for(unsigned j=0; j<nat; j++) {
      const auto deriv=out.derivs.getView(i,j);
      const double x=inpt.pos[j][0], y=inpt.pos[j][1], z=inpt.pos[j][2];
      v00-=x*deriv[0];
      v01-=x*deriv[1];
      v02-=x*deriv[2];
      v10-=y*deriv[0];
      v11-=y*deriv[1];
      v12-=y*deriv[2];
      v20-=z*deriv[0];
      v21-=z*deriv[1];
      v22-=z*deriv[2];
    }
    out.virial.set( i, Tensor(v00,v01,v02,v10,v11,v12,v20,v21,v22) );
  }
}

} // namespace colvar
} //namespace PLMD
#endif // __PLUMED_colvar_ColvarInput_h
