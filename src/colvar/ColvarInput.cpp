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
#include "ColvarInput.h"
#include "core/Colvar.h"

namespace PLMD {
namespace colvar {

ColvarInput ColvarInput::createColvarInput( unsigned m,
    const std::vector<Vector>& p,
    const Colvar* colv ) {
  return ColvarInput( m,
                      p.size(),
                      &p[0][0],
                      colv->getMasses().data(),
                      colv->getCharges(true).data(),
                      colv->getPbc() );
}

void ColvarInput::setBoxDerivativesNoPbc( const ColvarInput& inpt, ColvarOutput& out ) {
  unsigned nat=inpt.pos.size();
  for(unsigned i=0; i<out.ncomponents; ++i) {
      Tensor v; v.zero();
      for(unsigned j=0; j<nat; j++) {
        v-=Tensor(Vector(inpt.pos[j][0],inpt.pos[j][1],inpt.pos[j][2]),
                  out.derivs.getAtomDerivatives(i,j));
      }           
      out.virial.set( i, v );
    }
}

} // namespace colvar
} // namespace PLMD
