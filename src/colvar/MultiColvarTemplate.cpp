/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2023 The plumed team
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
#include "MultiColvarTemplate.h"
#include "core/Colvar.h"

namespace PLMD {
namespace colvar {

ColvarInput::ColvarInput( unsigned m, unsigned natoms, const double* p, const double* w, const double* q, const Pbc& box ) :
  mode(m),
  pbc(box),
  pos(p,natoms),
  mass(w,natoms),
  charges(q,natoms)
{
}

ColvarInput ColvarInput::createColvarInput( unsigned m, const std::vector<Vector>& p, const Colvar* colv ) {
  return ColvarInput( m, p.size(), &p[0][0], colv->getMasses().data(), colv->getCharges(true).data(), colv->getPbc() );
}

ColvarOutput::ColvarOutput( View<double,helpers::dynamic_extent>& v, std::size_t m, std::vector<double>& d ):
  ncomponents(v.size()),
  values(v),
  derivs(m,d),
  virial(m,d)
{
}

ColvarOutput ColvarOutput::createColvarOutput( std::vector<double>& v, std::vector<double>& d, Colvar* action ) {
  View<double,helpers::dynamic_extent> val(v.data(),v.size());
  d.resize( action->getNumberOfComponents()*action->getNumberOfDerivatives() );
  return ColvarOutput( val, action->getNumberOfDerivatives(), d );
}

void ColvarOutput::setBoxDerivativesNoPbc( const ColvarInput& inpt ) {
  unsigned nat=inpt.pos.size();
  for(unsigned i=0; i<ncomponents; ++i) {
    Tensor v; v.zero();
    for(unsigned j=0; j<nat; j++) v-=Tensor(Vector(inpt.pos[j][0],inpt.pos[j][1],inpt.pos[j][2]),derivs.getAtomDerivatives(i,j));
    virial.set( i, v );
  }
}

}
}
