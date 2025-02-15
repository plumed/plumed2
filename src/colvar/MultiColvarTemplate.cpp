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
  mass(natoms,w),
  charges(natoms,q)
{
}

ColvarInput ColvarInput::createColvarInput( unsigned m, const std::vector<Vector>& p, const Colvar* colv ) {
  return ColvarInput( m, p.size(), &p[0][0], colv->getMasses().data(), colv->getCharges(true).data(), colv->getPbc() );
}

ColvarOutput::ColvarOutput( std::vector<double>& v, Matrix<Vector>& d, std::vector<Tensor>& t ):
  values(v),
  virial(t),
  derivs(d)
{
}

ColvarOutput ColvarOutput::createColvarOutput( std::vector<double>& v, Matrix<Vector>& d, std::vector<Tensor>& t ) {
  return ColvarOutput( v, d, t );
}

void ColvarOutput::setBoxDerivativesNoPbc( const ColvarInput& inpt ) {
  unsigned nat=inpt.pos.size();
  for(unsigned i=0; i<virial.size(); ++i) {
    virial[i].zero(); for(unsigned j=0; j<nat; j++) virial[i]-=Tensor(Vector(inpt.pos[j][0],inpt.pos[j][1],inpt.pos[j][2]),derivs[i][j]);
  }
}

}
}
