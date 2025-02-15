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

ColvarInput::ColvarInput( unsigned m, unsigned natoms, const std::vector<Vector>& p, const double* w, const double* q, const Pbc& box ) :
  mode(m),
  pbc(box),
  pos(p),
  mass(natoms,w),
  charges(natoms,q)
{
}

ColvarInput ColvarInput::createColvarInput( unsigned m, const std::vector<Vector>& p, const Colvar* colv ) {
  return ColvarInput( m, p.size(), p, colv->getMasses().data(), colv->getCharges(true).data(), colv->getPbc() );
}

}
}
