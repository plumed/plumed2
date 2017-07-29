/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2017 The plumed team
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
#include "core/ActionWithValue.h"
#include "core/ActionAtomistic.h"
#include "MatrixElementPack.h"

namespace PLMD {
namespace adjmat {

MatrixElementPack::MatrixElementPack( MultiValue& vals, const bool comp, ActionWithValue const * myaction ):
myvals(vals),
natoms_in_base(0),
noderiv(myaction->doNotCalculateDerivatives()),
components(comp),
indices(2),
positions(2)
{
  ActionAtomistic const * aa = dynamic_cast<ActionAtomistic const *>( myaction );
  if(aa) natoms_in_base = aa->getNumberOfAtoms();

  if( components ){
      plumed_dbg_assert( myaction->getNumberOfComponents()==4 );
      w_index = (myaction->copyOutput(0))->getPositionInStream();
      x_index = (myaction->copyOutput(1))->getPositionInStream();
      y_index = (myaction->copyOutput(2))->getPositionInStream();
      z_index = (myaction->copyOutput(3))->getPositionInStream();
  } else {
      plumed_dbg_assert( myaction->getNumberOfComponents()==1 );
      w_index = (myaction->copyOutput(0))->getPositionInStream();
  }
}

}
}
