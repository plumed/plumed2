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
struct ColvarInput {
  unsigned mode;
  const Pbc& pbc;
  View2D<const double,helpers::dynamic_extent,3> pos;
  View<const double> mass;
  View<const double> charges;
  ColvarInput( unsigned m,
               unsigned natoms,
               const double* p,
               const double* w,
               const double* q,
               const Pbc& box ) :
    mode(m),
    pbc(box),
    pos(p,natoms),
    mass(w,natoms),
    charges(q,natoms) {
  }

  static ColvarInput createColvarInput( unsigned m, const std::vector<Vector>& p, const Colvar* colv );
#pragma acc routine seq
  static void setBoxDerivativesNoPbc( const ColvarInput& inpt, ColvarOutput& out );
  /// same as setBoxDerivativesNoPbc but with no extra memory allocations
#pragma acc routine seq
  static void setBoxDerivativesNoPbc_inplace( const ColvarInput& inpt, ColvarOutput& out );
};

} // namespace colvar
} //namespace PLMD
#endif // __PLUMED_colvar_ColvarInput_h
