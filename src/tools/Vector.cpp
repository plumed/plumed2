/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2019 The plumed team
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
#include "Vector.h"
#include "Exception.h"
#include <cmath>

namespace PLMD {

/// Small auxiliary class.
/// I use it to test a few things that I am scary of and could introduce bugs.
/// It checks at startup that Vector satifies some requirement so as to allow
/// accessing a vector of tensors as a 3 times longer array of doubles.
static class VectorChecks {
public:
  VectorChecks() {
    if( sizeof(VectorGeneric<2>)==2*sizeof(double)
        && sizeof(VectorGeneric<3>)==3*sizeof(double)
        && sizeof(VectorGeneric<4>)==4*sizeof(double)) return;
    plumed_merror("sizeof(VectorGeneric<x>)!=x*sizeof(double). PLUMED cannot work properly in these conditions.");
  }
} checks;

}



