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
#include "Colvar.h"
#include "tools/OpenMP.h"
#include <vector>

namespace PLMD {

Colvar::Colvar(const ActionOptions&ao):
  Action(ao),
  ActionAtomistic(ao),
  ActionWithValue(ao) {
}

void Colvar::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys );
  keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions when calculating distances");
}

void Colvar::requestAtoms(const std::vector<AtomNumber> & a) {
// Tell actionAtomistic what atoms we are getting
  ActionAtomistic::requestAtoms(a);
// Resize the derivatives of all atoms
  for(unsigned i=0; i<getNumberOfComponents(); ++i) {
    getPntrToComponent(i)->resizeDerivatives(3*a.size()+9);
  }
}

void Colvar::apply() {
  if( !checkForForces() ) {
    return ;
  }
  unsigned ind=0;
  if( getNumberOfAtoms()>0 ) {
    setForcesOnAtoms( getForcesToApply(), ind );
  } else {
    setForcesOnCell( getForcesToApply(), ind );
  }
}

void Colvar::setBoxDerivativesNoPbc(Value* v) {
  Tensor virial;
  unsigned nat=getNumberOfAtoms();
  for(unsigned i=0; i<nat; i++)
    virial-=Tensor(getPosition(i),
                   Vector(v->getDerivative(3*i+0),
                          v->getDerivative(3*i+1),
                          v->getDerivative(3*i+2)));
  setBoxDerivatives(v,virial);
}
}
