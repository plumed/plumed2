/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2020 The plumed team
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
#include "ActionWithVirtualAtom.h"
#include "Atoms.h"

namespace PLMD {

void ActionWithVirtualAtom::registerKeywords(Keywords& keys) {
  Action::registerKeywords(keys);
  ActionAtomistic::registerKeywords(keys);
  keys.add("atoms","ATOMS","the list of atoms which are involved the virtual atom's definition");
  keys.addOutputComponent("x","default","the x coordinate of the virtual atom");
  keys.addOutputComponent("y","default","the y coordinate of the virtual atom");
  keys.addOutputComponent("z","default","the z coordinate of the virtual atom");
  keys.addOutputComponent("mass","default","the mass of the virtual atom");
  keys.addOutputComponent("charge","default","the charge of the virtual atom");
}

ActionWithVirtualAtom::ActionWithVirtualAtom(const ActionOptions&ao):
  Action(ao),
  ActionAtomistic(ao),
  ActionWithValue(ao),
  forcesToApply(9)
{
  addComponentWithDerivatives("x"); componentIsNotPeriodic("x");
  addComponentWithDerivatives("y"); componentIsNotPeriodic("y");
  addComponentWithDerivatives("z"); componentIsNotPeriodic("z");
  // Store the derivatives with respect to the virial only even if there are no atoms
  for(unsigned i=0; i<3; ++i) getPntrToComponent(i)->resizeDerivatives(9);
  addComponent("mass"); componentIsNotPeriodic("mass"); addComponent("charge"); componentIsNotPeriodic("charge");
  atoms.addAtomValues( getLabel(), copyOutput(0), copyOutput(1), copyOutput(2), copyOutput(3), copyOutput(4) );
}

void ActionWithVirtualAtom::apply() {
  std::fill(forcesToApply.begin(),forcesToApply.end(),0); unsigned mm=0;
  if( getForcesFromValues( forcesToApply ) ) setForcesOnAtoms( forcesToApply, mm );
}

void ActionWithVirtualAtom::requestAtoms(const std::vector<AtomNumber> & a) {
  ActionAtomistic::requestAtoms(a); forcesToApply.resize(3*a.size()+9);
  for(unsigned i=0; i<3; ++i) getPntrToComponent(i)->resizeDerivatives(3*a.size()+9);
}

void ActionWithVirtualAtom::setBoxDerivatives(const std::vector<Tensor> &d) {
  plumed_assert(d.size()==3); unsigned nbase = 3*getNumberOfAtoms();
  for(unsigned i=0; i<3; ++i) {
      Value* myval = getPntrToComponent(i);
      for(unsigned j=0; j<3; ++j) {
          for(unsigned k=0; k<3; ++k) myval->setDerivative( nbase + 3*j + k, d[i][j][k] );
      }
  }
// Subtract the trivial part coming from a distorsion applied to the ghost atom first.
// Notice that this part alone should exactly cancel the already accumulated virial
// due to forces on this atom.
  Vector pos; for(unsigned i=0; i<3; ++i) pos[i] = getPntrToComponent(i)->get();
  for(unsigned i=0; i<3; i++) for(unsigned j=0; j<3; j++) getPntrToComponent(j)->addDerivative( nbase + 3*i + j, pos[i] );
}

void ActionWithVirtualAtom::setBoxDerivativesNoPbc() {
  std::vector<Tensor> bd(3);
  for(unsigned i=0; i<3; i++) for(unsigned j=0; j<3; j++) for(unsigned k=0; k<3; k++) {
// Notice that this expression is very similar to the one used in Colvar::setBoxDerivativesNoPbc().
// Indeed, we have the negative of a sum over dependent atoms (l) of the external product between positions
// and derivatives. Notice that this only works only when Pbc have not been used to compute
// derivatives.
        for(unsigned l=0; l<getNumberOfAtoms(); l++) {
          bd[k][i][j]-=getPosition(l)[i]*getPntrToComponent(k)->getDerivative(3*l+j);
        }
      }
  setBoxDerivatives(bd);
}

}
