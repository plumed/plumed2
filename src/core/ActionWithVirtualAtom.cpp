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
}

ActionWithVirtualAtom::ActionWithVirtualAtom(const ActionOptions&ao):
  Action(ao),
  ActionAtomistic(ao),
  boxDerivatives(3)
{
  index=atoms.addVirtualAtom(this);
  log.printf("  serial associated to this virtual atom is %u\n",index.serial());
}

ActionWithVirtualAtom::~ActionWithVirtualAtom() {
  atoms.removeVirtualAtom(this);
}

void ActionWithVirtualAtom::apply() {
  Vector & f(atoms.forces[index.index()]);
  for(unsigned i=0; i<getNumberOfAtoms(); i++) modifyForces()[i]=matmul(derivatives[i],f);
  Tensor & v(modifyVirial());
  for(unsigned i=0; i<3; i++) v+=boxDerivatives[i]*f[i];
  f.zero(); // after propagating the force to the atoms used to compute the vatom, we reset this to zero
  // this is necessary to avoid double counting if then one tries to compute the total force on the c.o.m. of the system.
  // notice that this is currently done in FIT_TO_TEMPLATE
}

void ActionWithVirtualAtom::requestAtoms(const std::vector<AtomNumber> & a) {
  ActionAtomistic::requestAtoms(a);
  derivatives.resize(a.size());
}

void ActionWithVirtualAtom::setGradients() {
  gradients.clear();
  for(unsigned i=0; i<getNumberOfAtoms(); i++) {
    AtomNumber an=getAbsoluteIndex(i);
    // this case if the atom is a virtual one
    if(atoms.isVirtualAtom(an)) {
      const ActionWithVirtualAtom* a=atoms.getVirtualAtomsAction(an);
      for(const auto & p : a->gradients) {
        gradients[p.first]+=matmul(derivatives[i],p.second);
      }
      // this case if the atom is a normal one
    } else {
      gradients[an]+=derivatives[i];
    }
  }
}

void ActionWithVirtualAtom::setBoxDerivatives(const std::vector<Tensor> &d) {
  boxDerivatives=d;
  plumed_assert(d.size()==3);
// Subtract the trivial part coming from a distorsion applied to the ghost atom first.
// Notice that this part alone should exactly cancel the already accumulated virial
// due to forces on this atom.
  Vector pos=atoms.positions[index.index()];
  for(unsigned i=0; i<3; i++) for(unsigned j=0; j<3; j++) boxDerivatives[j][i][j]+=pos[i];
}

void ActionWithVirtualAtom::setBoxDerivativesNoPbc() {
  std::vector<Tensor> bd(3);
  for(unsigned i=0; i<3; i++) for(unsigned j=0; j<3; j++) for(unsigned k=0; k<3; k++) {
// Notice that this expression is very similar to the one used in Colvar::setBoxDerivativesNoPbc().
// Indeed, we have the negative of a sum over dependent atoms (l) of the external product between positions
// and derivatives. Notice that this only works only when Pbc have not been used to compute
// derivatives.
        for(unsigned l=0; l<getNumberOfAtoms(); l++) {
          bd[k][i][j]-=getPosition(l)[i]*derivatives[l][j][k];
        }
      }
  setBoxDerivatives(bd);
}



void ActionWithVirtualAtom::setGradientsIfNeeded() {
  if(isOptionOn("GRADIENTS")) {
    setGradients() ;
  }
}

}
