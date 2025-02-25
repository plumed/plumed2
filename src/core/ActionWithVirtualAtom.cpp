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
#include "ActionWithVirtualAtom.h"
#include <array>
#include <iostream>

namespace PLMD {

void ActionWithVirtualAtom::registerKeywords(Keywords& keys) {
  Action::registerKeywords(keys);
  ActionAtomistic::registerKeywords(keys);
  keys.add("atoms","ATOMS","the list of atoms which are involved the virtual atom's definition");
  keys.addOutputComponent("x","default","scalar","the x coordinate of the virtual atom");
  keys.addOutputComponent("y","default","scalar","the y coordinate of the virtual atom");
  keys.addOutputComponent("z","default","scalar","the z coordinate of the virtual atom");
  keys.addOutputComponent("mass","default","scalar","the mass of the virtual atom");
  keys.addOutputComponent("charge","default","scalar","the charge of the virtual atom");
}

ActionWithVirtualAtom::ActionWithVirtualAtom(const ActionOptions&ao):
  Action(ao),
  ActionAtomistic(ao),
  ActionWithValue(ao) {
  addComponentWithDerivatives("x");
  componentIsNotPeriodic("x");
  addComponentWithDerivatives("y");
  componentIsNotPeriodic("y");
  addComponentWithDerivatives("z");
  componentIsNotPeriodic("z");
  // Store the derivatives with respect to the virial only even if there are no atoms
  for(unsigned i=0; i<3; ++i) {
    getPntrToComponent(i)->resizeDerivatives(9);
  }
  addComponent("mass");
  componentIsNotPeriodic("mass");
  getPntrToComponent("mass")->isConstant();
  addComponent("charge");
  componentIsNotPeriodic("charge");
  getPntrToComponent("charge")->isConstant();
}

void ActionWithVirtualAtom::requestAtoms(const std::vector<AtomNumber> & a) {
  ActionAtomistic::requestAtoms(a);
  for(unsigned i=0; i<3; ++i) {
    getPntrToComponent(i)->resizeDerivatives(3*a.size()+9);
  }
}

void ActionWithVirtualAtom::apply() {
  if( !checkForForces() ) {
    return ;
  }

  Value* xval=getPntrToComponent(0);
  Value* yval=getPntrToComponent(1);
  Value* zval=getPntrToComponent(2);
  if( !xval->forcesWereAdded() && !yval->forcesWereAdded() && !zval->forcesWereAdded() ) {
    return ;
  }
  if( xval->isConstant() && yval->isConstant() && zval->isConstant() ) {
    return;
  }

  for(unsigned i=0; i<value_depends.size(); ++i) {
    xpos[value_depends[i]]->hasForce = true;
    ypos[value_depends[i]]->hasForce = true;
    zpos[value_depends[i]]->hasForce = true;
  }
  unsigned k=0;
  double xf = xval->inputForce[0];
  double yf = yval->inputForce[0];
  double zf = zval->inputForce[0];
// for(const auto & a : atom_value_ind) {
//   std::size_t nn = a.first, kk = a.second;
//   xpos[nn]->inputForce[kk] += xf*xval->data[1+k] + yf*yval->data[1+k] + zf*zval->data[1+k]; k++;
//   ypos[nn]->inputForce[kk] += xf*xval->data[1+k] + yf*yval->data[1+k] + zf*zval->data[1+k]; k++;
//   zpos[nn]->inputForce[kk] += xf*xval->data[1+k] + yf*yval->data[1+k] + zf*zval->data[1+k]; k++;
// }
  for(const auto & a : atom_value_ind_grouped) {
    const auto nn=a.first;
    auto & xp=xpos[nn]->inputForce;
    auto & yp=ypos[nn]->inputForce;
    auto & zp=zpos[nn]->inputForce;
    for(const auto & kk : a.second) {
      xp[kk] += xf*xval->data[1+k] + yf*yval->data[1+k] + zf*zval->data[1+k];
      k++;
      yp[kk] += xf*xval->data[1+k] + yf*yval->data[1+k] + zf*zval->data[1+k];
      k++;
      zp[kk] += xf*xval->data[1+k] + yf*yval->data[1+k] + zf*zval->data[1+k];
      k++;
    }
  }

  std::array<double,9> virial;
  for(unsigned i=0; i<9; ++i) {
    virial[i] = xf*xval->data[1+k] + yf*yval->data[1+k] + zf*zval->data[1+k];
    k++;
  }
  unsigned ind = 0;
  setForcesOnCell( virial.data(), virial.size(), ind );
  // The above can be achieved using the two functions below.  The code above that is specialised for the ActionWithVirtualAtom
  // class runs faster than the general code below.
  // if( !checkForForces() ) return ;
  // unsigned mm=0; setForcesOnAtoms( getForcesToApply(), mm );
}

void ActionWithVirtualAtom::setBoxDerivatives(const std::vector<Tensor> &d) {
  plumed_assert(d.size()==3);
  unsigned nbase = 3*getNumberOfAtoms();
  for(unsigned i=0; i<3; ++i) {
    Value* myval = getPntrToComponent(i);
    for(unsigned j=0; j<3; ++j) {
      for(unsigned k=0; k<3; ++k) {
        myval->setDerivative( nbase + 3*j + k, d[i][j][k] );
      }
    }
  }
// Subtract the trivial part coming from a distorsion applied to the ghost atom first.
// Notice that this part alone should exactly cancel the already accumulated virial
// due to forces on this atom.
  Vector pos;
  for(unsigned i=0; i<3; ++i) {
    pos[i] = getPntrToComponent(i)->get();
  }
  for(unsigned i=0; i<3; i++)
    for(unsigned j=0; j<3; j++) {
      getPntrToComponent(j)->addDerivative( nbase + 3*i + j, pos[i] );
    }
}

void ActionWithVirtualAtom::setBoxDerivativesNoPbc() {
  std::vector<Tensor> bd(3);
  for(unsigned i=0; i<3; i++)
    for(unsigned j=0; j<3; j++)
      for(unsigned k=0; k<3; k++) {
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
