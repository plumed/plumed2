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
#include "Bias.h"


namespace PLMD {
namespace bias {

Bias::Bias(const ActionOptions&ao):
  Action(ao),
  ActionPilot(ao),
  ActionWithValue(ao),
  ActionWithArguments(ao),
  outputForces(getNumberOfArguments(),0.0) {
  addComponentWithDerivatives("bias");
  componentIsNotPeriodic("bias");
  valueBias=getPntrToComponent("bias");

  if(getStride()>1) {
    log<<"  multiple time step "<<getStride()<<" ";
    log<<cite("Ferrarotti, Bottaro, Perez-Villa, and Bussi, J. Chem. Theory Comput. 11, 139 (2015)")<<"\n";
  }
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    (getPntrToArgument(i)->getPntrToAction())->turnOnDerivatives();
  }

  ActionWithValue::turnOnDerivatives();
}

void Bias::registerKeywords( Keywords& keys ) {
  Action::registerKeywords(keys);
  ActionPilot::registerKeywords(keys);
  ActionWithValue::registerKeywords(keys);
  ActionWithArguments::registerKeywords(keys);
  keys.addInputKeyword("compulsory","ARG","scalar","the labels of the scalars on which the bias will act");
  keys.add("hidden","STRIDE","the frequency with which the forces due to the bias should be calculated.  This can be used to correctly set up multistep algorithms");
  keys.addOutputComponent("bias","default","scalar","the instantaneous value of the bias potential");
}

void Bias::apply() {
  const unsigned noa=getNumberOfArguments();
  const unsigned ncp=getNumberOfComponents();

  if(onStep()) {
    double gstr = static_cast<double>(getStride());
    for(unsigned i=0; i<noa; ++i) {
      getPntrToArgument(i)->addForce(gstr*outputForces[i]);
    }
  }

  accumForces.assign(noa,0.0);
  forces.resize(noa);

  bool at_least_one_forced=false;
  for(unsigned i=0; i<ncp; ++i) {
    if(getPntrToComponent(i)->applyForce(forces)) {
      at_least_one_forced=true;
      for(unsigned j=0; j<noa; j++) {
        accumForces[j]+=forces[j];
      }
      //it may be faster like this(needs some measurements):
      //(it produces ~20 less asm instructions)
      /*
      auto aF = accumForces.begin();
      auto f = forces.begin();
      //accumForces and forces have the same size:
      while (aF!=accumForces.end()) {
        *aF+=*f;
        ++aF;
        ++f;
      }
      */
    }
  }

  if(at_least_one_forced && !onStep()) {
    error("you are biasing a bias with an inconsistent STRIDE");
  }

  if(at_least_one_forced)
    for(unsigned i=0; i<noa; ++i) {
      getPntrToArgument(i)->addForce(accumForces[i]);
    }

}

}
}


