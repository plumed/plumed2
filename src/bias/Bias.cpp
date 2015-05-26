/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2015 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

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


namespace PLMD{
namespace bias{

Bias::Bias(const ActionOptions&ao):
Action(ao),
ActionPilot(ao),
ActionWithValue(ao),
ActionWithArguments(ao),
outputForces(getNumberOfArguments(),0.0)
{
  if(getStride()>1){
    log<<"  multiple time step "<<getStride()<<" ";
    log<<cite("Ferrarotti, Bottaro, Perez-Villa, and Bussi, J. Chem. Theory Comput. 11, 139 (2015)")<<"\n";
  }
  for(unsigned i=0;i<getNumberOfArguments();++i){
     (getPntrToArgument(i)->getPntrToAction())->turnOnDerivatives();
  }
}

void Bias::registerKeywords( Keywords& keys ){
  Action::registerKeywords(keys);
  ActionPilot::registerKeywords(keys);
  ActionWithValue::registerKeywords(keys);
  ActionWithArguments::registerKeywords(keys);
  keys.add("hidden","STRIDE","the frequency with which the forces due to the bias should be calculated.  This can be used to correctly set up multistep algorithms");
}

void Bias::turnOnDerivatives(){
  error("a bias cannot be used as a collective variable");
}

void Bias::apply(){
  if(onStep()) for(unsigned i=0;i<getNumberOfArguments();++i){
    getPntrToArgument(i)->addForce(getStride()*outputForces[i]);
  }
}

}
}


