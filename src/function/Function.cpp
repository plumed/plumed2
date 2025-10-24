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
#include "Function.h"
#include "tools/OpenMP.h"
#include "tools/Communicator.h"

namespace PLMD {
namespace function {

void Function::registerKeywords(Keywords& keys) {
  Action::registerKeywords(keys);
  ActionWithValue::registerKeywords(keys);
  ActionWithArguments::registerKeywords(keys);
  keys.addInputKeyword("compulsory","ARG","scalar","the labels of the values from which the function is calculated");
  keys.reserve("compulsory","PERIODIC","if the output of your function is periodic then you should specify the periodicity of the function.  If the output is not periodic you must state this using PERIODIC=NO");
}

Function::Function(const ActionOptions&ao):
  Action(ao),
  ActionWithValue(ao),
  ActionWithArguments(ao) {
}

void Function::addValueWithDerivatives(const std::vector<std::size_t>& /* shape */) {
//the shape argument is ignoded (shall I put here a warning?)
  plumed_massert( getNumberOfArguments()!=0, "for functions you must requestArguments before adding values");
  ActionWithValue::addValueWithDerivatives();
  getPntrToValue()->resizeDerivatives(getNumberOfDerivatives());
}

void Function::addComponentWithDerivatives( const std::string& compName, const std::vector<std::size_t>& /* shape */) {
//the shape argument is ignoded (shall I put here a warning?)
  plumed_massert( getNumberOfArguments()!=0, "for functions you must requestArguments before adding values");
  ActionWithValue::addComponentWithDerivatives(compName);
  getPntrToComponent(compName)->resizeDerivatives(getNumberOfDerivatives());
}

void Function::apply() {
  if( !checkForForces() ) {
    return;
  }
  unsigned ind=0;
  addForcesOnArguments( 0, getForcesToApply(), ind );
}

}
}
