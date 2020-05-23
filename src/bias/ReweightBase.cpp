/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2020 The plumed team
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
#include "ReweightBase.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"

namespace PLMD {
namespace bias {

void ReweightBase::registerKeywords(Keywords& keys) {
  Action::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys );
  keys.setComponentsIntroduction("This action calculates the logarithm of a weight for reweighting");
  keys.add("optional","TEMP","the system temperature.  This is not required if your MD code passes this quantity to PLUMED");
  keys.remove("NUMERICAL_DERIVATIVES");
}

ReweightBase::ReweightBase(const ActionOptions&ao):
  Action(ao),
  ActionWithValue(ao),
  ActionWithArguments(ao)
{
  simtemp=0.; parse("TEMP",simtemp);
  if(simtemp>0) simtemp*=plumed.getAtoms().getKBoltzmann();
  else simtemp=plumed.getAtoms().getKbT();
  if(simtemp==0) error("The MD engine does not pass the temperature to plumed so you have to specify it using TEMP");
  // Create something to hold the weight
  addValue(); setNotPeriodic();
}

void ReweightBase::calculate() {
  double weight = getLogWeight();
  setValue( weight );
}

}
}
