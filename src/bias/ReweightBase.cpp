/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2016 The plumed team
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
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/Atoms.h"
#include "ReweightBase.h"

namespace PLMD {
namespace bias {

void ReweightBase::registerKeywords(Keywords& keys){
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

void ReweightBase::retrieveAllBiases( const std::string& lab, std::vector<Value*>& vals ){
   std::vector<ActionWithValue*> all=plumed.getActionSet().select<ActionWithValue*>();
   if( all.empty() ) error("your input file is not telling plumed to calculate anything");
   log.printf("  using the following biases in reweighting ");
   for(unsigned j=0;j<all.size();j++){
       std::string flab; flab=all[j]->getLabel() + "." + lab;
       if( all[j]->exists(flab) ){
           vals.push_back( all[j]->copyOutput(flab) );
           log.printf(" %s", flab.c_str());
       }
   }
   log.printf("\n");
   if( !vals.empty() ) requestArguments( vals );   
}

void ReweightBase::calculate(){
  double weight = getLogWeight();
  setValue( weight );
}

}
}
