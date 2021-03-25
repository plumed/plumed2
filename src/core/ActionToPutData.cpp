/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2017-2020 The plumed team
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
#include "ActionToPutData.h"
#include "ActionRegister.h"
#include "PlumedMain.h"
#include "ActionSet.h"
#include "Atoms.h"

namespace PLMD {

PLUMED_REGISTER_ACTION(ActionToPutData,"PUT")

void ActionToPutData::registerKeywords(Keywords& keys){
  Action::registerKeywords(keys); ActionWithValue::registerKeywords( keys );
  keys.add("compulsory","SHAPE","0","the shape of the value that is being passed to PLUMED");
  keys.add("optional","FORCES_FOR_POTENTIAL","If your input quantity is an energy this lists the input actions that hold the forces.  These are rescaled");
  keys.addFlag("SUM_OVER_DOMAINS",false,"does this quantity need to be summed over domains");
  keys.addFlag("NOFORCE",false,"always set the forces on this value to zero");
}

ActionToPutData::ActionToPutData(const ActionOptions&ao):
Action(ao),
ActionWithValue(ao),
wasset(false),
mydata(DataPassingObject::create(plumed.getAtoms().getRealPrecision()))
{
   std::vector<unsigned> shape; parseVector("SHAPE",shape);
   if( shape.size()==1 && shape[0]==0 ) { shape.resize(0); addValue( shape ); }
   else { addValue( shape ); }    
   parseFlag("SUM_OVER_DOMAINS",sum_domains); parseFlag("NOFORCE", noforce);
   std::vector<std::string> toscale; parseVector("FORCES_FOR_POTENTIAL",toscale);
   for(unsigned i=0;i<toscale.size();++i) {
       plumed_assert( shape.size()==0 );
       ActionToPutData* ap=plumed.getActionSet().selectWithLabel< ActionToPutData*>(toscale[i]);
       plumed_assert(ap); forces_to_scale.push_back(ap);  
   }
}

void ActionToPutData::setUnit( const double& u ) {
   mydata->setUnit(u);
}

void ActionToPutData::setForceUnit( const double& u ) {
   mydata->setForceUnit(u);
}

void ActionToPutData::set_value(void* val ) {
   wasset=true; mydata->setValuePointer(val);
}

void ActionToPutData::set_force(void* val ) {
   mydata->setForcePointer(val);
}

void ActionToPutData::set_domain( const bool& periodic, const std::string& min, const std::string& max ) {
   if( periodic ) setPeriodic( min, max ); else setNotPeriodic();
}

void ActionToPutData::interpretDataLabel( const std::string& mystr, Action* myuser, unsigned& nargs, std::vector<Value*>& args ) {
   if( mystr.find("*")!=std::string::npos ) return ;
   ActionWithValue::interpretDataLabel( mystr, myuser, nargs, args );
}

void ActionToPutData::wait() {
   plumed_assert( wasset ); mydata->share_data( getPntrToValue() );
}

void ActionToPutData::apply() {
   if( getPntrToValue()->forcesWereAdded() && !noforce ) {
       for(unsigned i=0;i<forces_to_scale.size();++i) forces_to_scale[i]->rescaleForces( 1.- getPntrToValue()->getForce(0)); 
       mydata->add_force( getPntrToValue() );
   }
}

void ActionToPutData::rescaleForces( const double& alpha ) {
   if( noforce ) return;
   mydata->rescale_force( alpha, getPntrToValue() );
}

void ActionToPutData::writeBinary(std::ostream&o) {
  // if( fixed ) return; 
  getPntrToValue()->writeBinary(o);
}

void ActionToPutData::readBinary(std::istream&i) {
  // if( fixed ) return; 
  getPntrToValue()->readBinary(i);
}

}
