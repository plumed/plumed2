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
  keys.add("compulsory","UNIT","the unit of the quantity that is being passed to PLUMED through this value.  Can be either number, energy, length, mass or charge");
  keys.add("compulsory","FORCE_UNIT","default","the units to use for the force");
  keys.add("optional","FORCES_FOR_POTENTIAL","If your input quantity is an energy this lists the input actions that hold the forces.  These are rescaled");
  keys.add("compulsory","PERIODIC","if the value being passed to plumed is periodic then you should specify the periodicity of the function.  If the value "
                                   "is not periodic you must state this using PERIODIC=NO.  Positions are passed with PERIODIC=NO even though special methods are used "
                                   "to deal with pbc");
  keys.addFlag("SUM_OVER_DOMAINS",false,"does this quantity need to be summed over domains");
  keys.addFlag("NOFORCE",false,"always set the forces on this value to zero");
  keys.addFlag("SCATTERED",false,"is this vector scattered over the domains");
  keys.addFlag("CONSTANT",false,"does this quantity not depend on time");
}

ActionToPutData::ActionToPutData(const ActionOptions&ao):
Action(ao),
ActionWithValue(ao),
wasset(false),
firststep(true),
dataCanBeSet(true),
wasscaled(false),
mydata(DataPassingObject::create(plumed.getRealPrecision()))
{
   std::vector<unsigned> shape; parseVector("SHAPE",shape);
   if( shape.size()==1 && shape[0]==0 ) { shape.resize(0); addValue( shape ); }
   else { addValue( shape ); }    

   std::string unitstr; parse("UNIT",unitstr);
   if( unitstr=="number" ) unit=n;
   else if( unitstr=="energy" ) unit=e;
   else if( unitstr=="length" ) unit=l;
   else if( unitstr=="mass" ) unit=m;
   else if( unitstr=="charge" ) unit=q;
   else error( unitstr + " is not a valid input unit");
   std::string funitstr; parse("FORCE_UNIT",funitstr);
   if( funitstr=="default" ) funit=d;
   else if( funitstr=="energy" ) funit=eng;
   else error( funitstr + " is not a valid input force unit");

   // Now sort out period
   std::vector<std::string> period; parseVector("PERIODIC",period);
   if( period.size()==1 ) {
       if( period[0]!="NO") error("input to PERIODIC keyword does not make sense");
       setNotPeriodic();
   } else if( period.size()==2 ) setPeriodic( period[0], period[1] );    
   else  error("input to PERIODIC keyword does not make sense");

   parseFlag("SUM_OVER_DOMAINS",sum_domains); parseFlag("NOFORCE", noforce); 
   parseFlag("CONSTANT",fixed); if( fixed ) { noforce=true; getPntrToOutput(0)->setConstant(); } 
   parseFlag("SCATTERED",scattered);  if( scattered ) plumed_assert( shape.size()>0 );
   std::vector<std::string> toscale; parseVector("FORCES_FOR_POTENTIAL",toscale);
   for(unsigned i=0;i<toscale.size();++i) {
       plumed_assert( shape.size()==0 );
       ActionToPutData* ap=plumed.getActionSet().selectWithLabel< ActionToPutData*>(toscale[i]);
       plumed_assert(ap); forces_to_scale.push_back(ap); addDependency( ap );  
   }
}

void ActionToPutData::setStride( const unsigned& sss ) {
  mydata->setStride(sss);
}

void ActionToPutData::updateUnits() {
  // Don't need to do anythign if this is just a number
  if( unit==n ) return ; 

  double vunits; 
  const Units& MDUnits = plumed.getAtoms().getMDUnits();
  const Units& units = plumed.getAtoms().getUnits();
  if( unit==e ) vunits = MDUnits.getEnergy()/units.getEnergy();  
  else if( unit==l ) vunits = MDUnits.getLength()/units.getLength(); 
  else if( unit==m ) vunits = MDUnits.getMass()/units.getMass();
  else if( unit==q ) vunits = MDUnits.getCharge()/units.getCharge();
  mydata->setUnit(vunits); 
  if( funit==eng ) mydata->setForceUnit(units.getEnergy()/MDUnits.getEnergy());
  else if( funit==d ) mydata->setForceUnit((units.getEnergy()/MDUnits.getEnergy())*vunits);
}

void ActionToPutData::set_value(void* val ) {
   wasset=true; plumed_massert( dataCanBeSet, "set " + getLabel() + " cannot be set at this time");
   mydata->setValuePointer(val);
}

void ActionToPutData::set_force(void* val ) {
   plumed_massert( dataCanBeSet, "force on " + getLabel() + " cannot be set at this time");
   mydata->setForcePointer(val);
}

void ActionToPutData::share( const unsigned& j, const unsigned& k ) {
  plumed_assert( scattered ); if( wasset ) mydata->share_data( j, k, getPntrToValue() ); 
}

void ActionToPutData::share( const std::set<AtomNumber>&index, const std::vector<unsigned>& i ) {
   plumed_assert( scattered && getPntrToValue()->getRank()==1 ); if( wasset ) mydata->share_data( index, i,  getPntrToValue() );
}

void ActionToPutData::wait() {
   dataCanBeSet=false; if( fixed || scattered || !wasset ) { return; } plumed_assert( wasset ); 
   mydata->share_data( 0, getPntrToValue()->getNumberOfValues(), getPntrToValue() );
}

void ActionToPutData::apply() {
   if( getPntrToValue()->forcesWereAdded() && !noforce ) {
       if( forces_to_scale.size()>0 ) {
           for(unsigned i=0;i<forces_to_scale.size();++i) forces_to_scale[i]->rescaleForces( 1.- getPntrToValue()->getForce(0)); 
       } else if( !scattered ) {
           mydata->add_force( getPntrToValue() );
       } else if( wasscaled || (int(plumed.getAtoms().getGatindex().size())==plumed.getAtoms().getNatoms() && plumed.getAtoms().shuffledAtoms==0) ) {
           mydata->add_force( plumed.getAtoms().getGatindex(), getPntrToValue() );
       } else { mydata->add_force( plumed.getAtoms().unique, plumed.getAtoms().uniq_index, getPntrToValue() ); }
   }
}

void ActionToPutData::rescaleForces( const double& alpha ) {
   if( noforce ) return; wasscaled=true;
   if( scattered ) mydata->rescale_force( plumed.getAtoms().getGatindex().size(), alpha, getPntrToValue() );
   else mydata->rescale_force( getPntrToValue()->getNumberOfValues(), alpha, getPntrToValue() );
}

void ActionToPutData::writeBinary(std::ostream&o) {
  if(!fixed) getPntrToValue()->writeBinary(o);
}

void ActionToPutData::readBinary(std::istream&i) {
  if(!fixed) getPntrToValue()->readBinary(i);
}

}
