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

namespace PLMD {

PLUMED_REGISTER_ACTION(ActionToPutData,"PUT")

void ActionToPutData::registerKeywords(Keywords& keys){
  ActionForInterface::registerKeywords( keys );
  keys.add("compulsory","SHAPE","0","the shape of the value that is being passed to PLUMED");
  keys.add("compulsory","UNIT","the unit of the quantity that is being passed to PLUMED through this value.  Can be either number, energy, length, mass or charge");
  keys.add("compulsory","FORCE_UNIT","default","the units to use for the force");
  keys.add("compulsory","PERIODIC","if the value being passed to plumed is periodic then you should specify the periodicity of the function.  If the value "
                                   "is not periodic you must state this using PERIODIC=NO.  Positions are passed with PERIODIC=NO even though special methods are used "
                                   "to deal with pbc");
  keys.addFlag("CONSTANT",false,"does this quantity not depend on time");
}

ActionToPutData::ActionToPutData(const ActionOptions&ao):
Action(ao),
ActionForInterface(ao),
noforce(false),
fixed(getName()=="TIMESTEP"),
dataCanBeSet(true),
mydata(DataPassingObject::create(plumed.getRealPrecision()))
{
   if( getName()!="ENERGY" && getName()!="PBC" && getName()!="TIMESTEP" ) {
       std::vector<unsigned> shape; parseVector("SHAPE",shape);
       if( shape.size()==1 && shape[0]==0 ) { shape.resize(0); addValue( shape ); }
       else { addValue( shape ); }    

       std::string unitstr, funitstr; parse("UNIT",unitstr); 
       parse("FORCE_UNIT",funitstr); setUnit( unitstr, funitstr );

       // Now sort out period
       std::vector<std::string> period; parseVector("PERIODIC",period);
       if( period.size()==1 ) {
           if( period[0]!="NO") error("input to PERIODIC keyword does not make sense");
           setNotPeriodic();
       } else if( period.size()==2 ) setPeriodic( period[0], period[1] );    
       else  error("input to PERIODIC keyword does not make sense");

       parseFlag("CONSTANT",fixed); if( fixed ) { noforce=true; getPntrToOutput(0)->setConstant(); } 
   }
}

void ActionToPutData::setUnit( const std::string& unitstr, const std::string& funitstr ) {
   if( unitstr=="number" ) unit=n;
   else if( unitstr=="energy" ) unit=e;
   else if( unitstr=="length" ) unit=l;
   else if( unitstr=="mass" ) unit=m;
   else if( unitstr=="charge" ) unit=q;
   else if( unitstr=="time" ) unit=t;
   else error( unitstr + " is not a valid input unit");
   // Set the force units
   if( funitstr=="default" ) funit=d;
   else if( funitstr=="energy" ) funit=eng;
   else error( funitstr + " is not a valid input force unit");
}

void ActionToPutData::setStride( const std::string& name, const unsigned& sss ) {
  plumed_assert( name==getLabel() ); mydata->setStride(sss);
}

void ActionToPutData::updateUnits( const Units& MDUnits, const Units& units ) {
  // Don't need to do anythign if this is just a number
  if( unit==n ) return ; 

  double vunits; 
  if( unit==e ) vunits = MDUnits.getEnergy()/units.getEnergy();  
  else if( unit==l ) vunits = MDUnits.getLength()/units.getLength(); 
  else if( unit==m ) vunits = MDUnits.getMass()/units.getMass();
  else if( unit==q ) vunits = MDUnits.getCharge()/units.getCharge();
  else if( unit==t ) vunits = MDUnits.getTime()/units.getTime();
  mydata->setUnit(vunits); if( fixed ) transferFixedValue();
  if( funit==eng ) mydata->setForceUnit(units.getEnergy()/MDUnits.getEnergy());
  else if( funit==d ) mydata->setForceUnit((units.getEnergy()/MDUnits.getEnergy())*vunits);
}

double ActionToPutData::MD2double( void* val ) const {
  return mydata->MD2double(val);
}

bool ActionToPutData::setValuePointer( const std::string& name, void* val ) {
   if( name!=getLabel() ) return false;
   wasset=true; plumed_massert( dataCanBeSet, "set " + getLabel() + " cannot be set at this time");
   mydata->setValuePointer(val); return true;
}

bool ActionToPutData::setForcePointer( const std::string& name, void* val ) {
   if( name!=getLabel() ) return false;
   plumed_massert( dataCanBeSet, "force on " + getLabel() + " cannot be set at this time");
   mydata->setForcePointer(val); return true;
}

void ActionToPutData::transferFixedValue() {
   plumed_assert( fixed ); if( !wasset ) return;
   mydata->share_data( 0, getPntrToValue()->getNumberOfValues(), getPntrToValue() );
}

void ActionToPutData::wait() {
   dataCanBeSet=false; if( fixed || !wasset ) { return; } plumed_assert( wasset ); 
   mydata->share_data( 0, getPntrToValue()->getNumberOfValues(), getPntrToValue() );
}

void ActionToPutData::apply() {
   if( getPntrToValue()->forcesWereAdded() && !noforce ) {
       if( getName()=="ENERGY" || getDependencies().size()==0 ) mydata->add_force( getPntrToValue() );
   }
}

unsigned ActionToPutData::getNumberOfForcesToRescale() const {
   if( getName()!="ENERGY" || getDependencies().size()>0 ) return getPntrToOutput(0)->getNumberOfValues();
   plumed_assert( getDependencies().size()==1 ); ActionForInterface* ai = dynamic_cast<ActionForInterface*>( getDependencies()[0] );
   return ai->getNumberOfForcesToRescale();
}

void ActionToPutData::rescaleForces( const double& alpha ) {
   if( noforce ) return; wasscaled=true;
   mydata->rescale_force( getNumberOfForcesToRescale(), alpha, getPntrToValue() );
    
}

void ActionToPutData::writeBinary(std::ostream&o) {
  if(!fixed) getPntrToValue()->writeBinary(o);
}

void ActionToPutData::readBinary(std::istream&i) {
  if(!fixed) getPntrToValue()->readBinary(i);
}

}
