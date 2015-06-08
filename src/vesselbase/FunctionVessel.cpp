/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013,2014 The plumed team
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
#include "FunctionVessel.h"
#include "core/ActionWithValue.h"

namespace PLMD{
namespace vesselbase{

void FunctionVessel::registerKeywords( Keywords& keys ){
  Vessel::registerKeywords( keys );
}

FunctionVessel::FunctionVessel( const VesselOptions& da ):
Vessel(da),
norm(false),
usetol(false),
nderivatives(0)
{
  ActionWithValue* a=dynamic_cast<ActionWithValue*>( getAction() );
  plumed_massert(a,"cannot create passable values as base action does not inherit from ActionWithValue");
  int numval = getNumericalLabel();
  if( numval<0 ){   // This allows us to make multicolvars pretend to be colvars - this is used in AlphaRMSD etc
     plumed_massert( a->getNumberOfComponents()==0,"you can't multiple values with the action label");
     a->addValueWithDerivatives(); 
     a->setNotPeriodic();
  } else {
     plumed_massert( !a->exists(getAction()->getLabel() + "." + getLabel() ), "you can't create the name multiple times");
     a->addComponentWithDerivatives( getLabel() ); 
     a->componentIsNotPeriodic( getLabel() );
  }
  final_value=a->copyOutput( a->getNumberOfComponents()-1 );
  diffweight=getAction()->weightHasDerivatives;
}

std::string FunctionVessel::description(){
  if( final_value->getName()==getAction()->getLabel() ) return "value " + getAction()->getLabel() + " contains " + function_description();
  return "value " + getAction()->getLabel() + "." + getLabel() + " contains " + function_description();
}

void FunctionVessel::resize(){
  if( getAction()->derivativesAreRequired() ){
     nderivatives=getAction()->getNumberOfDerivatives();
     resizeBuffer( (1+nderivatives)*2 ); 
     final_value->resizeDerivatives( nderivatives );
     diffweight=getAction()->weightHasDerivatives;
  } else {
     nderivatives=0;
     resizeBuffer(2);
     diffweight=false;  // Don't need to worry about differentiable weights if no derivatives
  }
}

void FunctionVessel::setNumberOfDerivatives( const unsigned& nder ){
  nderivatives=nder;
  final_value->resizeDerivatives( nder );
}

bool FunctionVessel::calculate( const unsigned& current, MultiValue& myvals, std::vector<double>& buffer, std::vector<unsigned>& der_list ) const {
  double weight=myvals.get(0); 
  plumed_dbg_assert( weight>=getTolerance() );  

  // This deals with the value
  double dval, f=calcTransform( myvals.get(1), dval );

  if( norm ){
     if( usetol && weight<getTolerance() ) return false;
     buffer[bufstart+1+nderivatives] += weight;
     if( diffweight ) myvals.chainRule( 0, 1, 1, 0, 1.0, bufstart, buffer );
  }

  double contr=weight*f;
  if( usetol && contr<getTolerance() ) return false;
  buffer[bufstart] += contr;

  if( diffweight ) myvals.chainRule( 0, 0, 1, 0, f, bufstart, buffer ); 
  if( getAction()->derivativesAreRequired() && fabs(dval)>0.0 ) myvals.chainRule( 1, 0, 1, 0, weight*dval, bufstart, buffer );

  return true;
}

double FunctionVessel::calcTransform( const double& , double& ) const { 
  plumed_error(); return 1.0; 
}

void FunctionVessel::finish( const std::vector<double>& buffer ){
  if( norm && diffweight ){
      double dv, val=finalTransform( buffer[bufstart], dv), weight=buffer[bufstart+1+nderivatives];
      final_value->set( val / weight );
      for(unsigned i=0;i<nderivatives;++i){
         final_value->addDerivative( i, buffer[bufstart+1+i]/weight - val*buffer[bufstart+1+nderivatives+1+i]/(weight*weight) );
      }
  } else if( norm ){
     double dv, val=finalTransform( buffer[bufstart], dv), weight=buffer[bufstart+1+nderivatives];
     final_value->set( val / weight );
     for(unsigned i=0;i<nderivatives;++i) final_value->addDerivative( i, buffer[bufstart+1+i]/weight );
  } else {
     double dv, val=finalTransform( buffer[bufstart], dv); final_value->set( val );
     for(unsigned i=0;i<nderivatives;++i) final_value->addDerivative( i, dv*buffer[bufstart+1+i] );
  }
}

double FunctionVessel::finalTransform( const double& val, double& dv ){
  dv=1.0; return val;
}

bool FunctionVessel::applyForce( std::vector<double>& forces ){
  std::vector<double> tmpforce( forces.size() );
  forces.assign(forces.size(),0.0); bool wasforced=false;
  if( final_value->applyForce( tmpforce ) ){
      wasforced=true;
      for(unsigned j=0;j<forces.size();++j) forces[j]+=tmpforce[j];
  }
  return wasforced;
}

}
}

