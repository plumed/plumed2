/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2020 The plumed team
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
#include "FunctionVessel.h"
#include "core/ActionWithValue.h"

namespace PLMD {
namespace vesselbase {

void FunctionVessel::registerKeywords( Keywords& keys ) {
  ValueVessel::registerKeywords( keys );
}

FunctionVessel::FunctionVessel( const VesselOptions& da ):
  ValueVessel(da),
  norm(false),
  usetol(false)
{
  diffweight=getAction()->weightHasDerivatives;
}

void FunctionVessel::resize() {
  if( getAction()->derivativesAreRequired() ) {
    unsigned nderivatives=getAction()->getNumberOfDerivatives();
    getFinalValue()->resizeDerivatives( nderivatives );
    resizeBuffer( (1+nderivatives)*2 );
    diffweight=getAction()->weightHasDerivatives;
  } else {
    resizeBuffer(2);
    diffweight=false;  // Don't need to worry about differentiable weights if no derivatives
  }
}

void FunctionVessel::calculate( const unsigned& current, MultiValue& myvals, std::vector<double>& buffer, std::vector<unsigned>& der_list ) const {
  unsigned nderivatives=getFinalValue()->getNumberOfDerivatives();
  double weight=myvals.get(0);
  plumed_dbg_assert( weight>=getTolerance() );

  // This deals with the value
  double dval, f=calcTransform( myvals.get(mycomp), dval );

  if( norm ) {
    if( usetol && weight<getTolerance() ) return;
    buffer[bufstart+1+nderivatives] += weight;
    if( getAction()->derivativesAreRequired() && diffweight ) myvals.chainRule( 0, 1, 1, 0, 1.0, bufstart, buffer );
  }

  double contr=weight*f;
  if( usetol && contr<getTolerance() ) return;
  buffer[bufstart] += contr;

  if( diffweight ) myvals.chainRule( 0, 0, 1, 0, f, bufstart, buffer );
  if( getAction()->derivativesAreRequired() && fabs(dval)>0.0 ) myvals.chainRule( mycomp, 0, 1, 0, weight*dval, bufstart, buffer );

  return;
}

double FunctionVessel::calcTransform( const double&, double& ) const {
  plumed_error(); return 1.0;
}

void FunctionVessel::finish( const std::vector<double>& buffer ) {
  unsigned nderivatives=getFinalValue()->getNumberOfDerivatives();
  if( norm && diffweight ) {
    double dv, val=finalTransform( buffer[bufstart], dv), weight=buffer[bufstart+1+nderivatives];
    getFinalValue()->set( val / weight );
    for(unsigned i=0; i<nderivatives; ++i) {
      getFinalValue()->addDerivative( i, buffer[bufstart+1+i]/weight - val*buffer[bufstart+1+nderivatives+1+i]/(weight*weight) );
    }
  } else if( norm ) {
    double dv, val=finalTransform( buffer[bufstart], dv), weight=buffer[bufstart+1+nderivatives];
    getFinalValue()->set( val / weight );
    for(unsigned i=0; i<nderivatives; ++i) getFinalValue()->addDerivative( i, buffer[bufstart+1+i]/weight );
  } else {
    double dv, val=finalTransform( buffer[bufstart], dv); getFinalValue()->set( val );
    for(unsigned i=0; i<nderivatives; ++i) getFinalValue()->addDerivative( i, dv*buffer[bufstart+1+i] );
  }
}

double FunctionVessel::finalTransform( const double& val, double& dv ) {
  dv=1.0; return val;
}

}
}

