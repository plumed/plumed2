/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2019 The plumed team
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
#include "OrderingVessel.h"
#include "core/ActionWithValue.h"

namespace PLMD {
namespace vesselbase {


void OrderingVessel::registerKeywords( Keywords& keys ) {
  ValueVessel::registerKeywords( keys );
}

OrderingVessel::OrderingVessel( const VesselOptions& da ) :
  ValueVessel(da)
{
  mydata = getAction()->buildDataStashes( NULL );
  for(unsigned i=0; i<getAction()->getNumberOfVessels(); ++i) {
    if( getAction()->getPntrToVessel(i)->getName()==getName() )
      error("calculating lowest/highest value multiple times serves no purpose");
  }
}

void OrderingVessel::resize() {
  resizeBuffer( 0 );
  if( getAction()->derivativesAreRequired() ) getFinalValue()->resizeDerivatives( getAction()->getNumberOfDerivatives() );
}

void OrderingVessel::finish( const std::vector<double>& buffer ) {
  std::vector<double> values( getAction()->getNumberOfQuantities() );
  mydata->retrieveSequentialValue( 0, false, values );

  double min=values[mycomp]; unsigned mini=getAction()->getPositionInFullTaskList(0);
  for(unsigned i=1; i<mydata->getNumberOfStoredValues(); ++i) {
    mydata->retrieveSequentialValue( i, false, values );
    double newval = values[mycomp];
    if( compare( newval, min ) ) { min=newval; mini=getAction()->getPositionInFullTaskList(i); }
  }
  setOutputValue( min );

  if( getAction()->derivativesAreRequired() ) {
    MultiValue myvals( getAction()->getNumberOfQuantities(), getAction()->getNumberOfDerivatives() );
    mydata->retrieveDerivatives( mini, false, myvals ); Value* fval=getFinalValue();
    for(unsigned i=0; i<myvals.getNumberActive(); ++i) {
      unsigned ider=myvals.getActiveIndex(i);
      fval->setDerivative( ider, myvals.getDerivative(mycomp,ider) );
    }
  }
}

}
}
