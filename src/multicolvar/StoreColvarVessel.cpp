/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2014 The plumed team
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
#include "StoreColvarVessel.h"
#include "MultiColvarFunction.h"

namespace PLMD {
namespace multicolvar {

void StoreColvarVessel::registerKeywords( Keywords& keys ){
  StoreDataVessel::registerKeywords( keys );
}

StoreColvarVessel::StoreColvarVessel( const vesselbase::VesselOptions& da):
StoreDataVessel(da)
{
  completeSetup( 0, 1 );
}

void StoreColvarVessel::chainRuleForComponent( const unsigned& icolv, const unsigned& jout, const unsigned& base_cv_no, 
                                               const double& weight, MultiColvarFunction* funcout ){
  plumed_dbg_assert( getAction()->derivativesAreRequired() );  

  if( usingLowMem() ){
     unsigned ibuf = icolv*getAction()->getNumberOfDerivatives();
     for(unsigned ider=0;ider<getNumberOfDerivatives(icolv);++ider){
         funcout->addStoredDerivative( jout, base_cv_no, getStoredIndex( icolv, ider ), weight*getLocalDerivative(ibuf+ider) );
     } 
  } else {
     unsigned ibuf = icolv*getNumberOfDerivativeSpacesPerComponent() + 1;
     for(unsigned ider=0;ider<getNumberOfDerivatives(icolv);++ider){
         funcout->addStoredDerivative( jout, base_cv_no, getStoredIndex( icolv, ider ), weight*getBufferElement(ibuf+ider) );
     } 
  }
}

}
}
