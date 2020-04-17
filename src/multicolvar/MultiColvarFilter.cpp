/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2020 The plumed team
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
#include "MultiColvarFilter.h"

namespace PLMD {
namespace multicolvar {

void MultiColvarFilter::registerKeywords( Keywords& keys ) {
  BridgedMultiColvarFunction::registerKeywords( keys );
  if( keys.reserved("VMEAN") ) keys.use("VMEAN");
  keys.use("MEAN"); keys.use("MOMENTS"); keys.use("MIN"); keys.use("MAX");
  keys.use("ALT_MIN"); keys.use("LOWEST"); keys.use("HIGHEST");
}

MultiColvarFilter::MultiColvarFilter(const ActionOptions&ao):
  Action(ao),
  BridgedMultiColvarFunction(ao)
{
  if( getPntrToMultiColvar()->isDensity() ) error("filtering/transforming density makes no sense");

  if( getName().find("MFILTER")!=std::string::npos ) filter=true;
  else {
    plumed_assert( getName().find("MTRANSFORM")!=std::string::npos );
    filter=false;
  }

  readVesselKeywords();
}

void MultiColvarFilter::doJobsRequiredBeforeTaskList() {
  ActionWithValue::clearDerivatives();
  ActionWithVessel::doJobsRequiredBeforeTaskList();
}

void MultiColvarFilter::completeTask( const unsigned& curr, MultiValue& invals, MultiValue& outvals ) const {
  invals.copyValues( outvals );
  if( derivativesAreRequired() ) invals.copyDerivatives( outvals );

  // Retrieve the value of the multicolvar and apply filter
  double val=invals.get(1), df, weight=applyFilter( val, df );

  // Now propegate derivatives
  if( filter && !getPntrToMultiColvar()->weightHasDerivatives ) {
    outvals.setValue( 0, weight );
    if( derivativesAreRequired() ) {
      for(unsigned i=0; i<invals.getNumberActive(); ++i) {
        unsigned jder=invals.getActiveIndex(i);
        outvals.addDerivative( 0, jder, df*invals.getDerivative(1, jder ) );
      }
    }
  } else if( filter ) {
    double ww=outvals.get(0); outvals.setValue( 0, ww*weight );
    if( derivativesAreRequired() ) {
      for(unsigned i=0; i<outvals.getNumberActive(); ++i) {
        unsigned ider=outvals.getActiveIndex(i);
        outvals.setDerivative( 0, ider, weight*outvals.getDerivative(1,ider) + ww*df*outvals.getDerivative(0,ider) );
      }
    }
  } else {
    outvals.setValue( 1, weight );
    if( derivativesAreRequired() ) {
      for(unsigned i=0; i<invals.getNumberActive(); ++i) {
        unsigned jder=invals.getActiveIndex(i);
        outvals.setDerivative( 1, jder, df*invals.getDerivative(1, jder ) );
      }
    }
  }
}

void MultiColvarFilter::addBridgeForces( const std::vector<double>& bb ) {
  plumed_dbg_assert( bb.size()==0 );
}

}
}
