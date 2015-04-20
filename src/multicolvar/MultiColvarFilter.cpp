/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014 The plumed team
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
#include "MultiColvarFilter.h"

namespace PLMD {
namespace multicolvar {

void MultiColvarFilter::registerKeywords( Keywords& keys ){
  BridgedMultiColvarFunction::registerKeywords( keys );
  if( keys.reserved("VMEAN") ) keys.use("VMEAN");
  keys.use("MEAN"); keys.use("MOMENTS"); keys.use("MIN"); keys.use("MAX");
}

MultiColvarFilter::MultiColvarFilter(const ActionOptions&ao):
Action(ao),
BridgedMultiColvarFunction(ao)
{
  if( getPntrToMultiColvar()->isDensity() ) error("filtering density makes no sense");
  readVesselKeywords(); 
}

void MultiColvarFilter::doJobsRequiredBeforeTaskList(){
  ActionWithValue::clearDerivatives();
  ActionWithVessel::doJobsRequiredBeforeTaskList();
}

void MultiColvarFilter::completeTask(){
  MultiColvarBase* mcolv=getPntrToMultiColvar(); 
  // Copy the derivatives across
  mcolv->copyElementsToBridgedColvar( this );
  // Retrive the value of the multicolvar and apply filter
  double val=getElementValue(0), df, weight=applyFilter( val, df );

  // Now propegate derivatives
  if( !mcolv->weightHasDerivatives ){
     unsigned nstart=getNumberOfDerivatives(); setElementValue( 1, weight );
     for(unsigned i=0;i<mcolv->atoms_with_derivatives.getNumberActive();++i){
        unsigned n=mcolv->atoms_with_derivatives[i], nx=3*n;
        atoms_with_derivatives.activate(n);
        addElementDerivative( nstart+nx+0, df*getElementDerivative(nx+0) );
        addElementDerivative( nstart+nx+1, df*getElementDerivative(nx+1) );
        addElementDerivative( nstart+nx+2, df*getElementDerivative(nx+2) );
     }
     for(unsigned i=3*mcolv->getNumberOfAtoms();i<mcolv->getNumberOfDerivatives();++i){
        addElementDerivative( nstart+i, df*getElementDerivative(i) );
     }
  } else {
      unsigned nstart=getNumberOfDerivatives();
      double ww=mcolv->getElementValue(1); setElementValue( 1, ww*weight );
      for(unsigned i=0;i<mcolv->atoms_with_derivatives.getNumberActive();++i){
          unsigned n=mcolv->atoms_with_derivatives[i], nx=3*n;
          atoms_with_derivatives.activate(n);
          addElementDerivative( nstart+nx+0, weight*mcolv->getElementDerivative(nstart+nx+0) + ww*df*getElementDerivative(nx+0) );
          addElementDerivative( nstart+nx+1, weight*mcolv->getElementDerivative(nstart+nx+0) + ww*df*getElementDerivative(nx+1) );
          addElementDerivative( nstart+nx+2, weight*mcolv->getElementDerivative(nstart+nx+0) + ww*df*getElementDerivative(nx+2) );
     }
     for(unsigned i=3*mcolv->getNumberOfAtoms();i<mcolv->getNumberOfDerivatives();++i){
         addElementDerivative( nstart+i, weight*mcolv->getElementDerivative(nstart+i) + ww*df*getElementDerivative(i) );
     }
  }
}

void MultiColvarFilter::addBridgeForces( const std::vector<double>& bb ){ 
  plumed_dbg_assert( bb.size()==0 );
}

}
}
