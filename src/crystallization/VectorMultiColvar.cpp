/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2017 The plumed team
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
#include "VectorMultiColvar.h"
#include "multicolvar/BridgedMultiColvarFunction.h"

namespace PLMD {
namespace crystallization {

void VectorMultiColvar::registerKeywords( Keywords& keys ) {
  MultiColvarBase::registerKeywords( keys );
}

VectorMultiColvar::VectorMultiColvar(const ActionOptions& ao):
  Action(ao),
  MultiColvarBase(ao),
  store_director(true)
{
  setLowMemOption(true);
}

void VectorMultiColvar::setVectorDimensionality( const unsigned& ncomp ) {
  ncomponents = ncomp; resizeFunctions(); // This resize needs to be here to ensure buffers are set to correct size in base
}

void VectorMultiColvar::doNotCalculateDirector() {
  store_director=false;    // Need a sanity check in here  so that you don't use the same instance of Q4 to calcualte vectors and directors
}

double VectorMultiColvar::compute( const unsigned& taskIndex, multicolvar::AtomValuePack& myatoms ) const {
  // Now calculate the vector
  calculateVector( myatoms );
  // Sort out the active derivatives
  updateActiveAtoms( myatoms );

  // Now calculate the norm of the vector (this is what we return here)
  double norm=0;
  for(unsigned i=0; i<ncomponents; ++i) norm += myatoms.getValue(2+i)*myatoms.getValue(2+i);
  norm=sqrt(norm);

  if( !doNotCalculateDerivatives() ) {
    double inorm = 1.0 / norm; std::vector<double> dervec( ncomponents );
    for(unsigned i=0; i<ncomponents; ++i) dervec[i] = inorm*myatoms.getValue(2+i);

    MultiValue& myvals=myatoms.getUnderlyingMultiValue();
    for(unsigned j=0; j<myvals.getNumberActive(); ++j) {
      unsigned jder=myvals.getActiveIndex(j);
      for(unsigned i=0; i<ncomponents; ++i) myvals.addDerivative( 1, jder, dervec[i]*myvals.getDerivative( 2+i, jder ) );
    }
  }

  return norm;
}

void VectorMultiColvar::normalizeVector( std::vector<double>& vals ) const {
  double inorm = 1.0;
  if( vals[1]>epsilon ) inorm = 1.0 / vals[1];
  for(unsigned i=2; i<vals.size(); ++i) vals[i] = inorm*vals[i];
}

void VectorMultiColvar::normalizeVectorDerivatives( MultiValue& myvals ) const {
  double v = myvals.get(1), weight = 1.0 / v,  wdf = 1.0 / ( v*v*v );
  for(unsigned j=0; j<myvals.getNumberActive(); ++j) {
    double comp2=0.0; unsigned jder=myvals.getActiveIndex(j);
    for(unsigned jcomp=2; jcomp<myvals.getNumberOfValues(); ++jcomp) comp2 += myvals.get(jcomp)*myvals.getDerivative( jcomp, jder );
    for(unsigned jcomp=2; jcomp<myvals.getNumberOfValues(); ++jcomp) {
      myvals.setDerivative( jcomp, jder, weight*myvals.getDerivative( jcomp, jder ) - wdf*comp2*myvals.get(jcomp) );
    }
  }
}

void VectorMultiColvar::addForcesOnAtoms( const std::vector<double>& inforces ) {
  plumed_dbg_assert( inforces.size()==getNumberOfDerivatives() );
  std::vector<double> oldforces( getNumberOfDerivatives() );
  getForcesFromVessels( oldforces );
  for(unsigned i=0; i<getNumberOfDerivatives(); ++i) oldforces[i]+=inforces[i];
  setForcesOnAtoms( oldforces );
}

}
}
