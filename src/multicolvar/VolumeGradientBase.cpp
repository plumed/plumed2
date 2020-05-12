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
#include "VolumeGradientBase.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "CatomPack.h"

namespace PLMD {
namespace multicolvar {

void VolumeGradientBase::registerKeywords( Keywords& keys ) {
  BridgedMultiColvarFunction::registerKeywords( keys );
}

VolumeGradientBase::VolumeGradientBase(const ActionOptions&ao):
  Action(ao),
  BridgedMultiColvarFunction(ao)
{
}

void VolumeGradientBase::requestAtoms( const std::vector<AtomNumber>& atoms ) {
  ActionAtomistic::requestAtoms(atoms); bridgeVariable=3*atoms.size();
  std::map<std::string,bool> checklabs;
  for(const auto & p : getDependencies() ) checklabs.insert(std::pair<std::string,bool>(p->getLabel(),false));
  for(const auto & p : plumed.getActionSet() ) {
    if( p->getLabel()==getPntrToMultiColvar()->getLabel() ) break;
    if( checklabs.count(p->getLabel()) ) checklabs[p->getLabel()]=true;
  }
  for(const auto & p : checklabs ) {
    if( !p.second ) error("the input for the virtual atoms used in the input for this action must appear in the input file before the input multicolvar");
  }
  addDependency( getPntrToMultiColvar() );
  tmpforces.resize( 3*atoms.size()+9 );
}

void VolumeGradientBase::doJobsRequiredBeforeTaskList() {
  ActionWithValue::clearDerivatives();
  retrieveAtoms(); setupRegions();
  ActionWithVessel::doJobsRequiredBeforeTaskList();
}

void VolumeGradientBase::completeTask( const unsigned& curr, MultiValue& invals, MultiValue& outvals ) const {
  if( getPntrToMultiColvar()->isDensity() ) {
    outvals.setValue(0, 1.0); outvals.setValue(1, 1.0);
  } else {
    // Copy derivatives of the colvar and the value of the colvar
    invals.copyValues( outvals );
    if( derivativesAreRequired() ) invals.copyDerivatives( outvals );
  }
  calculateAllVolumes( curr, outvals );
}

void VolumeGradientBase::setNumberInVolume( const unsigned& ivol, const unsigned& curr, const double& weight,
    const Vector& wdf, const Tensor& virial, const std::vector<Vector>& refders,
    MultiValue& outvals ) const {
  MultiColvarBase* mcolv=getPntrToMultiColvar();
  if( !mcolv->weightHasDerivatives ) {
    outvals.setValue(ivol, weight );
    if( derivativesAreRequired() ) {
      CatomPack catom; mcolv->getCentralAtomPack( 0, curr, catom );
      for(unsigned i=0; i<catom.getNumberOfAtomsWithDerivatives(); ++i) {
        unsigned jatom=3*catom.getIndex(i);
        outvals.addDerivative( ivol, jatom+0, catom.getDerivative(i,0,wdf) );
        outvals.addDerivative( ivol, jatom+1, catom.getDerivative(i,1,wdf) );
        outvals.addDerivative( ivol, jatom+2, catom.getDerivative(i,2,wdf) );
      }
      unsigned nmder=getPntrToMultiColvar()->getNumberOfDerivatives();
      for(unsigned i=0; i<3; ++i) for(unsigned j=0; j<3; ++j) outvals.addDerivative( ivol, nmder-9+3*i+j, virial(i,j) );
      for(unsigned i=0; i<refders.size(); ++i) {
        unsigned iatom=nmder+3*i;

        outvals.addDerivative( ivol, iatom+0, refders[i][0] );
        outvals.addDerivative( ivol, iatom+1, refders[i][1] );
        outvals.addDerivative( ivol, iatom+2, refders[i][2] );
      }
    }
  } else if(ivol==0) {
    double ww=outvals.get(0); outvals.setValue(ivol,ww*weight);
    if( derivativesAreRequired() ) {
      plumed_merror("This needs testing");
      CatomPack catom; mcolv->getCentralAtomPack( 0, curr, catom );
      for(unsigned i=0; i<catom.getNumberOfAtomsWithDerivatives(); ++i) {
        unsigned jatom=3*catom.getIndex(i);
        outvals.addDerivative( ivol, jatom+0, weight*outvals.getDerivative(ivol,jatom+0) + ww*catom.getDerivative(i,0,wdf) );
        outvals.addDerivative( ivol, jatom+1, weight*outvals.getDerivative(ivol,jatom+1) + ww*catom.getDerivative(i,1,wdf) );
        outvals.addDerivative( ivol, jatom+2, weight*outvals.getDerivative(ivol,jatom+2) + ww*catom.getDerivative(i,2,wdf) );
      }
      unsigned nmder=getPntrToMultiColvar()->getNumberOfDerivatives();
      for(unsigned i=0; i<3; ++i) for(unsigned j=0; j<3; ++j) outvals.addDerivative( ivol, nmder-9+3*i+j, ww*virial(i,j) );
      for(unsigned i=0; i<refders.size(); ++i) {
        unsigned iatom=nmder+3*i;
        outvals.addDerivative( ivol, iatom+0, ww*refders[i][0] );
        outvals.addDerivative( ivol, iatom+1, ww*refders[i][1] );
        outvals.addDerivative( ivol, iatom+2, ww*refders[i][2] );
      }
    }
  } else {
    double ww=outvals.get(0); outvals.setValue(ivol,ww*weight);
    if( derivativesAreRequired() ) {
      plumed_merror("This needs testing");
      CatomPack catom; mcolv->getCentralAtomPack( 0, curr, catom );
      for(unsigned i=0; i<catom.getNumberOfAtomsWithDerivatives(); ++i) {
        unsigned jatom=3*catom.getIndex(i);
        outvals.addDerivative( ivol, jatom+0, ww*catom.getDerivative(i,0,wdf) );
        outvals.addDerivative( ivol, jatom+1, ww*catom.getDerivative(i,1,wdf) );
        outvals.addDerivative( ivol, jatom+2, ww*catom.getDerivative(i,2,wdf) );
      }
      unsigned nmder=getPntrToMultiColvar()->getNumberOfDerivatives();
      for(unsigned i=0; i<3; ++i) for(unsigned j=0; j<3; ++j) outvals.addDerivative( ivol, nmder-9+3*i+j, ww*virial(i,j) );
      for(unsigned i=0; i<refders.size(); ++i) {
        unsigned iatom=nmder+3*i;
        outvals.addDerivative( ivol, iatom+0, ww*refders[i][0] );
        outvals.addDerivative( ivol, iatom+1, ww*refders[i][1] );
        outvals.addDerivative( ivol, iatom+2, ww*refders[i][2] );
      }
    }
  }
}

void VolumeGradientBase::addBridgeForces( const std::vector<double>& bb ) {
  plumed_dbg_assert( bb.size()==tmpforces.size()-9 );
  // Forces on local atoms
  for(unsigned i=0; i<bb.size(); ++i) tmpforces[i]=bb[i];
  // Virial contribution is zero
  for(unsigned i=bb.size(); i<bb.size()+9; ++i) tmpforces[i]=0.0;
  setForcesOnAtoms( tmpforces, 0 );
}

}
}
