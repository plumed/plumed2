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
#include "VolumeGradientBase.h"

namespace PLMD {
namespace multicolvar {

void VolumeGradientBase::registerKeywords( Keywords& keys ){
  BridgedMultiColvarFunction::registerKeywords( keys );
}

VolumeGradientBase::VolumeGradientBase(const ActionOptions&ao):
Action(ao),
BridgedMultiColvarFunction(ao)
{
}

void VolumeGradientBase::requestAtoms( const std::vector<AtomNumber>& atoms ){
  ActionAtomistic::requestAtoms(atoms); bridgeVariable=3*atoms.size();
  addDependency( getPntrToMultiColvar() ); 
  tmpforces.resize( 3*atoms.size()+9 );
}

void VolumeGradientBase::doJobsRequiredBeforeTaskList(){
  ActionWithValue::clearDerivatives();
  retrieveAtoms(); setupRegions();
  ActionWithVessel::doJobsRequiredBeforeTaskList();
}

void VolumeGradientBase::completeTask(){
  if( getPntrToMultiColvar()->isDensity() ){ 
     setElementValue( 0, 1.0 );
  } else {
     // Copy derivatives of the colvar and the value of the colvar
     getPntrToMultiColvar()->copyElementsToBridgedColvar( this );
  }
  calculateAllVolumes(); 
}

void VolumeGradientBase::setNumberInVolume( const unsigned& ivol, const double& weight, const Vector& wdf ){
  MultiColvarBase* mcolv=getPntrToMultiColvar();
  if( !mcolv->weightHasDerivatives ){
      unsigned nstart=ivol*getNumberOfDerivatives();
      setElementValue( ivol, weight ); 
      for(unsigned i=0;i<mcolv->atomsWithCatomDer.getNumberActive();++i){
         unsigned n=mcolv->atomsWithCatomDer[i], nx=nstart + 3*n;
         atoms_with_derivatives.activate(n);
         addElementDerivative( nx+0, mcolv->getCentralAtomDerivative(n, 0, wdf ) );
         addElementDerivative( nx+1, mcolv->getCentralAtomDerivative(n, 1, wdf ) );
         addElementDerivative( nx+2, mcolv->getCentralAtomDerivative(n, 2, wdf ) );
      } 
  } else {
      unsigned nstart=ivol*getNumberOfDerivatives();
      double ww=mcolv->getElementValue(1); setElementValue( ivol, ww*weight );
      for(unsigned i=0;i<mcolv->atomsWithCatomDer.getNumberActive();++i){
          unsigned n=mcolv->atomsWithCatomDer[i], nx=nstart + 3*n;
          atoms_with_derivatives.activate(n);
          addElementDerivative( nx+0, ww*mcolv->getCentralAtomDerivative(n, 0, wdf ) );
          addElementDerivative( nx+1, ww*mcolv->getCentralAtomDerivative(n, 1, wdf ) );
          addElementDerivative( nx+2, ww*mcolv->getCentralAtomDerivative(n, 2, wdf ) );
     }
     unsigned nder=mcolv->getNumberOfDerivatives(); 
     for(unsigned i=0;i<mcolv->atoms_with_derivatives.getNumberActive();++i){
        unsigned n=mcolv->atoms_with_derivatives[i], nx=nder + 3*n, ny=nstart + 3*n;
        atoms_with_derivatives.activate(n);
        addElementDerivative( ny+0, weight*mcolv->getElementDerivative(nx+0) );
        addElementDerivative( ny+1, weight*mcolv->getElementDerivative(nx+1) );
        addElementDerivative( ny+2, weight*mcolv->getElementDerivative(nx+2) );
     }
     unsigned nwvir=mcolv->getNumberOfDerivatives()-9, nwvir2=nstart + 3*mcolv->getNumberOfAtoms();
     for(unsigned i=0;i<9;++i){
        addElementDerivative( nwvir2, weight*mcolv->getElementDerivative(nwvir) ); nwvir++; nwvir2++;
     }
  }
}

void VolumeGradientBase::addBridgeForces( const std::vector<double>& bb ){ 
  plumed_dbg_assert( bb.size()==tmpforces.size()-9 );
  // Forces on local atoms
  for(unsigned i=0;i<bb.size();++i) tmpforces[i]=bb[i];
  // Virial contribution is zero
  for(unsigned i=bb.size();i<bb.size()+9;++i) tmpforces[i]=0.0;
  setForcesOnAtoms( tmpforces, 0 );
}

}
}
