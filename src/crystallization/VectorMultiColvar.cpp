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
#include "multicolvar/MultiColvarFunction.h"
#include "multicolvar/BridgedMultiColvarFunction.h"
#include "VectorMultiColvar.h"

namespace PLMD {
namespace crystallization {

void VectorMultiColvar::registerKeywords( Keywords& keys ){
  MultiColvar::registerKeywords( keys );
}

VectorMultiColvar::VectorMultiColvar(const ActionOptions& ao):
PLUMED_MULTICOLVAR_INIT(ao),
firstcall(false),
vecs(NULL)
{
  setLowMemOption(true);
}

void VectorMultiColvar::setVectorDimensionality( const unsigned& ncomp, const bool& comp, const int& nat ){
  // Store number of derivatives and if vectors are complex
  ncomponents = ncomp; complexvec=comp; 
  if(complexvec) dervec.resize( 2*ncomponents );
  else dervec.resize( ncomponents );
  // Read in the atoms if we are using multicolvar reading
  int natoms=nat; readAtoms( natoms );
  // Create the store vector object
  std::string param; vesselbase::VesselOptions da("","",0,param,this);
  Keywords keys; StoreVectorsVessel::registerKeywords( keys );
  vesselbase::VesselOptions da2(da,keys);
  vecs = new StoreVectorsVessel(da2);
  // Add the vessel to the base
  addVessel(vecs);
  // Read in any vessels
  readVesselKeywords();
  // Resize a holder for the derivatives of the norm of the vector
}

void VectorMultiColvar::doNotCalculateDirector(){
  vecs->store_director=false;
}

double VectorMultiColvar::doCalculation(){
  // Now calculate the vector
  calculateVector();
  // Sort out the active derivatives
  updateActiveAtoms();

  // Now calculate the norm of the vector (this is what we return here)
  double norm=0, inorm;
  if(complexvec){
     for(unsigned i=0;i<ncomponents;++i) norm += getComponent(i)*getComponent(i) + getImaginaryComponent(i)*getImaginaryComponent(i); 
     norm=sqrt(norm); inorm = 1.0 / norm;
     for(unsigned i=0;i<ncomponents;++i){ dervec[i] = inorm*getComponent(i); dervec[ncomponents+i] = inorm*getImaginaryComponent(i); } 
  } else {
     for(unsigned i=0;i<ncomponents;++i) norm += getComponent(i)*getComponent(i);
     norm=sqrt(norm); inorm = 1.0 / norm;
     for(unsigned i=0;i<ncomponents;++i) dervec[i] = inorm*getComponent(i); 
  }

  if( !doNotCalculateDerivatives() ){
      if( usingLowMem() ){
         vecs->storeDerivativesLowMem( 0 );
         vecs->chainRule( 0, dervec );
      } else {
         vecs->storeDerivativesHighMem( getCurrentPositionInTaskList() );
         vecs->chainRule( getCurrentPositionInTaskList(), dervec );
      }

      // Add derivatives to base multicolvars
      Vector tmpd;
      for(unsigned i=0;i<atoms_with_derivatives.getNumberActive();++i){
           unsigned k=atoms_with_derivatives[i];
           tmpd[0]=vecs->getFinalDerivative(3*i+0); 
           tmpd[1]=vecs->getFinalDerivative(3*i+1); 
           tmpd[2]=vecs->getFinalDerivative(3*i+2); 
           MultiColvarBase::addAtomsDerivatives( 0, k, tmpd );
      }   
      unsigned vvbase=3*atoms_with_derivatives.getNumberActive(); Tensor tmpv;
      for(unsigned i=0;i<3;++i){
          for(unsigned j=0;j<3;++j){
              tmpv(i,j) = vecs->getFinalDerivative( vvbase+3*i+j ); 
          }   
      }   
      MultiColvarBase::addBoxDerivatives( 0, tmpv );
  }
  
  return norm;
}

vesselbase::StoreDataVessel* VectorMultiColvar::buildDataStashes( const bool& allow_wcutoff, const double& wtol ){
  // Build everyting for the multicolvar
  vesselbase::StoreDataVessel* vsv=MultiColvarBase::buildDataStashes( allow_wcutoff, wtol );
  if( allow_wcutoff ) vsv->setHardCutoffOnWeight( wtol );
  // Resize the variable
  vecs->resize();
  // And make sure we set up the vector storage correctly
  vv1.resize( 1 ); vv2.resize( getNumberOfQuantities() - 5 );
  // And return
  return vsv;
}

void VectorMultiColvar::getValueForTask( const unsigned& iatom, std::vector<double>& vals ){
  plumed_dbg_assert( vecs && vals.size()==(getNumberOfQuantities()-4) ); 
  MultiColvarBase::getValueForTask( iatom, vv1 ); vecs->getVector( iatom, vv2 );
  vals[0]=vv1[0]; for(unsigned i=0;i<vv2.size();++i) vals[i+1]=vv2[i];
}

void VectorMultiColvar::addWeightedValueDerivatives( const unsigned& iatom, const unsigned& base_cv_no, const double& weight, multicolvar::MultiColvarFunction* func ){
  if( usingLowMem() ){
     vecs->recompute( iatom, 1 ); 
     for(unsigned j=0;j<getNumberOfQuantities()-5;++j) vecs->chainRuleForComponent( 1, j, 5+j, base_cv_no, weight, func );
  } else {
     for(unsigned j=0;j<getNumberOfQuantities()-5;++j) vecs->chainRuleForComponent( iatom, j, 5+j, base_cv_no, weight, func );
  }
}

void VectorMultiColvar::finishWeightedAverageCalculation( multicolvar::MultiColvarFunction* func ){
  // And calculate the norm of the vector
  double norm=0, inorm; std::vector<unsigned> tmpindices( 1 + func->getNumberOfDerivatives() );
  if(complexvec){
     for(unsigned i=0;i<ncomponents;++i){
        // Calculate average vector
        func->quotientRule(5+i, 1, 5+i); func->quotientRule(5+ncomponents+i, 1, 5+ncomponents+i);
        // Calculate length of vector
        norm += func->getElementValue(5+i)*func->getElementValue(5+i) + func->getElementValue(5+ncomponents+i)*func->getElementValue(5+ncomponents+i);
     }
     norm=sqrt(norm); inorm = 1.0 / norm;
     for(unsigned i=0;i<ncomponents;++i){ 
        dervec[i] = inorm*func->getElementValue(5+i); dervec[ncomponents+i] = inorm*func->getElementValue(5+ncomponents+i); 
     }
     func->getIndexList( 1, 0, func->getNumberOfDerivatives(), tmpindices );
     unsigned nder = func->getNumberOfDerivatives();
     for(unsigned i=0;i<tmpindices[0];++i){
         unsigned ind = tmpindices[1+i];
         for(unsigned j=0;j<ncomponents;++j){
             func->addElementDerivative( ind, dervec[j]*func->getElementDerivative(nder*(5+j) + ind) );
             func->addElementDerivative( ind, dervec[ncomponents+j]*func->getElementDerivative(nder*(5+ncomponents+j) + ind) );
         }
     }
  } else {
     for(unsigned i=0;i<ncomponents;++i){
         // Calculate average vector
         func->quotientRule(5+i, 1, 5+i);
         // Calculate length of vector
         norm += func->getElementValue(5+i)*func->getElementValue(5+i);
     }
     norm=sqrt(norm); inorm = 1.0 / norm;
     for(unsigned i=0;i<ncomponents;++i) dervec[i] = inorm*func->getElementValue(5+i); 
     func->getIndexList( 1, 0, func->getNumberOfDerivatives(), tmpindices );
     // And set derivatives given magnitude of the vector
     unsigned nder = func->getNumberOfDerivatives();
     for(unsigned i=0;i<tmpindices[0];++i){
         unsigned ind = tmpindices[1+i];
         for(unsigned j=0;j<ncomponents;++j){
             func->addElementDerivative( ind, dervec[j]*func->getElementDerivative(nder*(5+j) + ind) );
         }
     }
  }
  func->setElementValue( 0, norm );
}

void VectorMultiColvar::addOrientationDerivativesToBase( const unsigned& iatom, const unsigned& jstore, const unsigned& base_cv_no, 
                                                         const std::vector<double>& der, multicolvar::MultiColvarFunction* func ){
  if( usingLowMem() ){
      if(jstore==1){
         if(firstcall){ vecs->recompute( iatom, jstore ); firstcall=false; }
         vecs->chainRuleForVector( jstore, 0, base_cv_no, der, func );
      } else {
         vecs->recompute( iatom, jstore );
         vecs->chainRuleForVector( jstore, 0, base_cv_no, der, func );
      }
  } else {
      vecs->chainRuleForVector( iatom, 0, base_cv_no, der, func );
  }
}

void VectorMultiColvar::addForcesOnAtoms( const std::vector<double>& inforces ){
  plumed_dbg_assert( inforces.size()==getNumberOfDerivatives() );
  std::vector<double> oldforces( getNumberOfDerivatives() ); 
  getForcesFromVessels( oldforces ); 
  for(unsigned i=0;i<getNumberOfDerivatives();++i) oldforces[i]+=inforces[i];
  setForcesOnAtoms( oldforces );
}

void VectorMultiColvar::copyElementsToBridgedColvar( multicolvar::BridgedMultiColvarFunction* func ){
  MultiColvarBase::copyElementsToBridgedColvar( func );

  for(unsigned icomp=5;icomp<getNumberOfQuantities();++icomp){
      func->setElementValue( icomp-4, getElementValue(icomp) );
      unsigned nbase =  icomp * getNumberOfDerivatives();
      unsigned nbasev = (icomp-4) * func->getNumberOfDerivatives();
      for(unsigned jatom=0;jatom<atoms_with_derivatives.getNumberActive();++jatom){
          unsigned n=atoms_with_derivatives[jatom], nx=nbase + 3*n, ny=nbasev + 3*n;
          func->addElementDerivative( ny+0, getElementDerivative(nx+0) );
          func->addElementDerivative( ny+1, getElementDerivative(nx+1) );
          func->addElementDerivative( ny+2, getElementDerivative(nx+2) );
     }
     unsigned nwvir=nbase + 3*getNumberOfAtoms(), nwvirv=nbasev + 3*getNumberOfAtoms();
     for(unsigned k=0;k<9;++k){
        func->addElementDerivative( nwvirv, getElementDerivative(nwvir) ); nwvir++; nwvirv++;
     }
  }
}

}
}
