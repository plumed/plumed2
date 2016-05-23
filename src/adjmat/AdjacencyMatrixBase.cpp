/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015,2016 The plumed team
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
#include "AdjacencyMatrixBase.h"
#include "multicolvar/BridgedMultiColvarFunction.h"
#include "multicolvar/AtomValuePack.h"
#include "multicolvar/CatomPack.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"

namespace PLMD {
namespace adjmat {

void AdjacencyMatrixBase::registerKeywords( Keywords& keys ){
  multicolvar::MultiColvarBase::registerKeywords( keys );
  keys.add("compulsory","WTOL","0.0","if the base multicolvars have weights then you must define a hard cutoff on those you want to consider explicitally");
  keys.remove("LOWMEM"); keys.use("HIGHMEM");
}

AdjacencyMatrixBase::AdjacencyMatrixBase(const ActionOptions& ao):
Action(ao),
MultiColvarBase(ao),
connect_id(0),
no_third_dim_accum(true)
{
  // Weight of this does have derivatives
  parse("WTOL",wtolerance);
  if(wtolerance>0) log.printf("  only considering those colvars with a weight greater than %f \n",wtolerance);
}

bool AdjacencyMatrixBase::parseAtomList(const std::string& key, const int& num, std::vector<AtomNumber>& t){
  std::vector<std::string> mlabs; 
  if( num<0 ) parseVector(key,mlabs);
  else parseNumberedVector(key,num,mlabs);

  if( mlabs.size()==0 ) return false;

  bool found_acts=interpretInputMultiColvars(mlabs,0.0);
  if( !found_acts ){
     ActionAtomistic::interpretAtomList( mlabs, t );
     log.printf("  involving atoms ");
     for(unsigned i=0;i<t.size();++i) log.printf("%d ",t[i].serial() );
     log.printf("\n");
  }
  return true;
}

void AdjacencyMatrixBase::parseConnectionDescriptions( const std::string& key, const unsigned& nrow_t ){
  if( getNumberOfNodeTypes()==1 || (getNumberOfNodeTypes()==2 && nrow_t==1) ){
      std::string sw; parse(key,sw);
      if(sw.length()==0) error("could not find " + key + " keyword");
      setupConnector( connect_id, 0, 0, sw );
  } else {
      unsigned nr, nc;
      if( nrow_t==0 ){
        nr=nc=getNumberOfNodeTypes();
      } else{
        nr=nrow_t; nc = getNumberOfNodeTypes() - nr;
      }
      for(unsigned i=0;i<nr;++i){
          // Retrieve the base number  
          unsigned ibase;
          if( nc<10 ){
             ibase=(i+1)*10;
          } else if ( nc<100 ){
             ibase=(i+1)*100;
          } else {
             error("wow this is an error I never would have expected");
          }

          for(unsigned j=i;j<nc;++j){
             std::string sw; parseNumbered(key,ibase+j+1,sw);
             if(sw.length()==0){
                std::string num; Tools::convert(ibase+j+1,num);
                error("could not find " + key + num + " keyword. Need one " + key + " keyword for each distinct base-multicolvar-pair type");
             }
             setupConnector( connect_id, i, j, sw );
          }
      }
  }
  connect_id++;
}

unsigned AdjacencyMatrixBase::getSizeOfInputVectors() const {
  if( mybasemulticolvars.size()==0 ) return 2; 

  unsigned nq = mybasemulticolvars[0]->getNumberOfQuantities();
  for(unsigned i=1;i<mybasemulticolvars.size();++i){
     if( mybasemulticolvars[i]->getNumberOfQuantities()!=nq ) error("mismatch between vectors in base colvars");
  }
  return nq;
}

unsigned AdjacencyMatrixBase::getNumberOfNodeTypes() const {
  unsigned size=mybasemulticolvars.size();
  if( size==0 ) return 1;
  return size;
}

void AdjacencyMatrixBase::requestAtoms( const std::vector<AtomNumber>& atoms, const bool& symmetric, const bool& true_square, const std::vector<unsigned>& dims ){
  unsigned icoef, jcoef, kcoef, kcount;
  // Create the task list
  ablocks.resize( dims.size() ); nblock=dims[0]; 
  for(unsigned i=1;i<dims.size();++i){
     if( dims[i]>nblock ) nblock=dims[i];
  }
  if( dims.size()==2 ){ 
     icoef=nblock; jcoef=1; kcoef=0; kcount=1;
  } else if( dims.size()==3 ){
     if( no_third_dim_accum ){
        icoef=nblock; jcoef=1; kcoef=0; kcount=1;  
     } else{
        icoef=nblock*nblock; jcoef=nblock; kcoef=1; kcount=dims[2];
     }
     ablocks[2].resize( dims[2] );
     if(symmetric || true_square ) {
        for(unsigned i=0;i<ablocks[2].size();++i) ablocks[2][i]=dims[0]+i;
     } else {
        for(unsigned i=0;i<ablocks[2].size();++i) ablocks[2][i]=dims[0]+dims[1]+i;
     }
  } else {
     plumed_error();
  }
  if( (symmetric || true_square) ){
     plumed_assert( !(symmetric && true_square) && dims[0]==dims[1] );
     resizeBookeepingArray( dims[0], dims[0] );
     ablocks[0].resize( dims[0] ); ablocks[1].resize( dims[0] );
     for(unsigned i=0;i<dims[0];++i) ablocks[0][i]=ablocks[1][i]=i; 
  } else {
     resizeBookeepingArray( dims[0], dims[1] ); ablocks[0].resize( dims[0] ); ablocks[1].resize( dims[1] );
     for(unsigned i=0;i<dims[0];++i) ablocks[0][i]=i;
     for(unsigned i=0;i<dims[1];++i) ablocks[1][i]=dims[0]+i;
  }
  if( symmetric ){
    for(unsigned i=1;i<dims[0];++i){
      for(unsigned j=0;j<i;++j){
        bookeeping(j,i).first=bookeeping(i,j).first=getFullNumberOfTasks();
        for(unsigned k=0;k<kcount;++k) addTaskToList( i*icoef + j*jcoef + k*kcoef );
        bookeeping(j,i).second=bookeeping(i,j).second=getFullNumberOfTasks();
      }
    }
  } else {
    for(unsigned i=0;i<dims[0];++i){
       for(unsigned j=0;j<dims[1];++j){
           bookeeping(i,j).first=getFullNumberOfTasks();
           for(unsigned k=0;k<kcount;++k) addTaskToList( i*icoef + j*jcoef + k*kcoef );
           bookeeping(i,j).second=getFullNumberOfTasks();
       }
    }
  }

  // Create the storeAdjacencyMatrixVessel
  std::string param;
  if( symmetric && dims[0]==dims[1] ) param="SYMMETRIC"; 
  if( !symmetric && dims[0]==dims[1] ) param="HBONDS";
 
  vesselbase::VesselOptions da("","",0,param,this);
  Keywords keys; AdjacencyMatrixVessel::registerKeywords( keys );
  vesselbase::VesselOptions da2(da,keys);
  mat = new AdjacencyMatrixVessel(da2);
  // Set a cutoff for clustering
  mat->setHardCutoffOnWeight( getTolerance() );
  // Add the vessel to the base
  addVessel( mat );
  setupMultiColvarBase( atoms );
}

// Maybe put this back GAT to check that it is returning an atom number that is one of the nodes
// and not a hydrogen if we are doing HBPAMM
// AtomNumber AdjacencyMatrixBase::getAbsoluteIndexOfCentralAtom(const unsigned& i) const {
//   plumed_dbg_assert( i<myinputdata.getFullNumberOfBaseTasks() );
//   return myinputdata.getAtomicIndex( i );
// } 

void AdjacencyMatrixBase::addOrientationDerivatives( const unsigned& ival, const unsigned& iatom, const std::vector<double>& der, multicolvar::AtomValuePack& myatoms ) const {
  unsigned jatom=myatoms.getIndex(iatom); plumed_dbg_assert( jatom<colvar_label.size() );
  MultiValue myder(0,0); unsigned mmc=colvar_label[ival]; plumed_assert( !mybasemulticolvars[mmc]->weightWithDerivatives() );
  plumed_dbg_assert( mybasedata[mmc]->storedValueIsActive( convertToLocalIndex(ival,mmc) ) );
  if( myder.getNumberOfValues()!=mybasemulticolvars[mmc]->getNumberOfQuantities() ||
      myder.getNumberOfDerivatives()!=mybasemulticolvars[mmc]->getNumberOfDerivatives() ){
          myder.resize( mybasemulticolvars[mmc]->getNumberOfQuantities(), mybasemulticolvars[mmc]->getNumberOfDerivatives() );
  }
  mybasedata[mmc]->retrieveDerivatives( convertToLocalIndex(ival,mmc), true, myder );

  // Get start of indices for this atom
  unsigned basen=0; for(unsigned i=0;i<mmc;++i) basen+=3*mybasemulticolvars[i]->getNumberOfAtoms();

  MultiValue& myvals=myatoms.getUnderlyingMultiValue();
  // Now get the start of the virial
  unsigned virbas = myvals.getNumberOfDerivatives()-9;
  for(unsigned j=0;j<myder.getNumberActive();++j){
     unsigned jder=myder.getActiveIndex(j);
     if( jder<3*mybasemulticolvars[mmc]->getNumberOfAtoms() ){
         unsigned kder=basen+jder;
         for(unsigned icomp=2;icomp<der.size();++icomp){
             myvals.addDerivative( ival, kder, der[icomp]*myder.getDerivative( icomp, jder ) );
         }
     } else {
         unsigned kder=virbas + (jder - 3*mybasemulticolvars[mmc]->getNumberOfAtoms());
         for(unsigned icomp=2;icomp<der.size();++icomp){
             myvals.addDerivative( ival, kder, der[icomp]*myder.getDerivative( icomp, jder ) );
         }
     }
  }
}

void AdjacencyMatrixBase::recalculateMatrixElement( const unsigned& myelem, MultiValue& myvals ){
  std::vector<unsigned> myatoms; decodeIndexToAtoms( getTaskCode( myelem ), myatoms );
  unsigned i=myatoms[0], j=myatoms[1];
  for(unsigned k=bookeeping(i,j).first;k<bookeeping(i,j).second;++k){
      if( !taskIsCurrentlyActive(k) ) continue;
      performTask( k, getTaskCode(k), myvals );  // This may not accumulate as we would like  GAT
  }
}

}
}
