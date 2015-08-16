/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2015 The plumed team
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
}

AdjacencyMatrixBase::AdjacencyMatrixBase(const ActionOptions& ao):
Action(ao),
MultiColvarBase(ao),
connect_id(0)
{
  // Weight of this does have derivatives
  parse("WTOL",wtolerance);
  if(wtolerance>0) log.printf("  only considering those colvars with a weight greater than %f \n",wtolerance);
}

void AdjacencyMatrixBase::parseAtomList(const std::string& key, const int& num, const bool& isnodes, std::vector<AtomNumber>& t){
  std::string newkey; t.resize(0);
  if( num<0 ){  
     newkey=key;
  } else {
     std::string snum;
     Tools::convert( num, snum ); newkey=key+snum;
  }

  if( isnodes ){
     std::vector<std::string> mlabs; parseVector(newkey,mlabs);
     if( mlabs.size()==0 ) return;
     myinputdata.setup( mlabs, plumed.getActionSet(), wtolerance, this );
     log.printf("  using colvars calculated by actions "); 
     for(unsigned i=0;i<mlabs.size();++i) log.printf("%s ",mlabs[i].c_str() );
     log.printf("\n"); 
  } else {
     ActionAtomistic::parseAtomList( key, num, t );
  }
}

void AdjacencyMatrixBase::parseConnectionDescriptions( const std::string& key, const unsigned& nrow_t ){
  if( (getNumberOfNodeTypes()==1 && nrow_t==0) || (getNumberOfNodeTypes()==2 && nrow_t==1) ){
      std::string sw; parse(key,sw);
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
  unsigned nq = myinputdata.getBaseColvar(0)->getNumberOfQuantities();
  for(unsigned i=1;i<myinputdata.getNumberOfBaseMultiColvars();++i){
     if( myinputdata.getBaseColvar(i)->getNumberOfQuantities()!=nq ) error("mismatch between vectors in base colvars");
  }
  return nq;
}

unsigned AdjacencyMatrixBase::getNumberOfNodeTypes() const {
  return myinputdata.getNumberOfBaseMultiColvars();
}

unsigned AdjacencyMatrixBase::getNumberOfNodes() const {
  return myinputdata.getFullNumberOfBaseTasks();
}

void AdjacencyMatrixBase::requestAtoms( const std::vector<AtomNumber>& atoms, const bool& symmetric, const unsigned& nrows ){
  unsigned icoef, jcoef, kcoef, kcount;
  // Create the task list
  if( atoms.size()==0 ){
      if( symmetric || nrows==0 ){
          nblock=getNumberOfNodes();
      } else {
          nblock=nrows;
          if( (getNumberOfNodes()-nrows)>nblock ) nblock = getNumberOfNodes()-nrows; 
      }
      ablocks.resize(2); icoef=nblock; jcoef=1; kcoef=0; kcount=1;
  } else {
      if( symmetric || nrows==0 ){
          nblock=getNumberOfNodes(); 
      } else {
          nblock=nrows;
          if( (getNumberOfNodes()-nrows)>nblock ) nblock = getNumberOfNodes()-nrows;
      } 
      if( kcount>nblock ) nblock=kcount; 

      kcount=atoms.size(); ablocks.resize(3);
      icoef=nblock*nblock; jcoef=nblock; kcoef=1; ablocks[2].resize( atoms.size() );
      for(unsigned i=0;i<ablocks[2].size();++i) ablocks[2][i]=getNumberOfAtoms() - atoms.size() + i;
  }
  
  if( symmetric && nrows==0 ){ 
     plumed_dbg_assert( nrows==0 );
     resizeBookeepingArray( getNumberOfNodes(), getNumberOfNodes() );
     ablocks[0].resize( getNumberOfNodes() ); ablocks[1].resize( getNumberOfNodes() ); 
     for(unsigned i=0;i<getNumberOfNodes();++i) ablocks[0][i]=ablocks[1][i]=i;
     for(unsigned i=1;i<getNumberOfNodes();++i){
        for(unsigned j=0;j<i;++j){
           bookeeping(i,j).first=getFullNumberOfTasks();
           for(unsigned k=0;k<kcount;++k) addTaskToList( i*icoef + j*jcoef + k*kcoef );
           bookeeping(i,j).second=getFullNumberOfTasks();
        }
     }
  } else if( nrows==0 ){
     resizeBookeepingArray( getNumberOfNodes(), getNumberOfNodes() );
     ablocks[0].resize( getNumberOfNodes() ); ablocks[1].resize( getNumberOfNodes() );
     for(unsigned i=0;i<getNumberOfNodes();++i) ablocks[0][i]=ablocks[1][i]=i;
     for(unsigned i=0;i<getNumberOfNodes();++i){
        for(unsigned j=0;j<getNumberOfNodes();++j){
           bookeeping(i,j).first=getFullNumberOfTasks();
           for(unsigned k=0;k<kcount;++k) addTaskToList( i*icoef + j*jcoef + k*kcoef );
           bookeeping(i,j).second=getFullNumberOfTasks();
        }
     }
  } else {
     unsigned nto = getNumberOfNodes() - nrows;
     resizeBookeepingArray( nrows, nto ); ablocks[0].resize( nrows ); ablocks[1].resize( nto );
     for(unsigned i=0;i<nrows;++i) ablocks[0][i]=i;
     for(unsigned i=0;i<nto;++i) ablocks[1][i] = nrows + i;
     for(unsigned i=0;i<nrows;++i){
         for(unsigned j=0;j<nto;++j){
             bookeeping(i,j).first=getFullNumberOfTasks();
             for(unsigned k=0;k<kcount;++k) addTaskToList( i*icoef + j*jcoef + k*kcoef );   
             bookeeping(i,j).second=getFullNumberOfTasks();
         }
     }
  }

  // Create the storeAdjacencyMatrixVessel
  std::string nr,ncols,param;
  Tools::convert( ablocks[0].size(), nr );
  Tools::convert( ablocks[1].size(), ncols ); 
  param = "NROWS=" + nr + " NCOLS=" + ncols;
  if( symmetric && nrows==0 ) param+=" SYMMETRIC"; 
  if( !symmetric && nrows==0 ) param+=" HBONDS";
 
  vesselbase::VesselOptions da("","",0,param,this);
  Keywords keys; AdjacencyMatrixVessel::registerKeywords( keys );
  vesselbase::VesselOptions da2(da,keys);
  mat = new AdjacencyMatrixVessel(da2);
  // Set a cutoff for clustering
  mat->setHardCutoffOnWeight( getTolerance() );
  // Add the vessel to the base
  addVessel( mat );
  // Request the data required
  myinputdata.makeDataRequests( atoms, this ); 
  setupMultiColvarBase();
}

void AdjacencyMatrixBase::calculate(){
  if( checkNumericalDerivatives() ) error("numerical derivatives currently broken");
  // Setup the linke cells
  setupLinkCells();
  // And run all tasks
  runAllTasks();
}

Vector AdjacencyMatrixBase::getPositionOfAtomForLinkCells( const unsigned& iatom ) const {
  if( iatom>=myinputdata.getFullNumberOfBaseTasks() ) return getPosition( iatom );
  return myinputdata.getPosition( iatom );
}

void AdjacencyMatrixBase::updateActiveAtoms( multicolvar::AtomValuePack& myatoms ) const {
  if( !doNotCalculateDerivatives() ) myatoms.updateDynamicList();
}

bool AdjacencyMatrixBase::isCurrentlyActive( const unsigned& bno, const unsigned& code ){
  if( code>=myinputdata.getFullNumberOfBaseTasks() ) return true; 
  return myinputdata.isCurrentlyActive( 0, code );
}

void AdjacencyMatrixBase::addAtomDerivatives( const unsigned& iatom, const Vector& der, multicolvar::AtomValuePack& myatoms ) const {
  unsigned jatom=myatoms.getIndex(iatom);

  if( jatom>myinputdata.getFullNumberOfBaseTasks() ){
      myatoms.addAtomsDerivatives( 1, jatom, der );
  } else {
      myinputdata.addComDerivatives( jatom, der, myatoms );
  }
}

void AdjacencyMatrixBase::addOrientationDerivatives( const unsigned& iatom, const std::vector<double>& der, multicolvar::AtomValuePack& myatoms ) const {
  unsigned jatom=myatoms.getIndex(iatom); plumed_dbg_assert( jatom<myinputdata.getFullNumberOfBaseTasks() );
  MultiValue myder(0,0); myinputdata.getVectorDerivatives( jatom, true, myder );
  myinputdata.mergeVectorDerivatives( 1, 2, der.size(), jatom, der, myder, myatoms );
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
