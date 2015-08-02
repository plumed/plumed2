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
MultiColvarBase(ao)
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

void AdjacencyMatrixBase::requestAtoms( const std::vector<AtomNumber>& atoms ){
  // Create the storeAdjacencyMatrixVessel
  std::string param; vesselbase::VesselOptions da("","",0,param,this);
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
