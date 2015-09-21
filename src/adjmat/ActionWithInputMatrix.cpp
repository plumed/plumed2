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
#include "ActionWithInputMatrix.h"
#include "AdjacencyMatrixVessel.h"
#include "AdjacencyMatrixBase.h"
#include "MatrixSummationBase.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"

namespace PLMD {
namespace adjmat {

void ActionWithInputMatrix::registerKeywords( Keywords& keys ){
  Action::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys );
  ActionWithVessel::registerKeywords( keys );
  keys.add("compulsory","MATRIX","the action that calcualtes the adjacency matrix vessel we would like to analyse"); 
  keys.addFlag("NOPBC",false,"don't use periodic boundary conditions when calculating separations");
}


ActionWithInputMatrix::ActionWithInputMatrix(const ActionOptions& ao):
Action(ao),
ActionWithValue(ao),
ActionAtomistic(ao),
vesselbase::ActionWithVessel(ao),
usepbc(true),
mymatrix(NULL)
{
  bool nopbc; parseFlag("NOPBC",nopbc); usepbc=!nopbc;
  // Find the object that calculates our adjacency matrix
  std::vector<std::string> matname(1); parse("MATRIX",matname[0]);
  ActionWithVessel* myvess = plumed.getActionSet().selectWithLabel<ActionWithVessel*>( matname[0] );
  if( !myvess ) error( matname[0] + " does not calculate an adjacency matrix");

  // Retrieve the adjacency matrix of interest
  for(unsigned i=0;i<myvess->getNumberOfVessels();++i){
      mymatrix = dynamic_cast<AdjacencyMatrixVessel*>( myvess->getPntrToVessel(i) );
      if( mymatrix ) break ;
  }
  if( !mymatrix ){
     MatrixSummationBase* mybase = dynamic_cast<MatrixSummationBase*>( myvess );
     if( !mybase ) error( matname[0] + " does not calculate an adjacency matrix");
     mymatrix = mybase->mymatrix; 
     // Now setup a data stash in the matrix summation object
     // myinputdata.setup( matname, plumed.getActionSet(), (mymatrix->function)->wtolerance, this );
  }
  log.printf("  using matrix calculated by action %s \n",(mymatrix->function)->getLabel().c_str() );

  // And get the atom requests
  ActionAtomistic* matoms = dynamic_cast<ActionAtomistic*>( myvess );
  plumed_assert( matoms ); requestAtoms( matoms->getAbsoluteIndexes() );
  // And add the dependency after the atom requst ( atom request resets dependences )
  addDependency( myvess );
}

void ActionWithInputMatrix::turnOnDerivatives(){
  ActionWithValue::turnOnDerivatives();
  needsDerivatives();
  forcesToApply.resize( getNumberOfDerivatives() );
} 

unsigned ActionWithInputMatrix::getNumberOfDerivatives() {
  return (mymatrix->function)->getNumberOfDerivatives();
}

unsigned ActionWithInputMatrix::getNumberOfNodes() const {
  return (mymatrix->function)->ablocks[0].size(); 
}

AdjacencyMatrixVessel* ActionWithInputMatrix::getAdjacencyVessel() const {
  return mymatrix;
}

AtomNumber ActionWithInputMatrix::getAbsoluteIndexOfCentralAtom(const unsigned& i) const {
  return (mymatrix->function)->getAbsoluteIndexOfCentralAtom(i);  
}

Vector ActionWithInputMatrix::getPosition( const unsigned& iatom ) const {
  return (mymatrix->function)->getPositionOfAtomForLinkCells(iatom); 
}

bool ActionWithInputMatrix::isCurrentlyActive( const unsigned& ind ) const {
  return (mymatrix->function)->isCurrentlyActive( 0, ind );  
}

void ActionWithInputMatrix::getVectorForTask( const unsigned& ind, const bool& normed, std::vector<double>& orient0 ) const {
  plumed_dbg_assert( isCurrentlyActive( ind ) ); 
  plumed_dbg_assert( ind<(mymatrix->function)->colvar_label.size() ); unsigned mmc=(mymatrix->function)->colvar_label[ind];
  plumed_dbg_assert( ((mymatrix->function)->mybasedata[mmc])->storedValueIsActive( (mymatrix->function)->convertToLocalIndex(ind,mmc) ) );
  ((mymatrix->function)->mybasedata[mmc])->retrieveValue( (mymatrix->function)->convertToLocalIndex(ind,mmc), normed, orient0 );
}

void ActionWithInputMatrix::getVectorDerivatives( const unsigned& ind, const bool& normed, MultiValue& myder ) const {
  plumed_dbg_assert( isCurrentlyActive( ind ) ); 
  plumed_dbg_assert( ind<(mymatrix->function)->colvar_label.size() ); unsigned mmc=(mymatrix->function)->colvar_label[ind];
  plumed_dbg_assert( ((mymatrix->function)->mybasedata[mmc])->storedValueIsActive( (mymatrix->function)->convertToLocalIndex(ind,mmc) ) );
  if( myder.getNumberOfValues()!=(mymatrix->function)->mybasemulticolvars[mmc]->getNumberOfQuantities() ||
      myder.getNumberOfDerivatives()!=(mymatrix->function)->mybasemulticolvars[mmc]->getNumberOfDerivatives() ){
          myder.resize( (mymatrix->function)->mybasemulticolvars[mmc]->getNumberOfQuantities(), (mymatrix->function)->mybasemulticolvars[mmc]->getNumberOfDerivatives() );
  }
  (mymatrix->function)->mybasedata[mmc]->retrieveDerivatives( (mymatrix->function)->convertToLocalIndex(ind,mmc), normed, myder );
}

Vector ActionWithInputMatrix::getSeparation( const Vector& vec1, const Vector& vec2 ) const {
  if(usepbc) return pbcDistance( vec1, vec2 );
  return delta( vec1, vec2 );
}

unsigned ActionWithInputMatrix::getNumberOfNodeTypes() const {
  unsigned size = (mymatrix->function)->mybasemulticolvars.size();
  if( size==0 ) return 1;
  return size; 
}

unsigned ActionWithInputMatrix::getNumberOfAtomsInGroup( const unsigned& igrp ) const {
 plumed_dbg_assert( igrp<(mymatrix->function)->mybasemulticolvars.size() );
 return (mymatrix->function)->mybasemulticolvars[igrp]->getFullNumberOfTasks(); 
}

multicolvar::MultiColvarBase* ActionWithInputMatrix::getBaseMultiColvar( const unsigned& igrp ) const {
 return (mymatrix->function)->mybasemulticolvars[igrp];
}

void ActionWithInputMatrix::apply(){
 if( getForcesFromVessels( forcesToApply ) ) setForcesOnAtoms( forcesToApply ); 
}

}
}
