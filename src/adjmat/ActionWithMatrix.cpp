/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2023 The plumed team
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
#include "ActionWithMatrix.h"

namespace PLMD {
namespace adjmat {

void ActionWithMatrix::registerKeywords( Keywords& keys ) {
  ActionWithVector::registerKeywords( keys ); keys.use("ARG");
}

ActionWithMatrix::ActionWithMatrix(const ActionOptions&ao):
  Action(ao),
  ActionWithVector(ao),
  matrix_to_do_after(NULL),
  doInnerLoop(false)
{
}

void ActionWithMatrix::setupStreamedComponents( unsigned& nquants, unsigned& nmat, unsigned& maxcol, unsigned& nbookeeping ) {
  ActionWithVector::setupStreamedComponents( nquants, nmat, maxcol, nbookeeping );

  for(int i=0; i<getNumberOfComponents(); ++i) {
      Value* myval=getPntrToComponent(i);
      if( myval->getRank()!=2 || myval->hasDerivatives() ) continue;
      myval->setPositionInMatrixStash(nmat); nmat++; 
      if( !myval->valueIsStored() ) continue;
      if( myval->getShape()[1]>maxcol ) maxcol=myval->getShape()[1];
      myval->setMatrixBookeepingStart(nbookeeping);
      nbookeeping += myval->getShape()[0]*( 1 + getNumberOfColumns() );
      myval->reshapeMatrixStore( getNumberOfColumns() );
  }
}

void ActionWithMatrix::finishChainBuild( ActionWithVector* act ) {
   ActionWithMatrix* am=dynamic_cast<ActionWithMatrix*>(act); if(!am) return;
   if( matrix_to_do_after ) matrix_to_do_after->finishChainBuild( act );
   else matrix_to_do_after=am;
}

void ActionWithMatrix::getTotalMatrixBookeeping( unsigned& nbookeeping ) const {
  for(int i=0; i<getNumberOfComponents(); ++i) {
      const Value* myval=getConstPntrToComponent(i);
      if( myval->getRank()!=2 || myval->hasDerivatives() || !myval->valueIsStored() ) continue;
      nbookeeping += myval->getShape()[0]*( 1 + getNumberOfColumns() );
  }
  if( matrix_to_do_after ) matrix_to_do_after->getTotalMatrixBookeeping( nbookeeping );
}

void ActionWithMatrix::calculate() {
  // Update all the neighbour lists
  updateAllNeighbourLists();
  // Setup the matrix indices
  unsigned nbookeeping=0; getTotalMatrixBookeeping( nbookeeping ); 
  if( matrix_bookeeping.size()!=nbookeeping ) matrix_bookeeping.resize( nbookeeping );
  std::fill( matrix_bookeeping.begin(), matrix_bookeeping.end(), 0 );
  // And run all the tasks
  runAllTasks();
}

void ActionWithMatrix::updateAllNeighbourLists() {
  updateNeighbourList(); if( matrix_to_do_after ) matrix_to_do_after->updateAllNeighbourLists();
}

void ActionWithMatrix::performTask( const unsigned& task_index, MultiValue& myvals ) const {
  std::vector<unsigned> & indices( myvals.getIndices() ); 
  if( !doInnerLoop && actionInChain() ) {
      plumed_dbg_assert( myvals.inVectorCall() ); 
      runEndOfRowJobs( task_index, indices, myvals ); 
      return;
  }
  setupForTask( task_index, indices, myvals ); 

  // Now loop over the row of the matrix
  unsigned ntwo_atoms = myvals.getSplitIndex();
  for(unsigned i=0;i<ntwo_atoms;++i) {
      // This does everything in the stream that is done with single matrix elements
      runTask( getLabel(), task_index, indices[i], myvals );
      // Now clear only elements that are not accumulated over whole row
      clearMatrixElements( myvals );
  }
  // This updates the jobs that need to be completed when we get to the end of a row of the matrix 
  runEndOfRowJobs( task_index, indices, myvals );
}

void ActionWithMatrix::runTask( const std::string& controller, const unsigned& current, const unsigned colno, MultiValue& myvals ) const {
  double outval=0; myvals.setTaskIndex(current); myvals.setSecondTaskIndex( colno );
  if( isActive() ) performTask( controller, current, colno, myvals );
  bool hasval=false;
  for(int i=0; i<getNumberOfComponents(); ++i) { 
      if( fabs(myvals.get( getConstPntrToComponent(i)->getPositionInStream()) )>0 ) { hasval=true; break; }
  }

  if( hasval ) {
      for(int i=0; i<getNumberOfComponents(); ++i) {
        const Value* myval=getConstPntrToComponent(i);
        if( myval->getRank()!=2 || myval->hasDerivatives() || !myval->valueIsStored() ) continue;
        unsigned matindex = myval->getPositionInMatrixStash(), matbook_start = myval->getMatrixBookeepingStart(), col_stash_index = colno;
        if( colno>=myval->getShape()[0] ) col_stash_index = colno - myval->getShape()[0];
        if( myval->forcesWereAdded() ) {
          // unsigned sind = myval->getPositionInStream(), find = col_stash_index;
          // if( myval->getNumberOfColumns()<myval->getShape()[1] ) find=myvals.getNumberOfStashedMatrixElements(matindex);
          // double fforce = myval->getForce( myvals.getTaskIndex()*getNumberOfColumns() + find );
          // for(unsigned j=0; j<myvals.getNumberActive(sind); ++j) {
          //   unsigned kindex = myvals.getActiveIndex(sind,j); myvals.addMatrixForce( matindex, kindex, fforce*myvals.getDerivative(sind,kindex ) );
          // }
        } else myvals.stashMatrixElement( matindex, matbook_start+current*(1+getNumberOfColumns()), col_stash_index, myvals.get( myval->getPositionInStream() ) );
      } 
  }
  if( matrix_to_do_after ) matrix_to_do_after->runTask( controller, current, colno, myvals );
}

void ActionWithMatrix::gatherThreads( const unsigned& nt, const unsigned& bufsize, const std::vector<double>& omp_buffer, std::vector<double>& buffer, MultiValue& myvals ) {
  ActionWithVector::gatherThreads( nt, bufsize, omp_buffer, buffer, myvals );
  for(unsigned i=0; i<matrix_bookeeping.size(); ++i) matrix_bookeeping[i] += myvals.getMatrixBookeeping()[i];
}

void ActionWithMatrix::gatherProcesses( std::vector<double>& buffer ) {
  ActionWithVector::gatherProcesses( buffer );
  if( matrix_bookeeping.size()>0 && !runInSerial() ) comm.Sum( matrix_bookeeping ); 
  unsigned nval=0; transferNonZeroMatrixElementsToValues( nval, matrix_bookeeping );  
  if( matrix_to_do_after ) matrix_to_do_after->transferNonZeroMatrixElementsToValues( nval, matrix_bookeeping );
}

void ActionWithMatrix::transferNonZeroMatrixElementsToValues( unsigned& nval, const std::vector<unsigned>& matbook ) {
  for(int i=0; i<getNumberOfComponents(); ++i) {
      Value* myval=getPntrToComponent(i);
      if( myval->getRank()!=2 || myval->hasDerivatives() || !myval->valueIsStored() ) continue;
      unsigned nelements = myval->getShape()[0]*( 1 + getNumberOfColumns() ); 
      for(unsigned j=0; j<nelements; ++j) myval->setMatrixBookeepingElement( j, matbook[nval+j] );
      nval += nelements;
  }
}

void ActionWithMatrix::gatherStoredValue( const unsigned& valindex, const unsigned& code, const MultiValue& myvals,
                                          const unsigned& bufstart, std::vector<double>& buffer ) const {
  if( getConstPntrToComponent(valindex)->getRank()==1 ) { ActionWithVector::gatherStoredValue( valindex, code, myvals, bufstart, buffer ); return; }
  const Value* myval=getConstPntrToComponent(valindex);
  unsigned ncols = getNumberOfColumns(), matind = myval->getPositionInMatrixStash();
  unsigned matbook_start = myval->getMatrixBookeepingStart(), vindex = bufstart + code*getNumberOfColumns();
  const std::vector<unsigned> & matbook( myvals.getMatrixBookeeping() ); unsigned nelements = matbook[matbook_start+code*(1+ncols)];
  if( ncols>=myval->getShape()[1] ) {
       // In this case we store the full matrix
       for(unsigned j=0; j<nelements; ++j) {
         unsigned jind = matbook[matbook_start+code*(1+ncols)+1+j]; 
         plumed_dbg_massert( vindex+j<buffer.size(), "failing in " + getLabel() + " on value " + myval->getName() );
         buffer[vindex + jind] += myvals.getStashedMatrixElement( matind, jind );
       } 
  } else {
      // This is for storing sparse matrices when we can
      for(unsigned j=0; j<nelements; ++j) {
        unsigned jind = matbook[matbook_start+code*(1+ncols)+1+j]; 
        plumed_dbg_massert( vindex+j<buffer.size(), "failing in " + getLabel() + " on value " + myval->getName() );
        buffer[vindex + j] += myvals.getStashedMatrixElement( matind, jind );
      } 
  }
}

void ActionWithMatrix::clearMatrixElements( MultiValue& myvals ) const { 
  if( isActive() ) {
    for(int i=0; i<getNumberOfComponents(); ++i) {
      const Value* myval=getConstPntrToComponent(i);
      if( myval->getRank()==2 && !myval->hasDerivatives() ) myvals.clearDerivatives( myval->getPositionInStream() );
    }                            
  }                              
  if( matrix_to_do_after ) matrix_to_do_after->clearMatrixElements( myvals );
}

}
}
