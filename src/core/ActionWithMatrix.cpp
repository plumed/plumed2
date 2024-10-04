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
#include "tools/Communicator.h"

namespace PLMD {

void ActionWithMatrix::registerKeywords( Keywords& keys ) {
  ActionWithVector::registerKeywords( keys );
  keys.use("ARG");
}

ActionWithMatrix::ActionWithMatrix(const ActionOptions&ao):
  Action(ao),
  ActionWithVector(ao),
  next_action_in_chain(NULL),
  matrix_to_do_before(NULL),
  matrix_to_do_after(NULL),
  clearOnEachCycle(true) {
}

ActionWithMatrix::~ActionWithMatrix() {
  if( matrix_to_do_before ) {
    matrix_to_do_before->matrix_to_do_after=NULL;
    matrix_to_do_before->next_action_in_chain=NULL;
  }
}

void ActionWithMatrix::getAllActionLabelsInMatrixChain( std::vector<std::string>& mylabels ) const {
  bool found=false;
  for(unsigned i=0; i<mylabels.size(); ++i) {
    if( getLabel()==mylabels[i] ) {
      found=true;
    }
  }
  if( !found ) {
    mylabels.push_back( getLabel() );
  }
  if( matrix_to_do_after ) {
    matrix_to_do_after->getAllActionLabelsInMatrixChain( mylabels );
  }
}

void ActionWithMatrix::setupStreamedComponents( const std::string& headstr, unsigned& nquants, unsigned& nmat, unsigned& maxcol, unsigned& nbookeeping ) {
  ActionWithVector::setupStreamedComponents( headstr, nquants, nmat, maxcol, nbookeeping );

  for(int i=0; i<getNumberOfComponents(); ++i) {
    Value* myval=getPntrToComponent(i);
    if( myval->getRank()!=2 || myval->hasDerivatives() ) {
      continue;
    }
    myval->setPositionInMatrixStash(nmat);
    nmat++;
    if( !myval->valueIsStored() ) {
      continue;
    }
    if( myval->getShape()[1]>maxcol ) {
      maxcol=myval->getShape()[1];
    }
    myval->setMatrixBookeepingStart(nbookeeping);
    nbookeeping += myval->getShape()[0]*( 1 + myval->getNumberOfColumns() );
  }
  // Turn off clearning of derivatives after each matrix run if there are no matrices in the output of this action
  clearOnEachCycle = false;
  for(int i=0; i<getNumberOfComponents(); ++i) {
    const Value* myval=getConstPntrToComponent(i);
    if( myval->getRank()==2 && !myval->hasDerivatives() ) {
      clearOnEachCycle = true;
      break;
    }
  }
  // Turn off clearing of derivatives if we have only the values of adjacency matrices
  if( doNotCalculateDerivatives() && isAdjacencyMatrix() ) {
    clearOnEachCycle = false;
  }
}

void ActionWithMatrix::finishChainBuild( ActionWithVector* act ) {
  ActionWithMatrix* am=dynamic_cast<ActionWithMatrix*>(act);
  if( !am || act==this ) {
    return;
  }
  // Build the list that contains everything we are going to loop over in getTotalMatrixBookeepgin and updateAllNeighbourLists
  if( next_action_in_chain ) {
    next_action_in_chain->finishChainBuild( act );
  } else {
    next_action_in_chain=am;
    // Build the list of things we are going to loop over in runTask
    if( am->isAdjacencyMatrix() || act->getName()=="VSTACK" ) {
      return ;
    }
    plumed_massert( !matrix_to_do_after, "cannot add " + act->getLabel() + " in " + getLabel() + " as have already added " + matrix_to_do_after->getLabel() );
    matrix_to_do_after=am;
    am->matrix_to_do_before=this;
  }
}

const ActionWithMatrix* ActionWithMatrix::getFirstMatrixInChain() const {
  if( !actionInChain() ) {
    return this;
  }
  return matrix_to_do_before->getFirstMatrixInChain();
}

void ActionWithMatrix::getTotalMatrixBookeeping( unsigned& nbookeeping ) {
  for(int i=0; i<getNumberOfComponents(); ++i) {
    Value* myval=getPntrToComponent(i);
    if( myval->getRank()!=2 || myval->hasDerivatives() || !myval->valueIsStored() ) {
      continue;
    }
    myval->reshapeMatrixStore( getNumberOfColumns() );
    nbookeeping += myval->getShape()[0]*( 1 + myval->getNumberOfColumns() );
  }
  if( next_action_in_chain ) {
    next_action_in_chain->getTotalMatrixBookeeping( nbookeeping );
  }
}

void ActionWithMatrix::calculate() {
  if( actionInChain() ) {
    return ;
  }
  // Update all the neighbour lists
  updateAllNeighbourLists();
  // Setup the matrix indices
  unsigned nbookeeping=0;
  getTotalMatrixBookeeping( nbookeeping );
  if( matrix_bookeeping.size()!=nbookeeping ) {
    matrix_bookeeping.resize( nbookeeping );
  }
  std::fill( matrix_bookeeping.begin(), matrix_bookeeping.end(), 0 );
  // And run all the tasks
  runAllTasks();
}

void ActionWithMatrix::updateAllNeighbourLists() {
  updateNeighbourList();
  if( next_action_in_chain ) {
    next_action_in_chain->updateAllNeighbourLists();
  }
}

void ActionWithMatrix::performTask( const unsigned& task_index, MultiValue& myvals ) const {
  std::vector<unsigned> & indices( myvals.getIndices() );
  if( matrix_to_do_before ) {
    plumed_dbg_assert( myvals.inVectorCall() );
    runEndOfRowJobs( task_index, indices, myvals );
    return;
  }
  setupForTask( task_index, indices, myvals );

  // Now loop over the row of the matrix
  unsigned ntwo_atoms = myvals.getSplitIndex();
  for(unsigned i=1; i<ntwo_atoms; ++i) {
    // This does everything in the stream that is done with single matrix elements
    runTask( getLabel(), task_index, indices[i], myvals );
    // Now clear only elements that are not accumulated over whole row
    clearMatrixElements( myvals );
  }
  // This updates the jobs that need to be completed when we get to the end of a row of the matrix
  runEndOfRowJobs( task_index, indices, myvals );
}

void ActionWithMatrix::runTask( const std::string& controller, const unsigned& current, const unsigned colno, MultiValue& myvals ) const {
  double outval=0;
  myvals.setTaskIndex(current);
  myvals.setSecondTaskIndex( colno );
  if( isActive() ) {
    performTask( controller, current, colno, myvals );
  }
  bool hasval = !isAdjacencyMatrix();
  for(int i=0; i<getNumberOfComponents(); ++i) {
    if( fabs(myvals.get( getConstPntrToComponent(i)->getPositionInStream()) )>0 ) {
      hasval=true;
      break;
    }
  }

  if( hasval ) {
    for(int i=0; i<getNumberOfComponents(); ++i) {
      const Value* myval=getConstPntrToComponent(i);
      if( myval->getRank()!=2 || myval->hasDerivatives() || !myval->valueIsStored() ) {
        continue;
      }
      unsigned matindex = myval->getPositionInMatrixStash(), matbook_start = myval->getMatrixBookeepingStart(), col_stash_index = colno;
      if( colno>=myval->getShape()[0] ) {
        col_stash_index = colno - myval->getShape()[0];
      }
      unsigned rowstart = matbook_start+current*(1+myval->getNumberOfColumns());
      if( myval->forcesWereAdded() ) {
        unsigned sind = myval->getPositionInStream(), find = myvals.getMatrixBookeeping()[rowstart];
        double fforce = myval->getForce( myvals.getTaskIndex()*myval->getNumberOfColumns() + find );
        if( getNumberOfColumns()>=myval->getShape()[1] ) {
          fforce = myval->getForce( myvals.getTaskIndex()*myval->getShape()[1] + col_stash_index );
        }
        for(unsigned j=0; j<myvals.getNumberActive(sind); ++j) {
          unsigned kindex = myvals.getActiveIndex(sind,j);
          myvals.addMatrixForce( matindex, kindex, fforce*myvals.getDerivative(sind,kindex ) );
        }
      }
      double finalval = myvals.get( myval->getPositionInStream() );
      if( fabs(finalval)>0 ) {
        myvals.stashMatrixElement( matindex, rowstart, col_stash_index, finalval );
      }
    }
  }
  if( matrix_to_do_after ) {
    matrix_to_do_after->runTask( controller, current, colno, myvals );
  }
}

void ActionWithMatrix::gatherThreads( const unsigned& nt, const unsigned& bufsize, const std::vector<double>& omp_buffer, std::vector<double>& buffer, MultiValue& myvals ) {
  ActionWithVector::gatherThreads( nt, bufsize, omp_buffer, buffer, myvals );
  for(unsigned i=0; i<matrix_bookeeping.size(); ++i) {
    matrix_bookeeping[i] += myvals.getMatrixBookeeping()[i];
  }
}

void ActionWithMatrix::gatherProcesses( std::vector<double>& buffer ) {
  ActionWithVector::gatherProcesses( buffer );
  if( matrix_bookeeping.size()>0 && !runInSerial() ) {
    comm.Sum( matrix_bookeeping );
  }
  unsigned nval=0;
  transferNonZeroMatrixElementsToValues( nval, matrix_bookeeping );
}

void ActionWithMatrix::transferNonZeroMatrixElementsToValues( unsigned& nval, const std::vector<unsigned>& matbook ) {
  for(int i=0; i<getNumberOfComponents(); ++i) {
    Value* myval=getPntrToComponent(i);
    if( myval->getRank()!=2 || myval->hasDerivatives() || !myval->valueIsStored() || getNumberOfColumns()>=myval->getShape()[1] ) {
      continue;
    }
    unsigned nelements = myval->getShape()[0]*( 1 + myval->getNumberOfColumns() );
    for(unsigned j=0; j<nelements; ++j) {
      myval->setMatrixBookeepingElement( j, matbook[nval+j] );
    }
    nval += nelements;
  }
  if( next_action_in_chain ) {
    next_action_in_chain->transferNonZeroMatrixElementsToValues( nval, matbook );
  }
}

void ActionWithMatrix::gatherStoredValue( const unsigned& valindex, const unsigned& code, const MultiValue& myvals,
    const unsigned& bufstart, std::vector<double>& buffer ) const {
  if( getConstPntrToComponent(valindex)->getRank()==1 ) {
    ActionWithVector::gatherStoredValue( valindex, code, myvals, bufstart, buffer );
    return;
  }
  const Value* myval=getConstPntrToComponent(valindex);
  unsigned ncols = myval->getNumberOfColumns(), matind = myval->getPositionInMatrixStash();
  unsigned matbook_start = myval->getMatrixBookeepingStart(), vindex = bufstart + code*myval->getNumberOfColumns();
  const std::vector<unsigned> & matbook( myvals.getMatrixBookeeping() );
  unsigned nelements = matbook[matbook_start+code*(1+ncols)];
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

bool ActionWithMatrix::checkForTaskForce( const unsigned& itask, const Value* myval ) const {
  if( myval->getRank()<2 ) {
    return ActionWithVector::checkForTaskForce( itask, myval );
  }
  unsigned nelements = myval->getRowLength(itask), startr = itask*myval->getNumberOfColumns();
  for(unsigned j=0; j<nelements; ++j ) {
    if( fabs( myval->getForce( startr + j ) )>epsilon ) {
      return true;
    }
  }
  return false;
}

void ActionWithMatrix::gatherForcesOnStoredValue( const Value* myval, const unsigned& itask, const MultiValue& myvals, std::vector<double>& forces ) const {
  if( myval->getRank()==1 ) {
    ActionWithVector::gatherForcesOnStoredValue( myval, itask, myvals, forces );
    return;
  }
  unsigned matind = myval->getPositionInMatrixStash();
  for(unsigned j=0; j<forces.size(); ++j) {
    forces[j] += myvals.getStashedMatrixForce( matind, j );
  }
}

void ActionWithMatrix::clearMatrixElements( MultiValue& myvals ) const {
  if( isActive() && clearOnEachCycle ) {
    for(int i=0; i<getNumberOfComponents(); ++i) {
      const Value* myval=getConstPntrToComponent(i);
      if( myval->getRank()==2 && !myval->hasDerivatives() ) {
        myvals.clearDerivatives( myval->getPositionInStream() );
      }
    }
  }
  if( matrix_to_do_after ) {
    matrix_to_do_after->clearMatrixElements( myvals );
  }
}

}
