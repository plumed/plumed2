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
  ActionWithVector::registerKeywords( keys ); keys.use("ARG");
}

ActionWithMatrix::ActionWithMatrix(const ActionOptions&ao):
  Action(ao),
  ActionWithVector(ao)
{
}

void ActionWithMatrix::calculate() {
  // Update all the neighbour lists
  updateNeighbourList();
  // Setup the matrix indices
  for(int i=0; i<getNumberOfComponents(); ++i) {
    Value* myval=getPntrToComponent(i);
    if( myval->getRank()!=2 || myval->hasDerivatives() ) continue;
    myval->reshapeMatrixStore( getNumberOfColumns() );
  }
  // And run all the tasks
  runAllTasks();
}

void ActionWithMatrix::performTask( const unsigned& task_index, MultiValue& myvals ) const {
  std::vector<unsigned> & indices( myvals.getIndices() );
  setupForTask( task_index, indices, myvals );

  // Reset the bookeeping elements for storage
  for(unsigned i=0; i<getNumberOfComponents(); ++i) {
    Value* myval = const_cast<Value*>( getConstPntrToComponent(i) ); unsigned ncols = myval->getNumberOfColumns();
    if( myval->getRank()!=2 || myval->hasDerivatives() || ncols>=myval->getShape()[1] ) continue;
    myval->matrix_bookeeping[task_index*(1+ncols)]=0;
  }

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
  double outval=0; myvals.setTaskIndex(current); myvals.setSecondTaskIndex( colno );
  performTask( controller, current, colno, myvals );
  bool hasval = !isAdjacencyMatrix();
  for(int i=0; i<getNumberOfComponents(); ++i) {
    if( std::fabs(myvals.get(i))>epsilon ) { hasval=true; break; }
  }

  if( hasval ) {
    double checkval = myvals.get( 0 );
    for(int i=0; i<getNumberOfComponents(); ++i) {
      const Value* myval=getConstPntrToComponent(i); unsigned ncols = myval->getNumberOfColumns();
      if( myval->getRank()!=2 || myval->hasDerivatives() ) continue;
      unsigned col_stash_index = colno;
      if( colno>=myval->getShape()[0] ) col_stash_index = colno - myval->getShape()[0];
      if( myval->forcesWereAdded() ) {
        double fforce;
        if( ncols<myval->getShape()[1] ) fforce = myval->getForce( myvals.getTaskIndex()*ncols + myval->matrix_bookeeping[current*(1+ncols)] );
        else fforce = myval->getForce( myvals.getTaskIndex()*myval->getShape()[1] + col_stash_index );
        for(unsigned j=0; j<myvals.getNumberActive(i); ++j) {
          unsigned kindex = myvals.getActiveIndex(i,j); myvals.addMatrixForce( kindex, fforce*myvals.getDerivative(i,kindex ) );
        }
      }
      double finalval = myvals.get( i ); if( !isAdjacencyMatrix() ) checkval=finalval;
      if( std::fabs(checkval)>epsilon || (!isAdjacencyMatrix() && getNumberOfMasks()>0) ) {
        Value* myv = const_cast<Value*>( myval );
        if( ncols<myval->getShape()[1] ) {
          myv->set( current*ncols + myval->matrix_bookeeping[current*(1+ncols)], finalval );
          myv->matrix_bookeeping[current*(1+ncols)]++; myv->matrix_bookeeping[current*(1+ncols)+myval->matrix_bookeeping[current*(1+ncols)]] = col_stash_index;
        } else myv->set( current*myval->getShape()[1] + col_stash_index, finalval );
      }
    }
  }
}

bool ActionWithMatrix::checkForTaskForce( const unsigned& itask, const Value* myval ) const {
  if( myval->getRank()<2 ) return ActionWithVector::checkForTaskForce( itask, myval );
  unsigned nelements = myval->getRowLength(itask), startr = itask*myval->getNumberOfColumns();
  for(unsigned j=0; j<nelements; ++j ) {
    if( fabs( myval->getForce( startr + j ) )>epsilon ) return true;
  }
  return false;
}

void ActionWithMatrix::gatherForces( const unsigned& itask, const MultiValue& myvals, std::vector<double>& forces ) const {
  if( checkComponentsForForce() ) {
    const std::vector<unsigned>& mat_indices( myvals.getMatrixRowDerivativeIndices() );
    for(unsigned i=0; i<myvals.getNumberOfMatrixRowDerivatives(); ++i) { unsigned kind = mat_indices[i]; forces[kind] += myvals.getStashedMatrixForce( kind ); }
  }
}

void ActionWithMatrix::clearMatrixElements( MultiValue& myvals ) const {
  for(int i=0; i<getNumberOfComponents(); ++i) {
    const Value* myval=getConstPntrToComponent(i);
    if( myval->getRank()==2 && !myval->hasDerivatives() ) myvals.clearDerivatives( i );
  }
}

}
