/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2017 The plumed team
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
#include "MatrixProductBase.h"
#include "AdjacencyMatrixBase.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"

namespace PLMD {
namespace adjmat {

void MatrixProductBase::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys ); ActionAtomistic::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys ); ActionWithValue::registerKeywords( keys );
  keys.remove("NUMERICAL_DERIVATIVES"); keys.use("ARG");
}

MatrixProductBase::MatrixProductBase(const ActionOptions& ao):
  Action(ao),
  ActionAtomistic(ao),
  ActionWithArguments(ao),
  ActionWithValue(ao)
{
  if( getNumberOfArguments()!=2 ) error("should only have two arguments");
  if( getPntrToArgument(0)->getShape()[1]!=getPntrToArgument(1)->getShape()[0] ) error("number of columns in first matrix is not equal to number of columns in second");
  for(unsigned i=0; i<2; ++i) {
      if( getPntrToArgument(i)->getRank()!=2 || getPntrToArgument(i)->hasDerivatives() ) error("arguments should be matrices");
  }
  std::vector<Value*> args( getArguments() ); requestArguments( args, false ); 
  // Create a list of tasks for this action - n.b. each task calculates one row of the matrix
  std::vector<unsigned> shape(2); 
  shape[0]=getPntrToArgument(0)->getShape()[0];
  shape[1]=getPntrToArgument(1)->getShape()[1];

  for(unsigned j=0; j<shape[0]; ++j ) addTaskToList(j);
  // And create the matrix to hold the dot products
  addValue( shape ); 

  // Now do some stuff for time series
  bool timeseries=getPntrToArgument(0)->isTimeSeries();
  if( timeseries ) {
      for(unsigned i=1;i<getNumberOfArguments();++i) {
          if( !getPntrToArgument(i)->isTimeSeries() ) error("all arguments should either be time series or not time series");
      }
      getPntrToOutput(0)->makeTimeSeries();
  } else {
      for(unsigned i=1;i<getNumberOfArguments();++i) {
          if( getPntrToArgument(i)->isTimeSeries() ) error("all arguments should either be time series or not time series");
      }
  }
}

unsigned MatrixProductBase::getNumberOfDerivatives() const {
  unsigned numargs = ( getPntrToArgument(0)->getShape()[0] + getPntrToArgument(1)->getShape()[1])*getPntrToArgument(0)->getShape()[1];
  if( getNumberOfAtoms()>0 ) return 3*getNumberOfAtoms() + 9 + numargs;
  return numargs;
}

void MatrixProductBase::getTasksForParent( const std::string& parent, std::vector<std::string>& actionsThatSelectTasks, std::vector<unsigned>& tflags ) {
  if( tflags.size()!=getFullNumberOfTasks() ) return;
  // Check if parent has already been added
  bool found=false;
  for(unsigned i=0;i<actionsThatSelectTasks.size();++i) {
      if( actionsThatSelectTasks[i]==parent ) { found=true; break; }
  }
  if( found ) return;
  // Get the flags for this chain
  std::vector<unsigned> lflags( tflags.size(), 0 ); std::vector<unsigned> pTaskList, pIndexList;
  unsigned n_active = setTaskFlags( lflags, pTaskList, pIndexList );
  // Check if anything has been deactivated downstream
  if( n_active==tflags.size() ) return;
  // And retrieve non zero elements of contact matrix we are multiplying this by
  AdjacencyMatrixBase* ab = dynamic_cast<AdjacencyMatrixBase*>( getActionThatCalculates() );
  if( ab ) {
      // If tasks are deactivated in this child we can deactivate things in parent
      actionsThatSelectTasks.push_back( parent );
      // Get the atoms so that we can setup the neightbor lists
      ab->retrieveAtoms();
      // Now prepare the input matrix for the task loop
      ab->prepareForTasks( n_active, pTaskList );
      // Get the neighbours of each of the active atoms
      std::vector<unsigned> indices( getFullNumberOfTasks() );
      for(unsigned i=0;i<n_active;++i) {
          unsigned nneigh = ab->retrieveNeighbours( pTaskList[i], indices );
          for(unsigned j=0;j<nneigh;++j) tflags[indices[j]] = 1;
      }
      unsigned nacc = 0;
      for(unsigned i=0; i<tflags.size(); ++i) {
        if( tflags[i]>0 ) nacc++;
      }
  }
}

void MatrixProductBase::lockRequests() {
  ActionWithArguments::lockRequests();
  ActionAtomistic::lockRequests();
}

void MatrixProductBase::unlockRequests() {
  ActionWithArguments::unlockRequests();
  ActionAtomistic::unlockRequests();
}

void MatrixProductBase::calculateNumericalDerivatives( ActionWithValue* a ) { plumed_error(); }

void MatrixProductBase::calculate() {
  if( actionInChain() || skipCalculate() ) return;
  runAllTasks();
}

void MatrixProductBase::update() {
  if( skipUpdate() ) return;
  plumed_dbg_assert( !actionInChain() );
  if( getFullNumberOfTasks()>0 ) runAllTasks();
}

void MatrixProductBase::runFinalJobs() {
  if( skipUpdate() ) return;
  plumed_dbg_assert( !actionInChain() );
  resizeForFinalTasks(); runAllTasks();
}

void MatrixProductBase::updateCentralMatrixIndex( const unsigned& ind, MultiValue& myvals ) const {
  if( doNotCalculateDerivatives() ) return;

  unsigned nmat = getPntrToOutput(0)->getPositionInMatrixStash(), nmat_ind = myvals.getNumberOfMatrixIndices( nmat );
  std::vector<unsigned>& matrix_indices( myvals.getMatrixIndices( nmat ) ); unsigned invals=getFullNumberOfTasks(); 

  unsigned nargs = getPntrToArgument(0)->getShape()[1];
  for(unsigned i=0; i<nargs; ++i) { matrix_indices[nmat_ind] = nargs*ind + i; nmat_ind++; }
  if( getNumberOfAtoms()>0 ) {
    unsigned numargs = ( getPntrToArgument(0)->getShape()[0] + getPntrToArgument(1)->getShape()[1])*getPntrToArgument(0)->getShape()[1];
    matrix_indices[nmat_ind+0]=numargs + 3*ind+0;
    matrix_indices[nmat_ind+1]=numargs + 3*ind+1;
    matrix_indices[nmat_ind+2]=numargs + 3*ind+2;
    nmat_ind+=3; unsigned virbase = numargs + 3*getNumberOfAtoms();
    for(unsigned i=0; i<9; ++i) matrix_indices[nmat_ind+i]=virbase+i;
    nmat_ind+=9;
  }
  myvals.setNumberOfMatrixIndices( nmat, nmat_ind );
}

unsigned MatrixProductBase::getNumberOfColumns() const { 
  plumed_massert( !actionInChain(), "I am not sure how to do this so I am not allowing it GAT");  
  return getPntrToOutput(0)->getShape()[1];
}

void MatrixProductBase::performTask( const unsigned& current, MultiValue& myvals ) const {
  if( actionInChain() ) {
    if( myvals.inVectorCall() ) updateCentralMatrixIndex( myvals.getTaskIndex(), myvals );
    return ;
  }

  // Now loop over all atoms in coordination sphere
  unsigned start_n = getFullNumberOfTasks();
  myvals.setNumberOfIndicesInFirstBlock( start_n );
  unsigned size_v = getPntrToOutput(0)->getShape()[1];
  for(unsigned i=0; i<size_v; ++i) {
    // This does everything in the stream that is done with single matrix elements
    runTask( getLabel(), myvals.getTaskIndex(), current, start_n + i, myvals );
    // Now clear only elements that are not accumulated over whole row
    clearMatrixElements( myvals );
  }
  // Update the matrix index for the central atom
  updateCentralMatrixIndex( myvals.getTaskIndex(), myvals );
}

bool MatrixProductBase::performTask( const std::string& controller, const unsigned& index1, const unsigned& index2, MultiValue& myvals ) const {
  unsigned ind2 = index2; if( index2>=getFullNumberOfTasks() ) ind2 = index2 - getFullNumberOfTasks();
  unsigned sss = getPntrToArgument(1)->getShape()[1], nargs = getPntrToArgument(0)->getShape()[1];
  std::vector<double> args1(nargs), args2(nargs), der1(nargs), der2(nargs);
  for(unsigned i=0; i<nargs; ++i) {
    args1[i] = getPntrToArgument(0)->get( index1*nargs + i );
    args2[i] = getPntrToArgument(1)->get( i*sss + ind2 ); 
  }
  double val = computeVectorProduct( index1, ind2, args1, args2, der1, der2, myvals );
  unsigned ostrn = getPntrToOutput(0)->getPositionInStream();
  myvals.setValue( ostrn, val );
  // Return after calculation of value if we do not need derivatives
  if( doNotCalculateDerivatives() ) return true;

  unsigned nmat = getPntrToOutput(0)->getPositionInMatrixStash();
  std::vector<unsigned>& matrix_indices( myvals.getMatrixIndices( nmat ) );
  plumed_dbg_assert( matrix_indices.size()>=getNumberOfDerivatives() );
  unsigned nmat_ind = myvals.getNumberOfMatrixIndices( nmat );
  unsigned jind_start = getPntrToArgument(0)->getShape()[0]*getPntrToArgument(0)->getShape()[1];
  for(unsigned i=0; i<nargs; ++i) {
    plumed_dbg_assert( nargs*index1 + i<myvals.getNumberOfDerivatives() );
    myvals.addDerivative( ostrn, nargs*index1 + i, der1[i] );
    myvals.updateIndex( ostrn, nargs*index1 + i );
    plumed_dbg_assert( jind_start + i*sss + ind2<myvals.getNumberOfDerivatives() );
    myvals.addDerivative( ostrn, jind_start + i*sss + ind2, der2[i] );
    myvals.updateIndex( ostrn, jind_start + i*sss + ind2 );
    matrix_indices[nmat_ind] = jind_start + i*sss + ind2;
    nmat_ind++;
  }
  if( getNumberOfAtoms()>0 ) {
    unsigned numargs = ( getPntrToArgument(0)->getShape()[0] + getPntrToArgument(1)->getShape()[1])*getPntrToArgument(0)->getShape()[1];
    matrix_indices[nmat_ind+0]=numargs + 3*ind2+0;
    matrix_indices[nmat_ind+1]=numargs + 3*ind2+1;
    matrix_indices[nmat_ind+2]=numargs + 3*ind2+2;
    nmat_ind+=3;
  }
  myvals.setNumberOfMatrixIndices( nmat, nmat_ind );
  return true;
}

void MatrixProductBase::apply() {
  if( doNotCalculateDerivatives() ) return;
  std::fill(forcesToApply.begin(),forcesToApply.end(),0); unsigned mm=0;
  if( getForcesFromValues( forcesToApply ) ) {
    setForcesOnAtoms( forcesToApply, mm );
    setForcesOnArguments( 0, forcesToApply, mm );
  }
}

}
}
