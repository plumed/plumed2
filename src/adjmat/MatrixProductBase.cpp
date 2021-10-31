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
#include "core/ActionSetup.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"

namespace PLMD {
namespace adjmat {

void MatrixProductBase::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys ); ActionAtomistic::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys ); ActionWithValue::registerKeywords( keys );
  keys.remove("NUMERICAL_DERIVATIVES"); keys.use("ARG");
  keys.addFlag("ELEMENTS_ON_DIAGONAL_ARE_ZERO",false,"set all diagonal elements equal to zero");
  ActionWithValue::useCustomisableComponents(keys);
}

MatrixProductBase::MatrixProductBase(const ActionOptions& ao):
  Action(ao),
  ActionAtomistic(ao),
  ActionWithArguments(ao),
  ActionWithValue(ao),
  skip_ieqj(false),
  diagonal(false),
  doInnerLoop(false)
{
  if( keywords.exists("DIAGONAL_ELEMENTS_ONLY") ) parseFlag("DIAGONAL_ELEMENTS_ONLY",diagonal);
}

void MatrixProductBase::readMatricesToMultiply( const bool& periodic, const std::string& min, const std::string& max ) {
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
  if( getPntrToArgument(i)->getRank()==0 || getPntrToArgument(i)->hasDerivatives() ) error("arguments should be matrices or vectors");
  }
  std::vector<unsigned> shape( getMatrixShapeForFinalTasks() ); unsigned noutput = getNumberOfArguments() / 2;
  bool zero_diag=false; if( getPntrToArgument(0)->getRank()==1 && getPntrToArgument(noutput)->getRank()==1 ) parseFlag("ELEMENTS_ON_DIAGONAL_ARE_ZERO",zero_diag);
  if( zero_diag && shape[0]!=shape[1] ) error("cannot set diagonal elements of matrix to zero if matrix is not square");
  else if( zero_diag ) skip_ieqj=true;
  if( skip_ieqj ) log.printf("  ignoring diagonal elements of matrix \n");
  // Create a list of tasks for this action - n.b. each task calculates one row of the matrix
  if( shape.size()==0 ) addTaskToList(0);
  else { for(unsigned j=0; j<shape[0]; ++j ) addTaskToList(j); }
  // And create the matrix to hold the dot products
  if( noutput>1 ) {
      for(unsigned i=0; i<noutput; ++i) {
          std::string name, fullname = getPntrToArgument(i)->getName();
          std::size_t dot = fullname.find_first_of(".");
          std::size_t und = fullname.find_first_of("_");
          if( fullname.find(".")!=std::string::npos ) name = fullname.substr(dot+1);
          else if( fullname.find("_")!=std::string::npos ) name = fullname.substr(und+1);
          else name = fullname;
          addComponent( name, shape ); 
          if( periodic ) componentIsPeriodic( name, min, max ); else componentIsNotPeriodic( name );
      }
  } else {
      if( shape.size()==0 ) addValueWithDerivatives( shape ); else addValue( shape ); 
      if( periodic ) setPeriodic( min, max ); else setNotPeriodic();  
  }

  for(unsigned nv=0; nv<noutput; ++nv ) {
      // This sets up matrix times vector calculations that are done without storing the input matrix
      if( !diagonal && !getPntrToArgument(nv)->dataAlwaysStored() && !getPntrToArgument(nv)->isHistoryDependent() && 
          getPntrToArgument(nv)->getRank()==2 && getPntrToArgument(noutput+nv)->getRank()==1 && 
          (getPntrToArgument(noutput+nv)->getPntrToAction())->getName()!="TRANSPOSE" ) {
          // Chain off the action that computes the matrix
          std::vector<std::string> alabels(1); alabels[0]=(getPntrToArgument(nv)->getPntrToAction())->getLabel();
          (getPntrToArgument(nv)->getPntrToAction())->addActionToChain( alabels, this );
      } else getPntrToArgument(nv)->buildDataStore( getLabel() );
      // Store the vector or second matrix
      getPntrToArgument(noutput+nv)->buildDataStore( getLabel() );
  }
  // Rerequest arguments 
  std::vector<Value*> args( getArguments() ); requestArguments( args, true );
}

bool MatrixProductBase::canBeAfterInChain( ActionWithValue* av ) {
  // Input argument is the one that comes first in the order
  // We can be after this if it is not a matrix product, or if we or if we are outputting a vector
  MatrixProductBase* mp = dynamic_cast<MatrixProductBase*>(av);
  if( !mp ) return true;
  // Now check if there is a mismatch between the number of tasks in the two actions
  if( mp->getNumberOfInnerTasks()!=getNumberOfInnerTasks() && getPntrToOutput(0)->getRank()!=1 ) doInnerLoop = true;
  // If it is a vector or we are a vector everything is fine
  if( getPntrToOutput(0)->getRank()<2 || (av->copyOutput(0))->getRank()<2 ) return true;
  // We can be after it if is an AdjacencyMatrix
  AdjacencyMatrixBase* ab=dynamic_cast<AdjacencyMatrixBase*>(av);
  if( ab ) return true;
  // We can't be after it if we are an Adjacency matrix but it isn't
  const AdjacencyMatrixBase* ab2=dynamic_cast<const AdjacencyMatrixBase*>(this);
  if( ab2 ) return false;
  // The skip_ieqj thing needs putting back in Gareth
  // if( ab2 || mp->skip_ieqj!=skip_ieqj ) return false;
  return true;
}

unsigned MatrixProductBase::getNumberOfDerivatives() const {
  unsigned noutput = getNumberOfComponents();
  if( actionInChain() && getPntrToOutput(0)->getRank()<2 ) return (getPntrToArgument(0)->getPntrToAction())->getNumberOfDerivatives() + getPntrToArgument(noutput)->getSize();
  unsigned numargs = 0; if( getNumberOfArguments()>0 ) numargs = getPntrToArgument(0)->getSize() + getPntrToArgument(noutput)->getSize(); 
  if( getNumberOfAtoms()>0 ) return 3*getNumberOfAtoms() + 9 + numargs;
  return numargs;
}

bool MatrixProductBase::mustBeTreatedAsDistinctArguments() {
  AdjacencyMatrixBase* ab=dynamic_cast<AdjacencyMatrixBase*>(this);
  if( ab ) return ActionWithArguments::mustBeTreatedAsDistinctArguments();
  return true;
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
  } else {
      if( actionsLabelsInChain.size()==0 ) getAllActionLabelsInChain( actionsLabelsInChain );
      bool ignore = checkUsedOutsideOfChain( actionsLabelsInChain, parent, actionsThatSelectTasks, tflags );
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
  if( getPntrToOutput(0)->getRank()==2 ) done_over_stream=false; 
  runAllTasks();
}

void MatrixProductBase::update() {
  if( actionInChain() || skipUpdate() ) return;
  plumed_dbg_assert( !actionInChain() );
  if( getFullNumberOfTasks()>0 ) runAllTasks();
}

void MatrixProductBase::runFinalJobs() {
  if( actionInChain() || skipUpdate() ) return;
  plumed_dbg_assert( !actionInChain() );
  resizeForFinalTasks();
  runAllTasks();
}

unsigned MatrixProductBase::getNumberOfFinalTasks() {
  return getMatrixShapeForFinalTasks()[0];
}

std::vector<unsigned> MatrixProductBase::getMatrixShapeForFinalTasks() {
  std::vector<unsigned> shape(2); unsigned noutput=getNumberOfArguments() / 2;
  if( getPntrToArgument(0)->getRank()==1 && getPntrToArgument(noutput)->getRank()==1 ) {
      shape[0]=getPntrToArgument(0)->getShape()[0]; shape[1]=getPntrToArgument(noutput)->getShape()[0];
      for(unsigned i=1;i<noutput;++i) {
          if( getPntrToArgument(i)->getRank()!=1 || getPntrToArgument(noutput+i)->getRank()!=1 ) error("all input arguments should have rank one");
          if( shape[0]!=getPntrToArgument(i)->getShape()[0] ) error("shapes of all output matrices should match");
          if( shape[1]!=getPntrToArgument(noutput+i)->getShape()[1] ) error("shapes of all output matrices should match");
      }
  } else if( getPntrToArgument(0)->getRank()==2 && getPntrToArgument(noutput)->getRank()==2 ) {
      if( getPntrToArgument(0)->getShape()[1]!=getPntrToArgument(noutput)->getShape()[0] ) error("number of columns in first matrix is not equal to number of columns in second");
      shape[0]=getPntrToArgument(0)->getShape()[0]; shape[1]=getPntrToArgument(noutput)->getShape()[1];
      for(unsigned i=1;i<noutput;++i) {
          if( getPntrToArgument(i)->getRank()!=2 || getPntrToArgument(noutput+i)->getRank()!=2 ) error("all input arguments should have rank two");
          if( getPntrToArgument(i)->getShape()[1]!=getPntrToArgument(noutput+i)->getShape()[0] ) error("number of columns in first matrix is not equal to number of columns in second");
          if( shape[0]!=getPntrToArgument(i)->getShape()[0] ) error("shapes of all output matrices should match");
          if( shape[1]!=getPntrToArgument(noutput+i)->getShape()[1] ) error("shapes of all output matrices should match"); 
      }
      if( noutput==1 && !diagonal ) {
          // Check if we are multiplying a matrix by its transpose (if we are doing this we know the diagonal elements are all 1 or something similarly boring)
          if( (getPntrToArgument(0)->getPntrToAction())->getName()=="TRANSPOSE" ) {
               ActionWithArguments* aa = dynamic_cast<ActionWithArguments*>( getPntrToArgument(0)->getPntrToAction() );
               if( (aa->getPntrToArgument(0))->getName()==getPntrToArgument(1)->getName() && (getPntrToArgument(1)->getPntrToAction())->getName().find("STACK")!=std::string::npos ) skip_ieqj=true;
          } else if( (getPntrToArgument(1)->getPntrToAction())->getName()=="TRANSPOSE" ) {
               ActionWithArguments* aa = dynamic_cast<ActionWithArguments*>( getPntrToArgument(1)->getPntrToAction() );
               if( (aa->getPntrToArgument(0))->getName()==getPntrToArgument(0)->getName() && (getPntrToArgument(0)->getPntrToAction())->getName().find("STACK")!=std::string::npos ) skip_ieqj=true;
          }
      } else if( diagonal ) { 
          if( shape[0]!=shape[1] ) error("can only use DIAGONAL_ELEMENTS_ONLY flag if output matrix is square");
          unsigned nn=shape[0]; shape.resize(1); shape[0]=nn; 
      }
  } else if( getPntrToArgument(0)->getRank()==2 && getPntrToArgument(noutput)->getRank()==1 ) {
      if( getPntrToArgument(0)->getShape()[1]!=getPntrToArgument(noutput)->getShape()[0] ) error("number of columns in first matrix is not equal to number of elements in vector");
      for(unsigned i=1;i<noutput;++i) {
          if( getPntrToArgument(i)->getRank()!=2 || getPntrToArgument(noutput+i)->getRank()!=1 ) error("all first input arguments should have rank two and all second should have rank one");
          if( getPntrToArgument(i)->getShape()[1]!=getPntrToArgument(noutput+i)->getShape()[0] ) error("number of columns in first matrix is not equal to number of elements in vector");
      }
      if( diagonal && getPntrToArgument(0)->getShape()[0]==1 ) {
         shape.resize(0);
         for(unsigned i=1;i<noutput;++i) {
             if( 1!=getPntrToArgument(i)->getShape()[0] ) error("shapes of all output matrices should match");
         }
      } else {
         shape.resize(1); shape[0] = getPntrToArgument(0)->getShape()[0];
         for(unsigned i=1;i<noutput;++i) {
             if( shape[0]!=getPntrToArgument(i)->getShape()[0] ) error("shapes of all output matrices should match");
         }
     }
  } else error("cannot do product of row vector with matrix transpose matrix in ARG2 and do use DOT with matrix as ARG1 vector as ARG2");
  return shape;
}

void MatrixProductBase::updateCentralMatrixIndex( const unsigned& ind, const std::vector<unsigned>& indices, MultiValue& myvals ) const {
  if( getPntrToOutput(0)->getRank()<2 ) {
      if( !actionInChain() ) return ;
      for(unsigned nv=0; nv<getNumberOfComponents(); ++nv ) {
          unsigned istrn = getPntrToArgument(nv)->getPositionInMatrixStash();
          std::vector<unsigned>& mat_indices( myvals.getMatrixIndices( istrn ) );
          for(unsigned i=0; i<myvals.getNumberOfMatrixIndices(istrn); ++i) {
              unsigned ostrn = getPntrToOutput(nv)->getPositionInStream();
              myvals.updateIndex( ostrn, mat_indices[i] );
          }
      }
  } else {
      for(unsigned k=0; k<getNumberOfComponents(); ++k ) {
          unsigned nmat = getPntrToOutput(k)->getPositionInMatrixStash(), nmat_ind = myvals.getNumberOfMatrixIndices( nmat );
          std::vector<unsigned>& matrix_indices( myvals.getMatrixIndices( nmat ) ); unsigned invals=getFullNumberOfTasks(); 

          unsigned numargs = 0;
          if( getNumberOfArguments()>0 ) {
              unsigned nargs = 1; if( getPntrToArgument(0)->getRank()==2 ) nargs = getPntrToArgument(0)->getShape()[1];
              for(unsigned i=0; i<nargs; ++i) { matrix_indices[nmat_ind] = nargs*ind + i; nmat_ind++; }
              numargs = getPntrToArgument(0)->getSize() + getPntrToArgument(1)->getSize();
          }
          if( getNumberOfAtoms()>0 ) {
            matrix_indices[nmat_ind+0]=numargs + 3*ind+0;
            matrix_indices[nmat_ind+1]=numargs + 3*ind+1;
            matrix_indices[nmat_ind+2]=numargs + 3*ind+2;
            nmat_ind+=3; 
            for(unsigned i=myvals.getSplitIndex(); i<myvals.getNumberOfIndices(); ++i) { 
                matrix_indices[nmat_ind+0]=numargs + 3*indices[i]+0;
                matrix_indices[nmat_ind+1]=numargs + 3*indices[i]+1;
                matrix_indices[nmat_ind+2]=numargs + 3*indices[i]+2;
                nmat_ind+=3;
            }
            unsigned virbase = numargs + 3*getNumberOfAtoms();
            for(unsigned i=0; i<9; ++i) matrix_indices[nmat_ind+i]=virbase+i;
            nmat_ind+=9; 
          }
          myvals.setNumberOfMatrixIndices( nmat, nmat_ind );
      }
  }
}

unsigned MatrixProductBase::getNumberOfColumns() const { 
  plumed_assert( !diagonal );
  return getPntrToOutput(0)->getShape()[1];
}

unsigned MatrixProductBase::getNumberOfInnerTasks() const {
  if( getPntrToOutput(0)->getRank()<2 ) return getPntrToArgument(0)->getShape()[1];
  return getPntrToOutput(0)->getShape()[1]; 
}

void MatrixProductBase::setupForTask( const unsigned& current, MultiValue& myvals, std::vector<unsigned> & indices, std::vector<Vector>& atoms ) const {
  if( diagonal && getPntrToOutput(0)->getRank()!=0 ) {
      if( indices.size()!=2 ) indices.resize( 2 ); 
      indices[0] = current; indices[1] = getFullNumberOfTasks() + current;
  } else {
      unsigned size_v = getNumberOfInnerTasks();
      if( skip_ieqj ) {
         if( indices.size()!=size_v ) indices.resize( size_v );
      } else if( indices.size()!=size_v+1 ) indices.resize( size_v + 1 );; 
      unsigned k=1; indices[0]=current; unsigned start_n = getFullNumberOfTasks();
      for(unsigned i=0; i<size_v; ++i) {
          if( skip_ieqj && myvals.getTaskIndex()==i ) continue;
          indices[k]=start_n + i; k++;
      }
  }
  myvals.setSplitIndex( indices.size() ); myvals.setNumberOfIndices( indices.size() );
}

void MatrixProductBase::performTask( const unsigned& current, MultiValue& myvals ) const {
  std::vector<unsigned> & indices( myvals.getIndices() );
  if( !doInnerLoop && actionInChain() ) {
    // If this is not an adjacency matrix then have done the relevant calculations during the first pass through the loop 
    if( !doNotCalculateDerivatives() && myvals.inVectorCall() ) updateCentralMatrixIndex( myvals.getTaskIndex(), indices, myvals );
    return ;
  }
  std::vector<Vector> & atoms( myvals.getFirstAtomVector() );
  setupForTask( current, myvals, indices, atoms );

  // Now loop over all atoms in coordination sphere
  unsigned ntwo_atoms = myvals.getSplitIndex();
  for(unsigned i=1; i<ntwo_atoms; ++i) {
    // This does everything in the stream that is done with single matrix elements
    runTask( getLabel(), myvals.getTaskIndex(), current, indices[i], myvals );
    // Now clear only elements that are not accumulated over whole row
    clearMatrixElements( myvals );
  }
  // Update the matrix index for the central atom
  if( !doNotCalculateDerivatives() ) updateCentralMatrixIndex( myvals.getTaskIndex(), indices, myvals );
}

void MatrixProductBase::updateAtomicIndices( const unsigned& index1, const unsigned& index2, MultiValue& myvals ) const {
  if( !doInnerLoop && getNumberOfComponents()>1 ) error("cannot calculate multiple outputs for this type of action");
  unsigned narg_derivatives = 0; if( getNumberOfArguments()>0 ) narg_derivatives = getPntrToArgument(0)->getSize() + getPntrToArgument(1)->getSize();
  unsigned w_ind = getPntrToOutput(0)->getPositionInStream();
  // Update dynamic list indices for central atom
  myvals.updateIndex( w_ind, narg_derivatives + 3*index1+0 ); myvals.updateIndex( w_ind, narg_derivatives + 3*index1+1 ); myvals.updateIndex( w_ind, narg_derivatives + 3*index1+2 );
  // Update dynamic list indices for atom forming this bond
  myvals.updateIndex( w_ind, narg_derivatives + 3*index2+0 ); myvals.updateIndex( w_ind, narg_derivatives + 3*index2+1 ); myvals.updateIndex( w_ind, narg_derivatives + 3*index2+2 );
  // Now look after all the atoms in the third block
  std::vector<unsigned> & indices( myvals.getIndices() );
  for(unsigned i=myvals.getSplitIndex(); i<myvals.getNumberOfIndices(); ++i) {
    myvals.updateIndex( w_ind, narg_derivatives + 3*indices[i]+0 ); myvals.updateIndex( w_ind, narg_derivatives + 3*indices[i]+1 ); myvals.updateIndex( w_ind, narg_derivatives + 3*indices[i]+2 );
  }
  // Update dynamic list indices for virial
  unsigned base = narg_derivatives + 3*getNumberOfAtoms(); for(unsigned j=0; j<9; ++j) myvals.updateIndex( w_ind, base+j );
  // Matrix indices
  if( !myvals.inMatrixRerun() ) {
      unsigned nmat = getPntrToOutput(0)->getPositionInMatrixStash(), nmat_ind = myvals.getNumberOfMatrixIndices( nmat );
      std::vector<unsigned>& matrix_indices( myvals.getMatrixIndices( nmat ) );
      matrix_indices[nmat_ind+0]=narg_derivatives + 3*index2+0; matrix_indices[nmat_ind+1]=narg_derivatives + 3*index2+1; matrix_indices[nmat_ind+2]=narg_derivatives + 3*index2+2;
      nmat_ind+=3; myvals.setNumberOfMatrixIndices( nmat, nmat_ind );
  }
}

bool MatrixProductBase::performTask( const std::string& controller, const unsigned& index1, const unsigned& index2, MultiValue& myvals ) const {
  // This makes sure other AdjacencyMatrixBase actions in the stream don't get their matrix elements calculated here
  unsigned noutput=getNumberOfComponents(); if( doInnerLoop ) { noutput=1; if( controller!=getLabel() ) return false; }
  // Now do the calculation
  unsigned ss0=0, ss1=0, nargs=0, ind2 = index2; if( index2>=getFullNumberOfTasks() ) ind2 = index2 - getFullNumberOfTasks();
  for(unsigned nv=0; nv<noutput; ++nv) {
      if( getNumberOfArguments()>0 ) {
          ss0=1; if( getPntrToArgument(nv)->getRank()==2 ) ss0=getPntrToArgument(nv)->getShape()[1];
          ss1=1; if( getPntrToArgument(noutput+nv)->getRank()==2 ) ss1=getPntrToArgument(noutput+nv)->getShape()[1];
          nargs=ss0; if( !diagonal && getPntrToOutput(nv)->getRank()<2 ) nargs=1;
      }
      std::vector<double> args1(nargs), args2(nargs), der1(nargs), der2(nargs);
      if( getPntrToOutput(nv)->getRank()==0 || (!diagonal && getPntrToOutput(nv)->getRank()<2) ) {
          if( actionInChain() ) {
              if( myvals.inMatrixRerun() ) return true;
              args1[0] = myvals.get( getPntrToArgument(nv)->getPositionInStream() );
          } else args1[0] = getPntrToArgument(nv)->get( index1*ss0 + ind2 );
          args2[0] = getPntrToArgument(noutput+nv)->get( ind2 );
      } else {
          for(unsigned i=0; i<nargs; ++i) {
            args1[i] = getPntrToArgument(nv)->get( index1*ss0 + i );
            args2[i] = getPntrToArgument(noutput+nv)->get( i*ss1 + ind2 );
          }
      }
      double val = computeVectorProduct( index1, index2, args1, args2, der1, der2, myvals );
      unsigned ostrn = getPntrToOutput(nv)->getPositionInStream();
      if( fabs(val)<epsilon && getPntrToOutput(nv)->getRank()==2 ) {
          if( !doNotCalculateDerivatives() ) {
              if( getNumberOfAtoms()>0 ) updateAtomicIndices( index1, index2, myvals );
              clearMatrixElements( myvals );
          }
          if( noutput==1 ) return false;    
          else continue;   
      }
      myvals.addValue( ostrn, val );
      // Return after calculation of value if we do not need derivatives
      if( doNotCalculateDerivatives() ) continue;

      if( getPntrToOutput(nv)->getRank()==0 || (!diagonal && getPntrToOutput(nv)->getRank()<2) ) {
          unsigned vstart = getPntrToArgument(nv)->getSize();
          if( actionInChain() ) {
              unsigned my_weight = getPntrToArgument(nv)->getPositionInStream();
              for(unsigned k=0; k<myvals.getNumberActive(my_weight); ++k) {
                  unsigned kind=myvals.getActiveIndex(my_weight,k);
                  myvals.addDerivative( ostrn, kind, der1[0]*myvals.getDerivative( my_weight, kind ) );
              }
              vstart = (getPntrToArgument(nv)->getPntrToAction())->getNumberOfDerivatives();
          } else { myvals.addDerivative( ostrn, index1*ss0 + ind2, der1[0] ); myvals.updateIndex( ostrn, index1*ss0 + ind2 ); }
          if( fabs(der2[0])>0 ) { myvals.addDerivative( ostrn, vstart + ind2, der2[0] ); myvals.updateIndex( ostrn, vstart + ind2 ); }
      } else {
          unsigned jind_start = 0; 
          if( getNumberOfArguments()>0 ) jind_start = getPntrToArgument(nv)->getSize(); 
          for(unsigned i=0; i<nargs; ++i) {
            plumed_dbg_assert( ss0*index1 + i<myvals.getNumberOfDerivatives() );
            myvals.addDerivative( ostrn, ss0*index1 + i, der1[i] );
            myvals.updateIndex( ostrn, ss0*index1 + i );
            plumed_dbg_assert( jind_start + i*ss1 + ind2<myvals.getNumberOfDerivatives() );
            myvals.addDerivative( ostrn, jind_start + i*ss1 + ind2, der2[i] );
            myvals.updateIndex( ostrn, jind_start + i*ss1 + ind2 );
          }
          if( !diagonal && !myvals.inMatrixRerun() ) {
              unsigned nmat = getPntrToOutput(nv)->getPositionInMatrixStash();
              std::vector<unsigned>& matrix_indices( myvals.getMatrixIndices( nmat ) );
              plumed_dbg_assert( matrix_indices.size()>=getNumberOfDerivatives() );
              unsigned nmat_ind = myvals.getNumberOfMatrixIndices( nmat );
              for(unsigned i=0; i<nargs; ++i) { matrix_indices[nmat_ind] = jind_start + i*ss1 + ind2; nmat_ind++; }
              myvals.setNumberOfMatrixIndices( nmat, nmat_ind );
          }
          if( getNumberOfAtoms()>0 ) updateAtomicIndices( index1, index2, myvals );
      }
  }
  return true;
}

void MatrixProductBase::apply() {
  if( doNotCalculateDerivatives() ) return;
  if( forcesToApply.size()!=getNumberOfDerivatives() ) forcesToApply.resize( getNumberOfDerivatives() );
  std::fill(forcesToApply.begin(),forcesToApply.end(),0); unsigned mm=0;
  if( getForcesFromValues( forcesToApply ) ) {
    setForcesOnAtoms( forcesToApply, mm );
    setForcesOnArguments( 0, forcesToApply, mm );
  }
}

}
}
