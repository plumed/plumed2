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
#include "tools/OpenMP.h"

namespace PLMD {

void ActionWithMatrix::registerKeywords( Keywords& keys ) {
  ActionWithVector::registerKeywords( keys );
  keys.addFlag("AVOID_THREAD_GATHER",false,"avoid race conditions in shared memory by doing more calculations rather than using more memory");
  keys.addFlag("DO_THREAD_GATHER",false,"do fewer calculations and use more memory");
}

ActionWithMatrix::ActionWithMatrix(const ActionOptions&ao):
  Action(ao),
  ActionWithVector(ao),
  maxcolsize(0),
  diagzero(false),
  gatherForceOnColumns(false) {

  if( keywords.exists("ELEMENTS_ON_DIAGONAL_ARE_ZERO") ) {
    parseFlag("ELEMENTS_ON_DIAGONAL_ARE_ZERO",diagzero);
    if( diagzero ) {
      log.printf("  setting diagonal elements equal to zero\n");
    }
  }
  if( getName().find("ACC")!=std::string::npos ) {
    bool ignore;
    parseFlag("AVOID_THREAD_GATHER",ignore);
    if( ignore ) {
      warning("AVOID_THREAD_GATHER flag is always activated when the calculation is run on the GPU as the calculation is much faster with this option activated");
    }
    bool dothreads;
    parseFlag("DO_THREAD_GATHER",dothreads);
    if( !dothreads ) {
      no_thread_gather=true;
    } else {
      no_thread_gather=false;
    }
  } else {
    parseFlag("AVOID_THREAD_GATHER",no_thread_gather);
    if( OpenMP::getNumThreads()==1 ) {
      no_thread_gather = false;
    }
    bool ignore;
    parseFlag("DO_THREAD_GATHER",ignore);
    if( ignore ) {
      warning("NO_THREAD_GATHER flag is always activated when the calculation is run on the CPU as the calculation is faster with this option activated");
    }
  }
  if( no_thread_gather ) {
    log.printf("  turning on low memory implementation for force gathering\n");
  }
}

void ActionWithMatrix::checkThreadGatherWithMPI( const bool serial ) {
  if( no_thread_gather && comm.Get_size()>1 && !serial ) {
    if( getName().find("ACC")!=std::string::npos ) {
      error("with more than one MPI process the GPU version needs DO_THREAD_GATHER (the low memory gathering of the forces, used by default, is not parallelised with MPI) or SERIAL");
    }
    error("AVOID_THREAD_GATHER cannot be used with more than one MPI process unless SERIAL is also used (the low memory gathering of the forces is not parallelised with MPI)");
  }
}


class RequiredMatrixElementsUpdater {
  RequiredMatrixElements& outmat;
public:
  RequiredMatrixElementsUpdater( RequiredMatrixElements& mat ) : outmat(mat) {}
  ~RequiredMatrixElementsUpdater() {
    outmat.update();
  }
};

void ActionWithMatrix::copyMatrixBookeepingFromFirstComponent( RequiredMatrixElements& outmat ) {
  RequiredMatrixElementsUpdater updater(outmat);
  Value* mycomp = getPntrToComponent(0);
  outmat.ncols = mycomp->getNumberOfColumns();
  outmat.resize( mycomp->matrix_bookeeping.size() );
  for(unsigned i=0; i<outmat.size(); ++i) {
    outmat[i] = mycomp->matrix_bookeeping[i];
  }
  if( !no_thread_gather || doNotCalculateDerivatives() ) {
    return;
  }
  std::vector<unsigned> column_totals( getConstPntrToComponent(0)->getShape()[1], 0 );
  for(unsigned i=0; i<mycomp->getShape()[0]; ++i) {
    for(unsigned j=0; j<mycomp->getRowLength(i); ++j) {
      column_totals[ mycomp->getRowIndex(i,j) ]++;
    }
  }
  unsigned n_nonzerocols = 0;
  maxcolsize = column_totals[0];
  if( maxcolsize>0 ) {
    n_nonzerocols = 1;
  }
  for(unsigned i=1; i<column_totals.size(); ++i) {
    if( column_totals[i]>maxcolsize ) {
      maxcolsize = column_totals[i];
    }
    if( column_totals[i]>0 ) {
      n_nonzerocols++;
    }
  }
  column_list.resize( n_nonzerocols );
  n_nonzerocols = 0;
  for(unsigned i=0; i<column_totals.size(); ++i) {
    if( column_totals[i]>0 ) {
      column_list[n_nonzerocols]=i;
      n_nonzerocols++;
    }
  }
}

void ActionWithMatrix::updateBookeepingArrays( RequiredMatrixElements& outmat ) {
  Value* myval = getPntrToComponent(0);
  unsigned lstart = myval->getShape()[0];
  if( getNumberOfMasks()>0 ) {
    Value* maskarg = getPntrToArgument(getNumberOfArguments()-getNumberOfMasks());
    for(unsigned i=0; i<getNumberOfComponents(); ++i) {
      getPntrToComponent(i)->reshapeMatrixStore( maskarg->getNumberOfColumns() );
    }
    for(unsigned i=0; i<lstart; ++i) {
      unsigned rstart = i*(1+myval->getNumberOfColumns());
      myval->setMatrixBookeepingElement( rstart, maskarg->getRowLength(i) );
      for(unsigned j=0; j<maskarg->getRowLength(i); ++j) {
        myval->setMatrixBookeepingElement( rstart + 1 + j, maskarg->getRowIndex(i, j) );
      }
    }
  } else if ( diagzero ) {
    for(unsigned i=0; i<getNumberOfArguments(); ++i) {
      if( getPntrToArgument(i)->getRank()==2 && getPntrToArgument(i)->getShape()[1]!=getPntrToArgument(i)->getNumberOfColumns() ) {
        error("DIAGZERO flag is incompatible with sparse matrices");
      }
    }
    for(unsigned i=0; i<getNumberOfComponents(); ++i) {
      getPntrToComponent(i)->reshapeMatrixStore( myval->getShape()[1]-1 );
    }
    for(unsigned i=0; i<lstart; ++i) {
      unsigned k=0, rstart = i*(1+myval->getNumberOfColumns());
      myval->setMatrixBookeepingElement( rstart, myval->getShape()[1]-1 );
      for(unsigned j=0; j<myval->getShape()[1]; ++j) {
        if( i!=j ) {
          myval->setMatrixBookeepingElement( rstart + 1 + k, j );
          k++;
        }
      }
    }
  } else {
    for(unsigned i=0; i<getNumberOfComponents(); ++i) {
      getPntrToComponent(i)->reshapeMatrixStore( myval->getShape()[1] );
    }
  }
  copyMatrixBookeepingFromFirstComponent( outmat );
  for(unsigned i=1; i<getNumberOfComponents(); ++i) {
    getPntrToComponent(i)->copyBookeepingArrayFromArgument( myval );
  }
}

void ActionWithMatrix::getColumnBookeepingArrays( RequiredMatrixElements& outmat ) {
  outmat.ncols = maxcolsize;
  Value* mycomp = getPntrToComponent(0);
  outmat.resize( mycomp->getShape()[1]*(1+outmat.ncols) );
  for(unsigned i=0; i<mycomp->getShape()[1]; ++i) {
    outmat[i*(1+outmat.ncols)] = 0;
  }
  for(unsigned i=0; i<mycomp->getShape()[0]; ++i) {
    for(unsigned j=0; j<mycomp->getRowLength(i); ++j) {
      unsigned cstart = mycomp->getRowIndex(i,j)*( 1 + outmat.ncols );
      outmat[cstart + 1 + outmat[cstart]] = i;
      outmat[cstart]++;
    }
  }
}

std::vector<unsigned>& ActionWithMatrix::getListOfActiveTasks() {
  if( gatherForceOnColumns ) {
    return column_list;
  } else {
    return ActionWithVector::getListOfActiveTasks();
  }
}

unsigned ActionWithMatrix::getSizeOfValueInStash( const Value* myval ) {
  unsigned nstored = myval->getNumberOfStoredValues();
  if( no_thread_gather ) {
    unsigned ncolstore = maxcolsize*myval->getShape()[1];
    if( ncolstore>nstored ) {
      return ncolstore;
    }
  }
  return nstored;
}

// The stash holds the elements of a row contiguously, with the components of
// an element interleaved ([row][column][component]); each component Value
// stores its rows contiguously ([row][column]). The copies write the storage
// of the Values directly (ActionWithMatrix is a friend of Value), the rows
// are independent so the loop over them is parallel.
template <typename T>
void ActionWithMatrix::matrixStashToValuesImpl( const std::vector<unsigned>& partialTaskList, const std::vector<T>& stash ) {
  const std::size_t ncomp = getNumberOfComponents();
  const std::size_t ncols = getPntrToComponent(0)->getNumberOfColumns();
  const unsigned nrows = partialTaskList.size();
  std::vector<double*> data( ncomp );
  std::vector<Value*> periodic;
  for(std::size_t k=0; k<ncomp; ++k) {
    Value* myval = getPntrToComponent(k);
    plumed_dbg_assert( !myval->hasDeriv );
    myval->value_set=true;
    data[k] = myval->data.data();
    if( myval->isPeriodic() ) {
      periodic.push_back( myval );
    }
  }
  const Value* first = getPntrToComponent(0);
#ifdef _OPENMP
  // without OpenMP nt would be unused (an error with -Werror)
  const unsigned nt = OpenMP::getGoodNumThreads( stash );
#endif
  #pragma omp parallel for num_threads(nt) schedule(static)
  for(unsigned i=0; i<nrows; ++i) {
    const std::size_t row = partialTaskList[i];
    const unsigned ncr = first->getRowLength(row);
#ifndef NDEBUG
    for(std::size_t k=1; k<ncomp; ++k) {
      plumed_assert( ncr == getPntrToComponent(k)->getRowLength(row) );
    }
#endif
    const T* src = stash.data() + ncomp*ncols*row;
    const std::size_t dst = row*ncols;
    for(std::size_t j=0; j<ncr; ++j) {
      for(std::size_t k=0; k<ncomp; ++k) {
        data[k][dst+j] = src[j*ncomp+k];
      }
    }
    // e.g. the matrix of torsions is periodic
    for(Value* myval : periodic) {
      for(std::size_t j=0; j<ncr; ++j) {
        myval->applyPeriodicity( dst+j );
      }
    }
  }
}

template <typename T>
void ActionWithMatrix::matrixForcesToStashImpl( const std::vector<unsigned>& partialTaskList, std::vector<T>& stash ) const {
  const std::size_t ncomp = getNumberOfComponents();
  const std::size_t ncols = getConstPntrToComponent(0)->getNumberOfColumns();
  const Value* first = getConstPntrToComponent(0);
  std::vector<const double*> force( ncomp );
  for(std::size_t k=0; k<ncomp; ++k) {
    force[k] = getConstPntrToComponent(k)->inputForce.data();
  }
  if( gatherForceOnColumns ) {
    // forces on the columns: element (i,j) of the matrix goes to the slot of
    // column j; the slots of a column are filled in the order of the rows
    const std::vector<unsigned>& rowTasks( getConstListOfActiveTasks() );
    const unsigned nrows = rowTasks.size();
    std::vector<unsigned> column_totals( first->getShape()[1], 0 );
    for(unsigned i=0; i<nrows; ++i) {
      const std::size_t row = rowTasks[i];
      const unsigned ncr = first->getRowLength(row);
      for(std::size_t j=0; j<ncr; ++j) {
        const std::size_t colno = first->getRowIndex(row,j);
        T* dst = stash.data() + ncomp*maxcolsize*colno + column_totals[colno]*ncomp;
        for(std::size_t k=0; k<ncomp; ++k) {
          plumed_dbg_assert( colno==getConstPntrToComponent(k)->getRowIndex(row,j) );
          dst[k] = force[k][row*ncols+j];
        }
        column_totals[colno]++;
      }
    }
    return;
  }
  const unsigned nrows = partialTaskList.size();
#ifdef _OPENMP
  // without OpenMP nt would be unused (an error with -Werror)
  const unsigned nt = OpenMP::getGoodNumThreads( stash );
#endif
  #pragma omp parallel for num_threads(nt) schedule(static)
  for(unsigned i=0; i<nrows; ++i) {
    const std::size_t row = partialTaskList[i];
    const unsigned ncr = first->getRowLength(row);
#ifndef NDEBUG
    for(std::size_t k=1; k<ncomp; ++k) {
      plumed_assert( ncr == getConstPntrToComponent(k)->getRowLength(row) );
    }
#endif
    T* dst = stash.data() + ncomp*ncols*row;
    const std::size_t src = row*ncols;
    for(std::size_t j=0; j<ncr; ++j) {
      for(std::size_t k=0; k<ncomp; ++k) {
        dst[j*ncomp+k] = force[k][src+j];
      }
    }
  }
}

void ActionWithMatrix::transferStashToValues( const std::vector<unsigned>& partialTaskList, const std::vector<double>& stash ) {
  matrixStashToValuesImpl( partialTaskList, stash );
}

void ActionWithMatrix::transferStashToValues( const std::vector<unsigned>& partialTaskList, const std::vector<float>& stash ) {
  matrixStashToValuesImpl( partialTaskList, stash );
}

void ActionWithMatrix::transferForcesToStash( const std::vector<unsigned>& partialTaskList, std::vector<double>& stash ) const {
  matrixForcesToStashImpl( partialTaskList, stash );
}

void ActionWithMatrix::transferForcesToStash( const std::vector<unsigned>& partialTaskList, std::vector<float>& stash ) const {
  matrixForcesToStashImpl( partialTaskList, stash );
}

} //namespace PLMD
