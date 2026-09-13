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
  keys.addFlag("AVOID_THREAD_GATHER",false,"avoid race conditions in shared memory by doing calculations rather than using more memory");
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
  parseFlag("AVOID_THREAD_GATHER",no_thread_gather);
  if( OpenMP::getNumThreads()==1 && getName().find("ACC")==std::string::npos ) {
      no_thread_gather = false;
  }
  if( no_thread_gather && comm.Get_size()>1 ) {
      error("AVOID_THREAD_GATHER keyword is incompatible with MPI - it should be possible to fix this. Email: gareth.tribello@gmail.com if you are interested");
  }
  if( no_thread_gather ) log.printf("  turning on low memory implementation for force gathering\n");
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
  Value* mycomp = getPntrToComponent(0);
  outmat.ncols = mycomp->getNumberOfColumns();
  outmat.resize( mycomp->matrix_bookeeping.size() );
  for(unsigned i=0; i<outmat.size(); ++i) {
    outmat[i] = mycomp->matrix_bookeeping[i];
  }
}

void ActionWithMatrix::updateBookeepingArrays( RequiredMatrixElements& outmat ) {
  RequiredMatrixElementsUpdater updater(outmat);
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
  findMaximumColumnLength();
}

void ActionWithMatrix::findMaximumColumnLength() {
  if( !no_thread_gather || doNotCalculateDerivatives() ) {
    return;
  }
  Value* mycomp = getPntrToComponent(0);
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

void ActionWithMatrix::getColumnBookeepingArrays( RequiredMatrixElements& outmat ) {
  outmat.ncols = maxcolsize;
  Value* mycomp = getPntrToComponent(0);
  outmat.resize( outmat.ncols*(1+mycomp->getShape()[1]) );
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

void ActionWithMatrix::transferStashToValues( const std::vector<unsigned>& partialTaskList, const std::vector<double>& stash ) {
  unsigned ncomp = getNumberOfComponents();
  unsigned ncols = getPntrToComponent(0)->getNumberOfColumns();
  unsigned nrows = partialTaskList.size();
  for(unsigned i=0; i<nrows; ++i) {
    unsigned ncr = getPntrToComponent(0)->getRowLength(partialTaskList[i]);
#ifndef NDEBUG
    for(unsigned k=1; k<ncomp; ++k) {
      plumed_assert( ncr == getPntrToComponent(k)->getRowLength(partialTaskList[i]) );
    }
#endif
    for(unsigned j=0; j<ncr; ++j) {
      for(unsigned k=0; k<ncomp; ++k) {
        getPntrToComponent(k)->set( partialTaskList[i]*ncols+j, stash[ncomp*ncols*partialTaskList[i]+j*ncomp+k] );
      }
    }
  }
}

void ActionWithMatrix::transferStashToValues( const std::vector<unsigned>& partialTaskList, const std::vector<float>& stash ) {
  unsigned ncomp = getNumberOfComponents();
  unsigned ncols = getPntrToComponent(0)->getNumberOfColumns();
  unsigned nrows = partialTaskList.size();
  for(unsigned i=0; i<nrows; ++i) {
    unsigned ncr = getPntrToComponent(0)->getRowLength(partialTaskList[i]);
#ifndef NDEBUG
    for(unsigned k=1; k<ncomp; ++k) {
      plumed_assert( ncr == getPntrToComponent(k)->getRowLength(partialTaskList[i]) );
    }
#endif
    for(unsigned j=0; j<ncr; ++j) {
      for(unsigned k=0; k<ncomp; ++k) {
        getPntrToComponent(k)->set( partialTaskList[i]*ncols+j, stash[ncomp*ncols*partialTaskList[i]+j*ncomp+k] );
      }
    }
  }
}

void ActionWithMatrix::transferForcesToStash( const std::vector<unsigned>& partialTaskList, std::vector<double>& stash ) const {
  unsigned ncomp = getNumberOfComponents();
  unsigned ncols = getConstPntrToComponent(0)->getNumberOfColumns();
  if( gatherForceOnColumns ) {
     const std::vector<unsigned>& rowTasks( getConstListOfActiveTasks() ); 
     unsigned nrows = rowTasks.size();
     std::vector<unsigned> column_totals( getConstPntrToComponent(0)->getShape()[1], 0 );
     for(unsigned i=0; i<nrows; ++i) {
       unsigned ncr = getConstPntrToComponent(0)->getRowLength(rowTasks[i]);
       for(unsigned j=0; j<ncr; ++j) {
           unsigned colno = getConstPntrToComponent(0)->getRowIndex(rowTasks[i],j);
           for(unsigned k=0; k<ncomp; ++k) {
               plumed_dbg_assert( colno==getConstPntrToComponent(k)->getRowIndex(rowTasks[i],j) );
               stash[ncomp*maxcolsize*colno + column_totals[colno]*ncomp + k] = getConstPntrToComponent(k)->getForce( rowTasks[i]*ncols+j );
           }
           column_totals[colno]++;
       }
     }
  } else {
     unsigned nrows = partialTaskList.size(); 
     for(unsigned i=0; i<nrows; ++i) {
       unsigned ncr = getConstPntrToComponent(0)->getRowLength(partialTaskList[i]);
#ifndef NDEBUG
       for(unsigned k=1; k<ncomp; ++k) {
         plumed_assert( ncr == getConstPntrToComponent(k)->getRowLength(partialTaskList[i]) );
       }
#endif
       for(unsigned j=0; j<ncr; ++j) {
         for(unsigned k=0; k<ncomp; ++k) {
           stash[ncomp*ncols*partialTaskList[i]+j*ncomp+k] = getConstPntrToComponent(k)->getForce( partialTaskList[i]*ncols+j );
         }
       }
     }
  }
}

void ActionWithMatrix::transferForcesToStash( const std::vector<unsigned>& partialTaskList, std::vector<float>& stash ) const {
  unsigned ncomp = getNumberOfComponents();
  unsigned ncols = getConstPntrToComponent(0)->getNumberOfColumns();
  if( gatherForceOnColumns ) {
     const std::vector<unsigned>& rowTasks( getConstListOfActiveTasks() );
     unsigned nrows = rowTasks.size(); 
     std::vector<unsigned> column_totals( getConstPntrToComponent(0)->getShape()[1], 0 );
     for(unsigned i=0; i<nrows; ++i) {
       unsigned ncr = getConstPntrToComponent(0)->getRowLength(rowTasks[i]);
       for(unsigned j=0; j<ncr; ++j) {
           unsigned colno = getConstPntrToComponent(0)->getRowIndex(rowTasks[i],j);
           for(unsigned k=0; k<ncomp; ++k) {
               plumed_dbg_assert( colno==getConstPntrToComponent(k)->getRowIndex(rowTasks[i],j) );
               stash[ncomp*maxcolsize*colno + column_totals[colno]*ncomp + k] = getConstPntrToComponent(k)->getForce( rowTasks[i]*ncols+j );
           }
           column_totals[colno]++;
       }
     }
  } else {
     unsigned nrows = partialTaskList.size();
     for(unsigned i=0; i<nrows; ++i) {
       unsigned ncr = getConstPntrToComponent(0)->getRowLength(partialTaskList[i]);
#ifndef NDEBUG
       for(unsigned k=1; k<ncomp; ++k) {
         plumed_assert( ncr == getConstPntrToComponent(k)->getRowLength(partialTaskList[i]) );
       }
#endif
       for(unsigned j=0; j<ncr; ++j) {
         for(unsigned k=0; k<ncomp; ++k) {
           stash[ncomp*ncols*partialTaskList[i]+j*ncomp+k] = getConstPntrToComponent(k)->getForce( partialTaskList[i]*ncols+j );
         }
       }
     }
  }
}

} //namespace PLMD
