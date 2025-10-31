/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2020 The plumed team
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
#ifndef __PLUMED_function_FunctionOfMatrix_h
#define __PLUMED_function_FunctionOfMatrix_h

#include "core/ActionWithVector.h"
#include "core/ParallelTaskManager.h"
#include "FunctionSetup.h"
#include "FunctionOfVector.h"
#include "Custom.h"

namespace PLMD {
namespace function {

template <class T>
class FunctionOfMatrix : public ActionWithVector {
public:
  using input_type = FunctionData<T>;
  using PTM = ParallelTaskManager<FunctionOfMatrix<T>>;
private:
/// The parallel task manager
  PTM taskmanager;
/// Set equal to one if we are doing EvaluateGridFunction
  unsigned argstart;
/// Number of scalars that appear in the input
  std::size_t nscalars;
/// Used to hold the list of tasks we are running
  std::vector<unsigned> active_tasks;
/// Get the number of arguments the function uses
  unsigned getNumberOfFunctionArguments() const ;
public:
  static void registerKeywords(Keywords&);
  explicit FunctionOfMatrix(const ActionOptions&);
/// Get the label to write in the graph
  std::string writeInGraph() const override ;
/// Get the number of derivatives for this action
  unsigned getNumberOfDerivatives() override ;
/// Resize the matrices
  void prepare() override ;
  void calculate() override ;
  void getNumberOfTasks( unsigned& ntasks ) override ;
  std::vector<unsigned>& getListOfActiveTasks( ActionWithVector* action ) override ;
  void getInputData( std::vector<double>& inputdata ) const override ;
  static void performTask( std::size_t task_index,
                           const FunctionData<T>& actiondata,
                           ParallelActionsInput& input,
                           ParallelActionsOutput& output );
/// Add some forces
  void applyNonZeroRankForces( std::vector<double>& outforces ) override ;
/// Get the indices of the forces
  static int getNumberOfValuesPerTask( std::size_t task_index, const FunctionData<T>& actiondata );
  static void getForceIndices( std::size_t task_index,
                               std::size_t colno,
                               std::size_t ntotal_force,
                               const FunctionData<T>& actiondata,
                               const ParallelActionsInput& input,
                               ForceIndexHolder force_indices );
};

template <class T>
void FunctionOfMatrix<T>::registerKeywords(Keywords& keys ) {
  ActionWithVector::registerKeywords(keys);
  std::string name = keys.getDisplayName();
  std::size_t und=name.find("_MATRIX");
  keys.setDisplayName( name.substr(0,und) );
  keys.addInputKeyword("compulsory","ARG","scalar/matrix","the labels of the scalar and matrices that on which the function is being calculated elementwise");
  keys.addInputKeyword("optional","MASK","matrix","a matrix that is used to used to determine which elements of the output matrix to compute");
  keys.add("hidden","NO_ACTION_LOG","suppresses printing from action on the log");
  keys.reserve("compulsory","PERIODIC","if the output of your function is periodic then you should specify the periodicity of the function.  If the output is not periodic you must state this using PERIODIC=NO");
  T::registerKeywords( keys );
  if( keys.getDisplayName()=="SUM" ) {
    keys.setValueDescription("scalar","the sum of all the elements in the input matrix");
  } else if( keys.getDisplayName()=="HIGHEST" ) {
    keys.setValueDescription("scalar","the largest element of the input matrix");
  } else if( keys.getDisplayName()=="LOWEST" ) {
    keys.setValueDescription("scalar","the smallest element in the input matrix");
  } else if( keys.outputComponentExists(".#!value") ) {
    keys.setValueDescription("matrix","the matrix obtained by doing an element-wise application of " + keys.getOutputComponentDescription(".#!value") + " to the input matrix");
  }
  PTM::registerKeywords( keys );
}

template <class T>
unsigned FunctionOfMatrix<T>::getNumberOfFunctionArguments() const {
  unsigned nargs=getNumberOfArguments();
  if( getNumberOfMasks()>0 ) {
    return nargs - getNumberOfMasks();
  }
  return nargs;
}

template <class T>
FunctionOfMatrix<T>::FunctionOfMatrix(const ActionOptions&ao):
  Action(ao),
  ActionWithVector(ao),
  taskmanager(this),
  argstart(0),
  nscalars(0) {
  // Check if first argument is grid
  if( getPntrToArgument(0)->getRank()>0 && getPntrToArgument(0)->hasDerivatives() ) {
    argstart=1;
  }
  if( getNumberOfArguments()==argstart ) {
    error("no arguments specified");
  }

  if( getPntrToArgument(argstart)->getRank()!=2 ) {
    error("first argument to this action must be a matrix");
  }

  // Get the number of arguments
  unsigned nargs = getNumberOfArguments();
  int nmasks = getNumberOfMasks();
  if( nargs>=static_cast<unsigned>(nmasks) && nmasks>0 ) {
    nargs = nargs - nmasks;
  }
  // Get the shape of the output
  std::vector<std::size_t> shape( 2 );
  shape[0] = getPntrToArgument(argstart)->getShape()[0];
  shape[1] = getPntrToArgument(argstart)->getShape()[1];
  for(unsigned i=argstart+1; i<nargs; ++i) {
    if( getPntrToArgument(i)->getRank()==0 ) {
      nscalars++;
    } else if( getPntrToArgument(i)->getRank()==2 ) {
      if( getPntrToArgument(i)->getShape()[0]!=shape[0] || getPntrToArgument(i)->getShape()[1]!=shape[1] ) {
        error("mismatch between sizes of input arguments");
      } else if( nscalars>0 ) {
        error("scalars should be specified in argument list after all matrices");
      }
    } else {
      error("input arguments should be matrices or scalars");
    }
  }
  if( nmasks>0 ) {
    if( getPntrToArgument(getNumberOfArguments()-nmasks)->getShape()[0]!=shape[0] ||
        getPntrToArgument(getNumberOfArguments()-nmasks)->getShape()[1]!=shape[1] ) {
      error("input mask has wrong size");
    }
  }

  // Check if the output matrix is symmetric
  bool symmetric=true;
  for(unsigned i=argstart; i<getNumberOfArguments(); ++i) {
    if( getPntrToArgument(i)->getRank()==2 ) {
      if( !getPntrToArgument(i)->isSymmetric() ) {
        symmetric=false;
      }
    }
  }
  // Setup the values
  // Get the function data from the parallel task manager, to avoid copies
  auto & myfunc = taskmanager.getActionInput();
  myfunc.argstart = argstart;
  myfunc.nscalars = nscalars;
  FunctionData<T>::setup( myfunc.f, keywords.getOutputComponents(), shape, false, this );
  // Copy the fact that this is a symmetric matrix if the input matrices are all symmetric
  for(unsigned i=0; i<getNumberOfComponents(); ++i) {
    getPntrToComponent(i)->setSymmetric( symmetric );
  }
  taskmanager.setupParallelTaskManager( getNumberOfFunctionArguments() - argstart,
                                        nscalars );
}

template <class T>
std::string FunctionOfMatrix<T>::writeInGraph() const {
  std::size_t und = getName().find_last_of("_");
  return getName().substr(0,und);
}

template <class T>
unsigned FunctionOfMatrix<T>::getNumberOfDerivatives() {
  unsigned nder=0;
  for(unsigned i=argstart; i<getNumberOfFunctionArguments(); ++i) {
    nder += getPntrToArgument(i)->getNumberOfStoredValues();
  }
  return nder;
}

template <class T>
void FunctionOfMatrix<T>::prepare() {
  std::vector<std::size_t> shape(getPntrToArgument(argstart)->getShape());
  for(unsigned i=0; i<getNumberOfComponents(); ++i) {
    Value* myval = getPntrToComponent(i);
    if( myval->getRank()==2 && (myval->getShape()[0]!=shape[0] || myval->getShape()[1]!=shape[1]) ) {
      myval->setShape(shape);
    }
  }
  ActionWithVector::prepare();
  active_tasks.resize(0);
}

template <class T>
void FunctionOfMatrix<T>::getNumberOfTasks( unsigned& ntasks ) {
  ntasks=getPntrToComponent(0)->getNumberOfStoredValues();
}

template <class T>
std::vector<unsigned>& FunctionOfMatrix<T>::getListOfActiveTasks( ActionWithVector* action ) {
  if( active_tasks.size()>0 ) {
    return active_tasks;
  }

  Value* myarg = NULL;
  if( getNumberOfMasks()>0 ) {
    myarg = getPntrToArgument(getNumberOfArguments()-getNumberOfMasks());
  } else {
    myarg = getPntrToArgument(argstart);
  }
  unsigned atsize = 0;
  unsigned nrows = myarg->getShape()[0];
  for(unsigned i=0; i<nrows; ++i) {
    atsize += myarg->getRowLength(i);
  }
  active_tasks.resize( atsize );

  for(unsigned i=0, base=0,k=0; i<nrows; ++i) {
    unsigned ncols = myarg->getRowLength(i);
    for(unsigned j=0; j<ncols; ++j) {
      active_tasks[k] = base+j;
      ++k;
    }
    base += myarg->getNumberOfColumns();
  }
  if( getNumberOfMasks()>0 && doNotCalculateDerivatives() ) {
    return active_tasks;
  }
// All the matrices input have to have the same sparsity pattern.
// I can do everything I want to do with this limitation.  If
// anyone wants to make this smarter in the future they can
#ifndef DNDEBUG
  for(unsigned k=argstart; k<getNumberOfArguments(); ++k) {
    if( getPntrToArgument(k)->getRank()!=2 ) {
      continue ;
    }
    if( getNumberOfMasks()>0 && getPntrToArgument(k)->isConstant() ) {
      continue ;
    }
    for(unsigned i=0; i<nrows; ++i) {
      unsigned ncols = myarg->getRowLength(i);
      if( getNumberOfMasks()>0 && ncols==0 ) {
        continue;
      }
      plumed_massert( ncols==getPntrToArgument(k)->getRowLength(i), "failing in " + getLabel() );
      for(unsigned j=0; j<ncols; ++j) {
        plumed_assert( myarg->getRowIndex(i,j)==getPntrToArgument(k)->getRowIndex(i,j) );
      }
    }
  }
#endif
  return active_tasks;
}

template <class T>
void FunctionOfMatrix<T>::getInputData( std::vector<double>& inputdata ) const {
  int nmasks = getNumberOfMasks();
  unsigned nargs = getNumberOfFunctionArguments();

  const Value* myval = getConstPntrToComponent(0);
  std::size_t ntasks = myval->getNumberOfStoredValues();
  std::size_t ndata = static_cast<std::size_t>(nargs-argstart)*ntasks;
  if( inputdata.size()!=ndata ) {
    inputdata.resize( ndata );
  }

  for(unsigned j=argstart; j<nargs; ++j) {
    const Value* jarg =  getPntrToArgument(j);
    if( jarg->getRank()==0 ) {
      double val = jarg->get();
      for(unsigned i=0; i<myval->getShape()[0]; ++i) {
        unsigned colbase=i*myval->getNumberOfColumns();
        for(unsigned k=0; k<myval->getRowLength(i); ++k) {
          inputdata[(nargs-argstart)*(colbase+k) + j-argstart] = val;
        }
      }
    } else if( nmasks>0 ) {
      for(unsigned i=0; i<myval->getShape()[0]; ++i) {
        unsigned jcolbase = i*jarg->getShape()[1];
        unsigned vcolbase = i*myval->getNumberOfColumns();
        for(unsigned k=0; k<myval->getRowLength(i); ++k) {
          inputdata[(nargs-argstart)*(vcolbase+k) + j-argstart] = jarg->get(jcolbase+myval->getRowIndex(i,k),true);
        }
      }
    } else {
      for(unsigned i=0; i<jarg->getShape()[0]; ++i) {
        unsigned colbase=i*jarg->getNumberOfColumns();
        for(unsigned k=0; k<jarg->getRowLength(i); ++k) {
          inputdata[(nargs-argstart)*(colbase+k) + j-argstart] = jarg->get(colbase+k,false);
        }
      }
    }
  }
}

template <class T>
void FunctionOfMatrix<T>::calculate() {
  Value* myarg = NULL;
  if( getNumberOfMasks()>0 ) {
    myarg = getPntrToArgument(getNumberOfArguments()-getNumberOfMasks());
  } else {
    myarg = getPntrToArgument(argstart);
  }
  // Copy bookeeping arrays from input matrices to output matrices
  for(unsigned i=0; i<getNumberOfComponents(); ++i) {
    getPntrToComponent(i)->copyBookeepingArrayFromArgument( myarg );
  }
  taskmanager.setupParallelTaskManager( getNumberOfFunctionArguments()-argstart, nscalars );
  taskmanager.runAllTasks();
}

template <class T>
void FunctionOfMatrix<T>::performTask( std::size_t task_index,
                                       const FunctionData<T>& actiondata,
                                       ParallelActionsInput& input,
                                       ParallelActionsOutput& output ) {
  auto funcout = FunctionOutput::create( input.ncomponents,
                                         output.values.data(),
                                         input.nderivatives_per_scalar,
                                         output.derivatives.data() );
  T::calc( actiondata.f,
           input.noderiv,
           View<const double>( input.inputdata + task_index*input.nderivatives_per_scalar,
                               input.nderivatives_per_scalar ),
           funcout );
}

template <class T>
void FunctionOfMatrix<T>::applyNonZeroRankForces( std::vector<double>& outforces ) {
  taskmanager.applyForces( outforces );
}

template <class T>
int FunctionOfMatrix<T>::getNumberOfValuesPerTask( std::size_t task_index, const FunctionData<T>& actiondata ) {
  return 1;
}

template <class T>
void FunctionOfMatrix<T>::getForceIndices( std::size_t task_index,
    std::size_t colno,
    std::size_t ntotal_force,
    const FunctionData<T>& actiondata,
    const ParallelActionsInput& input,
    ForceIndexHolder force_indices ) {
  // The force indices are found in the same way as in FunctionOfVector so we reuse that function here
  FunctionOfVector<T>::getForceIndices( task_index, colno, ntotal_force, actiondata, input, force_indices );
}

}
}
#endif
