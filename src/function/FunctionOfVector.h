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
#ifndef __PLUMED_function_FunctionOfVector_h
#define __PLUMED_function_FunctionOfVector_h

#include "core/ActionWithVector.h"
#include "core/ParallelTaskManager.h"
#include "core/ActionSetup.h"
#include "FunctionSetup.h"
#include "Custom.h"

namespace PLMD {
namespace function {

template <class T>
class FunctionOfVector : public ActionWithVector {
public:
  using input_type = FunctionData<T>;
  using PTM = ParallelTaskManager<FunctionOfVector<T>>;
private:
/// The parallel task manager
  PTM taskmanager;
/// Set equal to one if we are doing EvaluateGridFunction
  unsigned argstart;
public:
  static void registerKeywords(Keywords&);
  explicit FunctionOfVector(const ActionOptions&);
  ~FunctionOfVector() {}
  std::string getOutputComponentDescription( const std::string& cname, const Keywords& keys ) const override ;
/// Get the number of derivatives for this action
  unsigned getNumberOfDerivatives() override ;
/// Resize vectors that are the wrong size
  void prepare() override ;
/// Get the label to write in the graph
  std::string writeInGraph() const override ;
/// This builds the task list for the action
  void calculate() override;
/// Add some forces
  void applyNonZeroRankForces( std::vector<double>& outforces ) override ;
/// Get the input data
  void getInputData( std::vector<double>& inputdata ) const override ;
/// Calculate the function
  static void performTask( std::size_t task_index,
                           const FunctionData<T>& actiondata,
                           ParallelActionsInput& input,
                           ParallelActionsOutput& output );
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
void FunctionOfVector<T>::registerKeywords(Keywords& keys ) {
  Action::registerKeywords(keys);
  ActionWithValue::registerKeywords(keys);
  ActionWithArguments::registerKeywords(keys);
  std::string name = keys.getDisplayName();
  std::size_t und=name.find("_VECTOR");
  keys.setDisplayName( name.substr(0,und) );
  keys.addInputKeyword("compulsory","ARG","scalar/vector","the labels of the scalar and vector that on which the function is being calculated elementwise");
  keys.reserve("compulsory","PERIODIC","if the output of your function is periodic then you should specify the periodicity of the function.  If the output is not periodic you must state this using PERIODIC=NO");
  keys.add("hidden","NO_ACTION_LOG","suppresses printing from action on the log");
  T::registerKeywords( keys );
  if( keys.getDisplayName()=="SUM" ) {
    keys.setValueDescription("scalar","the sum of all the elements in the input vector");
  } else if( keys.getDisplayName()=="MEAN" ) {
    keys.setValueDescription("scalar","the mean of all the elements in the input vector");
  } else if( keys.getDisplayName()=="HIGHEST" ) {
    keys.setValueDescription("scalar/vector","the largest element of the input vector if one vector specified.  If multiple vectors of the same size specified the largest elements of these vector computed elementwise.");
    keys.addInputKeyword("optional","MASK","vector","the label for a sparse vector that should be used to determine which elements of the vector should be computed");
  } else if( keys.getDisplayName()=="LOWEST" ) {
    keys.setValueDescription("scalar/vector","the smallest element in the input vector if one vector specified.  If multiple vectors of the same size specified the largest elements of these vector computed elementwise.");
    keys.addInputKeyword("optional","MASK","vector","the label for a sparse vector that should be used to determine which elements of the vector should be computed");
  } else if( keys.getDisplayName()=="SORT" ) {
    keys.setValueDescription("vector","a vector that has been sorted into ascending order");
    keys.addInputKeyword("optional","MASK","vector","the label for a sparse vector that should be used to determine which elements of the vector should be computed");
  } else if( keys.outputComponentExists(".#!value") ) {
    keys.addInputKeyword("optional","MASK","vector","the label for a sparse vector that should be used to determine which elements of the vector should be computed");
    keys.setValueDescription("vector","the vector obtained by doing an element-wise application of " + keys.getOutputComponentDescription(".#!value") + " to the input vectors");
  }
  PTM::registerKeywords( keys );
}

template <class T>
FunctionOfVector<T>::FunctionOfVector(const ActionOptions&ao):
  Action(ao),
  ActionWithVector(ao),
  taskmanager(this),
  argstart(0) {
  // Check if first argument is grid
  if( getPntrToArgument(0)->getRank()>0 && getPntrToArgument(0)->hasDerivatives() ) {
    argstart=1;
  }
  if( getNumberOfArguments()==argstart ) {
    error("no arguments specified");
  }

  if( getPntrToArgument(argstart)->getRank()!=1 ) {
    error("first argument to this action must be a vector");
  }

  // Get the number of arguments
  unsigned nargs = getNumberOfArguments();
  int nmasks = getNumberOfMasks();
  if( nargs>=static_cast<unsigned>(nmasks) && nmasks>0 ) {
    nargs = nargs - nmasks;
  }

  // Get the shape of the output
  std::size_t nscalars = 0;
  std::vector<std::size_t> shape(1);
  shape[0]=getPntrToArgument(argstart)->getShape()[0];
  for(unsigned i=argstart+1; i<nargs; ++i) {
    if( getPntrToArgument(i)->getRank()==0 ) {
      nscalars++;
    } else if( getPntrToArgument(i)->getRank()==1 ) {
      if( getPntrToArgument(i)->getShape()[0]!=shape[0] ) {
        error("mismatch between sizes of input arguments");
      } else if( nscalars>0 ) {
        error("scalars should be specified in argument list after all vectors");
      }
    } else {
      error("input arguments should be vectors or scalars");
    }
  }
  if( nmasks>0 && getPntrToArgument(getNumberOfArguments()-nmasks)->getShape()[0]!=shape[0] ) {
    error("input mask has wrong size");
  }

  // Setup the function and the values
  // Get the function data from the parallel task manager, to avoid copies
  auto & myfunc = taskmanager.getActionInput();
  myfunc.argstart = argstart;
  myfunc.nscalars = nscalars;
  FunctionData<T>::setup( myfunc.f, keywords.getOutputComponents(), shape, false, this );
  // Setup the parallel task manager
  taskmanager.setupParallelTaskManager( nargs-argstart, nscalars );
}

template <class T>
std::string FunctionOfVector<T>::getOutputComponentDescription( const std::string& cname, const Keywords& keys ) const {
  if( getName().find("SORT")==std::string::npos ) {
    return ActionWithValue::getOutputComponentDescription( cname, keys );
  }
  return "the " + cname + "th largest element in the input vectors";
}

template <class T>
std::string FunctionOfVector<T>::writeInGraph() const {
  std::size_t und = getName().find_last_of("_");
  return getName().substr(0,und);
}

template <class T>
unsigned FunctionOfVector<T>::getNumberOfDerivatives() {
  unsigned nder = 0;
  for(unsigned i=argstart; i<getNumberOfArguments(); ++i) {
    nder += getPntrToArgument(i)->getNumberOfStoredValues();
  }
  return nder;
}

template <class T>
void FunctionOfVector<T>::prepare() {
  std::vector<std::size_t> shape(1);
  for(unsigned i=argstart; i<getNumberOfArguments(); ++i) {
    if( getPntrToArgument(i)->getRank()==1 ) {
      shape[0] = getPntrToArgument(i)->getShape()[0];
      break;
    }
  }
  for(unsigned i=0; i<getNumberOfComponents(); ++i) {
    Value* myval = getPntrToComponent(i);
    if( myval->getRank()==1 && myval->getShape()[0]!=shape[0] ) {
      myval->setShape(shape);
    }
  }
  ActionWithVector::prepare();
}

template <class T>
void FunctionOfVector<T>::getInputData( std::vector<double>& inputdata ) const {
  unsigned nargs = getNumberOfArguments();
  int nmasks = getNumberOfMasks();
  if( nargs>=static_cast<unsigned>(nmasks) && nmasks>0 ) {
    nargs = nargs - nmasks;
  }

  std::size_t ntasks = 0;
  for(unsigned i=argstart; i<nargs; ++i) {
    if( getPntrToArgument(i)->getRank()==1 ) {
      ntasks = getPntrToArgument(i)->getShape()[0];
      break;
    }
  }

  std::size_t ndata = static_cast<std::size_t>(nargs-argstart)*ntasks;
  if( inputdata.size()!=ndata ) {
    inputdata.resize( ndata );
  }

  for(unsigned j=argstart; j<nargs; ++j) {
    const Value* myarg =  getPntrToArgument(j);
    if( myarg->getRank()==0 ) {
      double val = myarg->get();
      for(unsigned i=0; i<ntasks; ++i) {
        inputdata[(nargs-argstart)*i + j-argstart] = val;
      }
    } else {
      for(unsigned i=0; i<ntasks; ++i) {
        inputdata[(nargs-argstart)*i + j-argstart] = myarg->get(i);
      }
    }
  }
}

template <class T>
void FunctionOfVector<T>::calculate() {
  // This is done if we are calculating a function of multiple cvs
  taskmanager.runAllTasks();
}

template <class T>
void FunctionOfVector<T>::performTask( std::size_t task_index,
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
void FunctionOfVector<T>::applyNonZeroRankForces( std::vector<double>& outforces ) {
  taskmanager.applyForces( outforces );
}

template <class T>
int FunctionOfVector<T>::getNumberOfValuesPerTask( std::size_t task_index, const FunctionData<T>& actiondata ) {
  return 1;
}

template <class T>
void FunctionOfVector<T>::getForceIndices( std::size_t task_index,
    std::size_t colno,
    std::size_t ntotal_force,
    const FunctionData<T>& actiondata,
    const ParallelActionsInput& input,
    ForceIndexHolder force_indices ) {

  unsigned vector_end = actiondata.argstart + input.nderivatives_per_scalar - actiondata.nscalars;
  for(unsigned j=0; j<input.ncomponents; ++j) {
    for(unsigned k=actiondata.argstart; k<vector_end; ++k) {
      force_indices.indices[j][k-actiondata.argstart] = input.argstarts[k] + task_index;
    }
    for(unsigned k=vector_end; k<vector_end+actiondata.nscalars; ++k) {
      force_indices.indices[j][k-actiondata.argstart] = input.argstarts[k];
    }
    force_indices.threadsafe_derivatives_end[j] = input.nderivatives_per_scalar-actiondata.nscalars;
    force_indices.tot_indices[j] = input.nderivatives_per_scalar;
  }
}

}
}
#endif
