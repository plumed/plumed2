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
#ifndef __PLUMED_matrixtools_OuterProduct_h
#define __PLUMED_matrixtools_OuterProduct_h

#include "core/ActionWithMatrix.h"
#include "core/ParallelTaskManager.h"

namespace PLMD {
namespace matrixtools {

template <class T>
class OuterProductInput {
public:
  T funcinput;
  RequiredMatrixElements outmat;
};

template <class T>
class OuterProductBase : public ActionWithMatrix {
public:
  using input_type = OuterProductInput<T>;
  using PTM = ParallelTaskManager<OuterProductBase<T>>;
private:
  bool isproduct;
  PTM taskmanager;
public:
  static void registerKeywords( Keywords& keys );
  explicit OuterProductBase(const ActionOptions&);
  unsigned getNumberOfDerivatives() override;
  void prepare() override ;
  int checkTaskIsActive( const unsigned& itask ) const override ;
  void calculate() override ;
  void applyNonZeroRankForces( std::vector<double>& outforces ) override ;
  static void performTask( std::size_t task_index,
                           const OuterProductInput<T>& actiondata,
                           ParallelActionsInput& input,
                           ParallelActionsOutput& output );
  static int getNumberOfValuesPerTask( std::size_t task_index,
                                       const OuterProductInput<T>& actiondata );
  static void getForceIndices( std::size_t task_index,
                               std::size_t colno,
                               std::size_t ntotal_force,
                               const OuterProductInput<T>& actiondata,
                               const ParallelActionsInput& input,
                               ForceIndexHolder force_indices );
};

template <class T>
void OuterProductBase<T>::registerKeywords( Keywords& keys ) {
  ActionWithMatrix::registerKeywords(keys);
  T::registerKeywords( keys );
  PTM::registerKeywords( keys );
  keys.add("hidden","NO_ACTION_LOG","suppresses printing from action on the log");
}

template <class T>
OuterProductBase<T>::OuterProductBase(const ActionOptions&ao):
  Action(ao),
  ActionWithMatrix(ao),
  isproduct(false),
  taskmanager(this) {
  unsigned nargs=getNumberOfArguments();
  if( getNumberOfMasks()>0 ) {
    nargs = nargs - getNumberOfMasks();
  }
  if( nargs%2!=0 ) {
    error("should be an even number of arguments to this action, they should all be vectors");
  }
  std::size_t nvals = nargs / 2;
  for(unsigned i=0; i<nvals; ++i) {
    if( getPntrToArgument(i)->getRank()!=1 || getPntrToArgument(i)->hasDerivatives() ) {
      error("first argument to this action should be a vector");
    }
    if( getPntrToArgument(0)->getShape()[0]!=getPntrToArgument(i)->getShape()[0] ) {
      error("mismatch between sizes of input vectors");
    }
    if( getPntrToArgument(nvals+i)->getRank()!=1 || getPntrToArgument(nvals+i)->hasDerivatives() ) {
      error("first argument to this action should be a vector");
    }
    if( getPntrToArgument(nvals)->getShape()[0]!=getPntrToArgument(nvals+i)->getShape()[0] ) {
      error("mismatch between sizes of input vectors");
    }
  }
  if( getNumberOfMasks()==1 ) {
    if( getPntrToArgument(nargs)->getRank()!=2 || getPntrToArgument(nargs)->hasDerivatives() ) {
      error("mask argument should be a matrix");
    }
    for(unsigned i=0; i<nvals; ++i) {
      if( getPntrToArgument(nargs)->getShape()[0]!=getPntrToArgument(i)->getShape()[0] ) {
        error("mask argument has wrong size");
      }
      if( getPntrToArgument(nargs)->getShape()[1]!=getPntrToArgument(nvals+i)->getShape()[0] ) {
        error("mask argument has wrong size");
      }
    }
  }

  std::vector<std::size_t> shape(2);
  shape[0]=getPntrToArgument(0)->getShape()[0];
  shape[1]=getPntrToArgument(nvals)->getShape()[0];

  std::string func;
  if( keywords.exists("FUNC") ) {
    parse("FUNC",func);
    isproduct=(func=="x*y");
  }
  OuterProductInput<T> actiondata;
  actiondata.funcinput.setup( shape, func, this );

  if( getNumberOfComponents()==0 ) {
    addValue( shape );
    setNotPeriodic();
  }
  if( getPntrToArgument(0)->isDerivativeZeroWhenValueIsZero() || getPntrToArgument(nvals)->isDerivativeZeroWhenValueIsZero() ) {
    for(unsigned i=0; i<getNumberOfComponents(); ++i) {
      getPntrToComponent(i)->setDerivativeIsZeroWhenValueIsZero();
    }
  }
  taskmanager.setActionInput( actiondata );
}

template <class T>
unsigned OuterProductBase<T>::getNumberOfDerivatives() {
  unsigned nc = getNumberOfComponents();
  return nc*(getPntrToArgument(0)->getNumberOfStoredValues() + getPntrToArgument(nc)->getNumberOfStoredValues());
}

template <class T>
void OuterProductBase<T>::prepare() {
  ActionWithVector::prepare();
  Value* myval=getPntrToComponent(0);
  if( myval->getShape()[0]==getPntrToArgument(0)->getShape()[0] && myval->getShape()[1]==getPntrToArgument(getNumberOfComponents())->getShape()[0] ) {
    return;
  }
  std::vector<std::size_t> shape(2);
  shape[0] = getPntrToArgument(0)->getShape()[0];
  shape[1] = getPntrToArgument(getNumberOfComponents())->getShape()[0];
  for(unsigned i=0; i<getNumberOfComponents(); ++i) {
    getPntrToComponent(i)->setShape( shape );
  }
}

template <class T>
int OuterProductBase<T>::checkTaskIsActive( const unsigned& itask ) const {
  if( getNumberOfMasks()>0 || !isproduct ) {
    return ActionWithVector::checkTaskIsActive( itask );
  }
  if( fabs( getPntrToArgument(0)->get(itask))>epsilon ) {
    return 1;
  }
  return -1;
}

template <class T>
void OuterProductBase<T>::calculate() {
  updateBookeepingArrays( taskmanager.getActionInput().outmat );
  taskmanager.setupParallelTaskManager( 2*getNumberOfComponents(), getNumberOfComponents()*getPntrToComponent(0)->getShape()[1] );
  taskmanager.setWorkspaceSize( 2*getNumberOfComponents() );
  taskmanager.runAllTasks();
}

template <class T>
void OuterProductBase<T>::performTask( std::size_t task_index,
                                       const OuterProductInput<T>& actiondata,
                                       ParallelActionsInput& input,
                                       ParallelActionsOutput& output ) {
  auto args = output.buffer.subview(0, 2*input.ncomponents);
  for(unsigned i=0; i<input.ncomponents; ++i) {
    args[i] = input.inputdata[input.argstarts[i] + task_index];
  }
  unsigned fstart = task_index*(1+actiondata.outmat.ncols);
  unsigned nelements = actiondata.outmat[fstart];
  for(unsigned i=0; i<nelements; ++i) {
    std::size_t argpos = actiondata.outmat[fstart+1+i];
    for(unsigned j=0; j<input.ncomponents; ++j) {
      args[input.ncomponents+j] = input.inputdata[input.argstarts[input.ncomponents+j] + argpos];
    }
    MatrixElementOutput matout( input.ncomponents,
                                2*input.ncomponents,
                                output.values.data()+i*input.ncomponents,
                                output.derivatives.data() + 2*i*input.ncomponents*input.ncomponents );
    T::calculate( input.noderiv, actiondata.funcinput, {args.data(),args.size()}, matout );
  }
}

template <class T>
void OuterProductBase<T>::applyNonZeroRankForces( std::vector<double>& outforces ) {
  taskmanager.applyForces( outforces );
}

template <class T>
int OuterProductBase<T>::getNumberOfValuesPerTask( std::size_t task_index,
    const OuterProductInput<T>& actiondata ) {
  unsigned fstart = task_index*(1+actiondata.outmat.ncols);
  return actiondata.outmat[fstart];
}

template <class T>
void OuterProductBase<T>::getForceIndices( std::size_t task_index,
    std::size_t colno,
    std::size_t ntotal_force,
    const OuterProductInput<T>& actiondata,
    const ParallelActionsInput& input,
    ForceIndexHolder force_indices ) {
  unsigned fstart = task_index*(1+actiondata.outmat.ncols);
  for(unsigned j=0; j<input.ncomponents; ++j) {
    for(unsigned k=0; k<input.ncomponents; ++k) {
      force_indices.indices[j][k] = input.argstarts[k] + task_index;
      force_indices.indices[j][input.ncomponents+k] = input.argstarts[input.ncomponents+k] + actiondata.outmat[fstart+1+colno];
    }
    force_indices.threadsafe_derivatives_end[j] = input.ncomponents;
    force_indices.tot_indices[j] = 2*input.ncomponents;
  }
}

}
}
#endif
