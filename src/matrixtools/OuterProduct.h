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
#ifndef __PLUMED_matrxtools_OuterProduct_h
#define __PLUMED_matrxtools_OuterProduct_h

#include "core/ActionWithMatrix.h"
#include "core/ParallelTaskManager.h"

namespace PLMD {
namespace matrixtools {

template <class T>
class OuterProductInput {
public:
  T funcinput;
  MatrixView outmat;
};

template <class T>
class OuterProductBase : public ActionWithMatrix {
public:
  using input_type = OuterProductInput<T>;
  using PTM = ParallelTaskManager<OuterProductBase<T>>;
private:
  PTM taskmanager;
public:
  static void registerKeywords( Keywords& keys );
  explicit OuterProductBase(const ActionOptions&);
  unsigned getNumberOfDerivatives();
  void prepare() override ;
  void calculate() override ;
  void applyNonZeroRankForces( std::vector<double>& outforces ) override ;
  static void performTask( std::size_t task_index,
                           const OuterProductInput<T>& actiondata,
                           ParallelActionsInput& input,
                           ParallelActionsOutput& output );
  static void gatherForces( std::size_t task_index,
                            const OuterProductInput<T>& actiondata,
                            const ParallelActionsInput& input,
                            const ForceInput& fdata,
                            ForceOutput forces );
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
  taskmanager.setupParallelTaskManager( 1, 2, getPntrToArgument(nvals)->getNumberOfStoredValues() );
  taskmanager.setActionInput( actiondata );
}

template <class T>
unsigned OuterProductBase<T>::getNumberOfDerivatives() {
  return getPntrToArgument(0)->getNumberOfStoredValues() + getPntrToArgument(getNumberOfComponents())->getNumberOfStoredValues();
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
void OuterProductBase<T>::calculate() {
  updateBookeepingArrays( taskmanager.getActionInput().outmat );
  OuterProductInput<T>& myinp = taskmanager.getActionInput();
  unsigned ncols = getPntrToComponent(0)->getNumberOfColumns();
  taskmanager.setupParallelTaskManager( 1, 2*getNumberOfComponents()*ncols, getNumberOfComponents()*getPntrToArgument(getNumberOfComponents())->getNumberOfStoredValues(), getNumberOfComponents()*ncols );
  taskmanager.runAllTasks();
}

template <class T>
void OuterProductBase<T>::performTask( std::size_t task_index,
                                       const OuterProductInput<T>& actiondata,
                                       ParallelActionsInput& input,
                                       ParallelActionsOutput& output ) {
  std::vector<double> args(2*input.ncomponents);
  for(unsigned i=0; i<input.ncomponents; ++i) {
    args[i] = input.inputdata[input.args[i].start + task_index];
  }
  unsigned fstart = task_index*(1+actiondata.outmat.ncols);
  unsigned nelements = actiondata.outmat.bookeeping[fstart];
  for(unsigned i=0; i<nelements; ++i) {
    std::size_t argpos = actiondata.outmat.bookeeping[fstart+1+i];
    for(unsigned j=0; j<input.ncomponents; ++j) {
      args[input.ncomponents+j] = input.inputdata[input.args[input.ncomponents+j].start + argpos];
    }
    MatrixElementOutput matout( input.ncomponents, 2*input.ncomponents, output.values.data()+i*input.ncomponents, output.derivatives.data() + 2*i*input.ncomponents*input.ncomponents );
    T::calculate( input.noderiv, actiondata.funcinput, args, matout );
  }
}

template <class T>
void OuterProductBase<T>::applyNonZeroRankForces( std::vector<double>& outforces ) {
  taskmanager.applyForces( outforces );
}

template <class T>
void OuterProductBase<T>::gatherForces( std::size_t task_index,
                                        const OuterProductInput<T>& actiondata,
                                        const ParallelActionsInput& input,
                                        const ForceInput& fdata,
                                        ForceOutput forces ) {
  unsigned fstart = task_index*(1+actiondata.outmat.ncols);
  unsigned nelements = actiondata.outmat.bookeeping[fstart];
  for(unsigned i=0; i<nelements; ++i) {
    std::size_t argpos = actiondata.outmat.bookeeping[fstart+1+i];
    View2D<const double,helpers::dynamic_extent,helpers::dynamic_extent> derivs( fdata.deriv.data() + 2*i*input.ncomponents*input.ncomponents, input.ncomponents, 2*input.ncomponents );
    for(unsigned j=0; j<input.ncomponents; ++j) {
      double force = fdata.force[i*input.ncomponents+j];
      for(unsigned k=0; k<input.ncomponents; ++k) {
        forces.thread_unsafe[input.args[k].start + task_index] += force*derivs[j][k];
        forces.thread_safe[k*input.args[input.ncomponents].shape[0] + argpos] += force*derivs[j][input.ncomponents+k];
      }
    }
  }
}

}
}
#endif
