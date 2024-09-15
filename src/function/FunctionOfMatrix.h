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
#include "FunctionOfVector.h"
#include "Sum.h"
#include "tools/Matrix.h"

namespace PLMD {
namespace function {

template <class T>
class FunctionOfMatrix : public ActionWithVector {
private:
/// The function that is being computed
  T myfunc;
/// Used to hold the list of tasks we are running
  std::vector<unsigned> active_tasks;
/// Get the shape of the output matrix
  std::vector<unsigned> getValueShapeFromArguments();
/// Get a pointer to the first matrix argument
  Value* getPntrToFirstMatrixArgument() const ;
public:
  static void registerKeywords(Keywords&);
  explicit FunctionOfMatrix(const ActionOptions&);
/// Get the label to write in the graph
  std::string writeInGraph() const override { return myfunc.getGraphInfo( getName() ); }
/// Make sure the derivatives are turned on
  void turnOnDerivatives() override;
/// Get the number of derivatives for this action
  unsigned getNumberOfDerivatives() override ;
/// Resize the matrices
  void prepare() override ;
  void calculate() override ;
  std::vector<unsigned>& getListOfActiveTasks( ActionWithVector* action ) override ;
  void performTask( const unsigned& current, MultiValue& myvals ) const override ;
};

template <class T>
void FunctionOfMatrix<T>::registerKeywords(Keywords& keys ) {
  ActionWithVector::registerKeywords(keys); keys.use("ARG"); std::string name = keys.getDisplayName();
  std::size_t und=name.find("_MATRIX"); keys.setDisplayName( name.substr(0,und) ); keys.use("MASK");
  keys.add("hidden","NO_ACTION_LOG","suppresses printing from action on the log");
  keys.reserve("compulsory","PERIODIC","if the output of your function is periodic then you should specify the periodicity of the function.  If the output is not periodic you must state this using PERIODIC=NO");
  T tfunc; tfunc.registerKeywords( keys );
  if( keys.getDisplayName()=="SUM" ) {
    keys.setValueDescription("the sum of all the elements in the input matrix");
  } else if( keys.getDisplayName()=="HIGHEST" ) {
    keys.setValueDescription("the largest element of the input matrix");
  } else if( keys.getDisplayName()=="LOWEST" ) {
    keys.setValueDescription("the smallest element in the input matrix");
  } else if( keys.outputComponentExists(".#!value") ) {
    keys.setValueDescription("the matrix obtained by doing an element-wise application of " + keys.getOutputComponentDescription(".#!value") + " to the input matrix");
  }
}

template <class T>
FunctionOfMatrix<T>::FunctionOfMatrix(const ActionOptions&ao):
  Action(ao),
  ActionWithVector(ao)
{
  // Get the shape of the output
  std::vector<unsigned> shape( getValueShapeFromArguments() );
  // Check if the output matrix is symmetric
  bool symmetric=true; unsigned argstart=myfunc.getArgStart();
  for(unsigned i=argstart; i<getNumberOfArguments(); ++i) {
    if( getPntrToArgument(i)->getRank()==2 ) {
      if( !getPntrToArgument(i)->isSymmetric() ) { symmetric=false;  }
    }
  }
  unsigned nargs = getNumberOfArguments(); if( getNumberOfMasks()>0 ) nargs = nargs - getNumberOfMasks();
  // Read the input and do some checks
  myfunc.read( this );
  // Check we are not calculating a sum
  if( myfunc.zeroRank() ) shape.resize(0);
  // Get the names of the components
  std::vector<std::string> components( keywords.getOutputComponents() );
  // Create the values to hold the output
  std::vector<std::string> str_ind( myfunc.getComponentsPerLabel() );
  for(unsigned i=0; i<components.size(); ++i) {
    if( str_ind.size()>0 ) {
      std::string compstr = components[i]; if( components[i]==".#!value" ) compstr = "";
      for(unsigned j=0; j<str_ind.size(); ++j) {
        if( myfunc.zeroRank() ) {
          addComponentWithDerivatives( compstr + str_ind[j], shape );
        } else {
          addComponent( compstr + str_ind[j], shape );
          getPntrToComponent(i*str_ind.size()+j)->setSymmetric( symmetric );
        }
      }
    } else if( components[i]==".#!value" && myfunc.zeroRank() ) {
      addValueWithDerivatives( shape );
    } else if( components[i]==".#!value" ) {
      addValue( shape ); getPntrToComponent(0)->setSymmetric( symmetric );
    } else if( components[i].find_first_of("_")!=std::string::npos ) {
      if( nargs-argstart==1 ) { addValue( shape ); getPntrToComponent(0)->setSymmetric( symmetric ); }
      else {
        for(unsigned j=argstart; j<nargs; ++j) {
          addComponent( getPntrToArgument(j)->getName() + components[i], shape );
          getPntrToComponent(i*(nargs-argstart)+j-argstart)->setSymmetric( symmetric );
        }
      }
    } else { addComponent( components[i], shape ); getPntrToComponent(i)->setSymmetric( symmetric ); }
  }
  // Check if this can be sped up
  if( myfunc.getDerivativeZeroIfValueIsZero() )  {
    for(int i=0; i<getNumberOfComponents(); ++i) getPntrToComponent(i)->setDerivativeIsZeroWhenValueIsZero();
  }
  // Set the periodicities of the output components
  myfunc.setPeriodicityForOutputs( this );
  // Now setup the action in the chain if we can
  unsigned nderivatives = buildArgumentStore(myfunc.getArgStart());
}

template <class T>
void FunctionOfMatrix<T>::turnOnDerivatives() {
  if( !myfunc.derivativesImplemented() ) error("derivatives have not been implemended for " + getName() );
  ActionWithValue::turnOnDerivatives(); myfunc.setup(this);
}

template <class T>
std::vector<unsigned> FunctionOfMatrix<T>::getValueShapeFromArguments() {
  unsigned argstart=myfunc.getArgStart(); std::vector<unsigned> shape(2); shape[0]=shape[1]=0;
  unsigned nargs = getNumberOfArguments(); if( getNumberOfMasks()>0 ) nargs = nargs - getNumberOfMasks();
  for(unsigned i=argstart; i<nargs; ++i) {
    plumed_assert( getPntrToArgument(i)->getRank()==2 || getPntrToArgument(i)->getRank()==0 );
    if( getPntrToArgument(i)->getRank()==2 ) {
      if( shape[0]>0 && (getPntrToArgument(i)->getShape()[0]!=shape[0] || getPntrToArgument(i)->getShape()[1]!=shape[1]) ) error("all matrices input should have the same shape");
      else if( shape[0]==0 ) { shape[0]=getPntrToArgument(i)->getShape()[0]; shape[1]=getPntrToArgument(i)->getShape()[1]; }
      plumed_assert( !getPntrToArgument(i)->hasDerivatives() );
    }
  }
  myfunc.setPrefactor( this, 1.0 ); return shape;
}


template <class T>
unsigned FunctionOfMatrix<T>::getNumberOfDerivatives() {
  unsigned nder=0, argstart = myfunc.getArgStart();
  unsigned nargs = getNumberOfArguments(); if( getNumberOfMasks()>0 ) nargs = nargs - getNumberOfMasks();
  for(unsigned i=argstart; i<nargs; ++i) nder += getPntrToArgument(i)->getNumberOfStoredValues();
  return nder;
}

template <class T>
void FunctionOfMatrix<T>::prepare() {
  unsigned argstart = myfunc.getArgStart(); std::vector<unsigned> shape(2);
  for(unsigned i=argstart; i<getNumberOfArguments(); ++i) {
    if( getPntrToArgument(i)->getRank()==2 ) {
      shape[0] = getPntrToArgument(i)->getShape()[0];
      shape[1] = getPntrToArgument(i)->getShape()[1];
      break;
    }
  }
  for(unsigned i=0; i<getNumberOfComponents(); ++i) {
    Value* myval = getPntrToComponent(i);
    if( myval->getRank()==2 && (myval->getShape()[0]!=shape[0] || myval->getShape()[1]!=shape[1]) ) myval->setShape(shape);
  }
  ActionWithVector::prepare(); active_tasks.resize(0);
}

template <class T>
Value* FunctionOfMatrix<T>::getPntrToFirstMatrixArgument() const {
  unsigned argstart=myfunc.getArgStart();
  for(unsigned i=argstart; i<getNumberOfArguments(); ++i) {
    if( getPntrToArgument(i)->getRank()==2 ) return getPntrToArgument(i);
  }
  plumed_merror("could not find matrix argument");
  return NULL;
}

template <class T>
std::vector<unsigned>& FunctionOfMatrix<T>::getListOfActiveTasks( ActionWithVector* action ) {
  if( active_tasks.size()>0 ) return active_tasks;

  Value* myarg = NULL;
  if( getNumberOfMasks()>0 ) myarg = getPntrToArgument(getNumberOfArguments()-getNumberOfMasks());
  else myarg = getPntrToFirstMatrixArgument();

  unsigned base=0; unsigned nrows = myarg->getShape()[0];
  for(unsigned i=0; i<nrows; ++i) {
    unsigned ncols = myarg->getRowLength(i);
    for(unsigned j=0; j<ncols; ++j) active_tasks.push_back(base+j);
    base += myarg->getNumberOfColumns();
  }
  if( getNumberOfMasks()>0 && doNotCalculateDerivatives() ) return active_tasks;
// All the matrices input have to have the same sparsity pattern.
// I can do everything I want to do with this limitation.  If
// anyone wants to make this smarter in the future they can
#ifndef DNDEBUG
  unsigned argstart=myfunc.getArgStart();
  for(unsigned k=argstart; k<getNumberOfArguments(); ++k) {
    if( getPntrToArgument(k)->getRank()!=2 ) continue ;
    for(unsigned i=0; i<nrows; ++i) {
      unsigned ncols = myarg->getRowLength(i);
      if( getNumberOfMasks()>0 && ncols==0 ) continue;
      plumed_massert( ncols==getPntrToArgument(k)->getRowLength(i), "failing in " + getLabel() );
      for(unsigned j=0; j<ncols; ++j) plumed_assert( myarg->getRowIndex(i,j)==getPntrToArgument(k)->getRowIndex(i,j) );
    }
  }
#endif
  return active_tasks;
}

template <class T>
void FunctionOfMatrix<T>::performTask( const unsigned& taskno, MultiValue& myvals) const {
  unsigned nargs = getNumberOfArguments(); if( getNumberOfMasks()>0 ) nargs = nargs - getNumberOfMasks();
  unsigned argstart=myfunc.getArgStart(); std::vector<double> args( nargs - argstart );

  // Retrieve the arguments
  if( getNumberOfMasks()>0 ) {
    Value* maskarg = getPntrToArgument(getNumberOfArguments()-getNumberOfMasks());
    unsigned index1 = std::floor( taskno / maskarg->getNumberOfColumns() );
    unsigned index2 = taskno - index1*maskarg->getNumberOfColumns();
    unsigned tind2 = maskarg->getRowIndex( index1, index2 );
    unsigned maskncol = maskarg->getRowLength(index1);
    for(unsigned i=argstart; i<nargs; ++i) {
      if( getPntrToArgument(i)->getRank()==2 && getPntrToArgument(i)->getRowLength(index1)>maskncol ) args[i-argstart]=getPntrToArgument(i)->get( index1*getPntrToArgument(i)->getShape()[1] + tind2 );
      else if( getPntrToArgument(i)->getRank()==2 ) {
        plumed_dbg_assert( maskncol==getPntrToArgument(i)->getRowLength(index1) );
        args[i-argstart]=getPntrToArgument(i)->get( taskno, false );
      } else args[i-argstart] = getPntrToArgument(i)->get();
    }
  } else {
    for(unsigned i=argstart; i<nargs; ++i) {
      if( getPntrToArgument(i)->getRank()==2 ) args[i-argstart]=getPntrToArgument(i)->get( taskno, false );
      else args[i-argstart] = getPntrToArgument(i)->get();
    }
  }
  // Calculate the function and its derivatives
  std::vector<double> vals( getNumberOfComponents() ); Matrix<double> derivatives( getNumberOfComponents(), nargs-argstart );
  myfunc.calc( this, args, vals, derivatives );

  // And set the values
  for(unsigned i=0; i<vals.size(); ++i) myvals.addValue( i, vals[i] );

  // Return if we are not computing derivatives
  if( doNotCalculateDerivatives() ) return;

  unsigned base=0;
  for(unsigned j=argstart; j<nargs; ++j) {
    if( getPntrToArgument(j)->getRank()==2 ) {
      for(int i=0; i<getNumberOfComponents(); ++i) {
        myvals.addDerivative( i, base + taskno, derivatives(i,j) ); myvals.updateIndex( i, base + taskno );
      }
    } else {
      for(int i=0; i<getNumberOfComponents(); ++i) {
        myvals.addDerivative( i, base, derivatives(i,j) );
        myvals.updateIndex( i, base );
      }
    }
    base += getPntrToArgument(j)->getNumberOfStoredValues();
  }
}

template <class T>
void FunctionOfMatrix<T>::calculate() {
  Value* myarg = NULL;
  if( getNumberOfMasks()>0 ) myarg = getPntrToArgument(getNumberOfArguments()-getNumberOfMasks());
  else myarg = getPntrToFirstMatrixArgument();
  // Copy bookeeping arrays from input matrices to output matrices
  for(unsigned i=0; i<getNumberOfComponents(); ++i) {
    if( getPntrToComponent(i)->getRank()==2 ) getPntrToComponent(i)->copyBookeepingArrayFromArgument( myarg );
    else getPntrToComponent(i)->resizeDerivatives( getNumberOfDerivatives() );
  }
  runAllTasks();
}

}
}
#endif
