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
//#include "core/CollectFrames.h"
#include "core/ActionSetup.h"
#include "tools/Matrix.h"
#include "Sum.h"

namespace PLMD {
namespace function {

template <class T>
class FunctionOfVector : public ActionWithVector {
private:
/// Do the calculation at the end of the run
  bool doAtEnd;
/// The function that is being computed
  T myfunc;
public:
  static void registerKeywords(Keywords&);
/// This method is used to run the calculation with functions such as highest/lowest and sort.
/// It is static so we can reuse the functionality in FunctionOfMatrix
  static void runSingleTaskCalculation( const Value* arg, ActionWithValue* action, T& f );
  explicit FunctionOfVector(const ActionOptions&);
  ~FunctionOfVector() {}
  std::string getOutputComponentDescription( const std::string& cname, const Keywords& keys ) const override ;
/// Get the size of the task list at the end of the run
  unsigned getNumberOfFinalTasks();
/// Check if derivatives are available
  void turnOnDerivatives() override;
/// Get the number of derivatives for this action
  unsigned getNumberOfDerivatives() override ;
/// Resize vectors that are the wrong size
  void prepare() override ;
/// Get the label to write in the graph
  std::string writeInGraph() const override {
    return myfunc.getGraphInfo( getName() );
  }
/// This builds the task list for the action
  void calculate() override;
/// Calculate the function
  void performTask( const unsigned& current, MultiValue& myvals ) const override ;
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
  T tfunc;
  tfunc.registerKeywords( keys );
  if( keys.getDisplayName()=="SUM" ) {
    keys.setValueDescription("scalar","the sum of all the elements in the input vector");
  } else if( keys.getDisplayName()=="MEAN" ) {
    keys.setValueDescription("scalar","the mean of all the elements in the input vector");
  } else if( keys.getDisplayName()=="HIGHEST" ) {
    keys.setValueDescription("scalar/vector","the largest element of the input vector if one vector specified.  If multiple vectors of the same size specified the largest elements of these vector computed elementwise.");
  } else if( keys.getDisplayName()=="LOWEST" ) {
    keys.setValueDescription("scalar/vector","the smallest element in the input vector if one vector specified.  If multiple vectors of the same size specified the largest elements of these vector computed elementwise.");
  } else if( keys.getDisplayName()=="SORT" ) {
    keys.setValueDescription("vector","a vector that has been sorted into ascending order");
  } else if( keys.outputComponentExists(".#!value") ) {
    keys.add("optional","MASK","the label for a sparse matrix that should be used to determine which elements of the matrix should be computed");
    keys.setValueDescription("vector","the vector obtained by doing an element-wise application of " + keys.getOutputComponentDescription(".#!value") + " to the input vectors");
  }
}

template <class T>
FunctionOfVector<T>::FunctionOfVector(const ActionOptions&ao):
  Action(ao),
  ActionWithVector(ao),
  doAtEnd(true) {
  // Get the shape of the output
  std::vector<std::size_t> shape(1);
  shape[0]=getNumberOfFinalTasks();
  // Read the input and do some checks
  myfunc.read( this );
  // Create the task list
  if( myfunc.doWithTasks() ) {
    doAtEnd=false;
  } else if( getNumberOfArguments()!=1 ) {
    error("number of arguments should be equal to one");
  }
  // Get the names of the components
  std::vector<std::string> components( keywords.getOutputComponents() );
  // Create the values to hold the output
  std::vector<std::string> str_ind( myfunc.getComponentsPerLabel() );
  for(unsigned i=0; i<components.size(); ++i) {
    if( str_ind.size()>0 ) {
      std::string strcompn = components[i];
      if( components[i]==".#!value" ) {
        strcompn = "";
      }
      for(unsigned j=0; j<str_ind.size(); ++j) {
        if( myfunc.zeroRank() ) {
          addComponentWithDerivatives( strcompn + str_ind[j] );
        } else {
          addComponent( strcompn + str_ind[j], shape );
        }
      }
    } else if( components[i].find_first_of("_")!=std::string::npos ) {
      if( getNumberOfArguments()==1 && myfunc.zeroRank() ) {
        addValueWithDerivatives();
      } else if( getNumberOfArguments()==1 ) {
        addValue( shape );
      } else {
        unsigned argstart=myfunc.getArgStart();
        for(unsigned i=argstart; i<getNumberOfArguments(); ++i) {
          if( myfunc.zeroRank() ) {
            addComponentWithDerivatives( getPntrToArgument(i)->getName() + components[i] );
          } else {
            addComponent( getPntrToArgument(i)->getName() + components[i], shape );
          }
        }
      }
    } else if( components[i]==".#!value" && myfunc.zeroRank() ) {
      addValueWithDerivatives();
    } else if( components[i]==".#!value" ) {
      addValue(shape);
    } else if( myfunc.zeroRank() ) {
      addComponentWithDerivatives( components[i] );
    } else {
      addComponent( components[i], shape );
    }
  }
  // Check if we can turn off the derivatives when they are zero
  if( myfunc.getDerivativeZeroIfValueIsZero() )  {
    for(int i=0; i<getNumberOfComponents(); ++i) {
      getPntrToComponent(i)->setDerivativeIsZeroWhenValueIsZero();
    }
  }
  // Check if this is a timeseries
  unsigned argstart=myfunc.getArgStart();
  // Set the periodicities of the output components
  myfunc.setPeriodicityForOutputs( this );
}

template <class T>
std::string FunctionOfVector<T>::getOutputComponentDescription( const std::string& cname, const Keywords& keys ) const {
  if( getName().find("SORT")==std::string::npos ) {
    return ActionWithValue::getOutputComponentDescription( cname, keys );
  }
  if( getNumberOfArguments()==1 ) {
    return "the " + cname + "th largest element of the vector " + getPntrToArgument(0)->getName();
  }
  return "the " + cname + "th largest element in the input vectors";
}

template <class T>
void FunctionOfVector<T>::turnOnDerivatives() {
  if( !getPntrToComponent(0)->isConstant() && !myfunc.derivativesImplemented() ) {
    error("derivatives have not been implemended for " + getName() );
  }
  ActionWithValue::turnOnDerivatives();
  myfunc.setup(this );
}

template <class T>
unsigned FunctionOfVector<T>::getNumberOfDerivatives() {
  unsigned nder = 0, argstart = myfunc.getArgStart();
  for(unsigned i=argstart; i<getNumberOfArguments(); ++i) {
    nder += getPntrToArgument(i)->getNumberOfStoredValues();
  }
  return nder;
}

template <class T>
void FunctionOfVector<T>::prepare() {
  unsigned argstart = myfunc.getArgStart();
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
void FunctionOfVector<T>::performTask( const unsigned& current, MultiValue& myvals ) const {
  unsigned nargs=getNumberOfArguments();
  if( getNumberOfMasks()>0 ) {
    nargs = nargs - getNumberOfMasks();
  }
  unsigned argstart=myfunc.getArgStart();
  std::vector<double> args( nargs-argstart);
  for(unsigned i=argstart; i<nargs; ++i) {
    if( getPntrToArgument(i)->getRank()==1 ) {
      args[i-argstart]=getPntrToArgument(i)->get(current);
    } else {
      args[i-argstart] = getPntrToArgument(i)->get();
    }
  }
  // Calculate the function and its derivatives
  std::vector<double> vals( getNumberOfComponents() );
  Matrix<double> derivatives( getNumberOfComponents(), args.size() );
  myfunc.calc( this, args, vals, derivatives );
  // And set the values
  for(unsigned i=0; i<vals.size(); ++i) {
    myvals.addValue( i, vals[i] );
  }
  // Return if we are not computing derivatives
  if( doNotCalculateDerivatives() ) {
    return;
  }

  unsigned base=0;
  for(unsigned j=0; j<args.size(); ++j) {
    if( getPntrToArgument(argstart+j)->getRank()==1 ) {
      for(int i=0; i<getNumberOfComponents(); ++i) {
        myvals.addDerivative( i, base+current, derivatives(i,j) );
        myvals.updateIndex( i, base+current );
      }
    } else {
      for(int i=0; i<getNumberOfComponents(); ++i) {
        myvals.addDerivative( i, base, derivatives(i,j) );
        myvals.updateIndex( i, base );
      }
    }
    base += getPntrToArgument(argstart+j)->getNumberOfValues();
  }
}

template <class T>
unsigned FunctionOfVector<T>::getNumberOfFinalTasks() {
  unsigned nelements=0, argstart=myfunc.getArgStart();
  for(unsigned i=argstart; i<getNumberOfArguments(); ++i) {
    plumed_assert( getPntrToArgument(i)->getRank()<2 );
    if( getPntrToArgument(i)->getRank()==1 ) {
      if( nelements>0 ) {
        if(getPntrToArgument(i)->getShape()[0]!=nelements ) {
          error("all vectors input should have the same length");
        }
      } else if( nelements==0 ) {
        nelements=getPntrToArgument(i)->getShape()[0];
      }
      plumed_assert( !getPntrToArgument(i)->hasDerivatives() );
    }
  }
  // The prefactor for average and sum is set here so the number of input scalars is guaranteed to be correct
  myfunc.setPrefactor( this, 1.0 );
  return nelements;
}

template <class T>
void FunctionOfVector<T>::runSingleTaskCalculation( const Value* arg, ActionWithValue* action, T& f ) {
  // This is used if we are doing sorting actions on a single vector
  unsigned nv = arg->getNumberOfValues();
  std::vector<double> args( nv );
  for(unsigned i=0; i<nv; ++i) {
    args[i] = arg->get(i);
  }
  std::vector<double> vals( action->getNumberOfComponents() );
  Matrix<double> derivatives( action->getNumberOfComponents(), nv );
  ActionWithArguments* aa=dynamic_cast<ActionWithArguments*>(action);
  plumed_assert( aa );
  f.calc( aa, args, vals, derivatives );
  for(unsigned i=0; i<vals.size(); ++i) {
    action->copyOutput(i)->set( vals[i] );
  }
  // Return if we are not computing derivatives
  if( action->doNotCalculateDerivatives() ) {
    return;
  }
  // Now set the derivatives
  for(unsigned j=0; j<nv; ++j) {
    for(unsigned i=0; i<vals.size(); ++i) {
      action->copyOutput(i)->setDerivative( j, derivatives(i,j) );
    }
  }
}

template <class T>
void FunctionOfVector<T>::calculate() {
  // This is done if we are calculating a function of multiple cvs
  if( !doAtEnd ) {
    runAllTasks();
  }
  // This is used if we are doing sorting actions on a single vector
  else if( !myfunc.doWithTasks() ) {
    runSingleTaskCalculation( getPntrToArgument(0), this, myfunc );
  }
}

}
}
#endif
