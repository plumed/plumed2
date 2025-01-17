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
/// Is this the first time we are doing the calc
  bool firststep;
/// The function that is being computed
  T myfunc;
/// The number of derivatives for this action
  unsigned nderivatives;
/// A vector that tells us if we have stored the input value
  std::vector<bool> stored_arguments;
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
/// Check if all he actions are required
  void areAllTasksRequired( std::vector<ActionWithVector*>& task_reducing_actions );
/// Get the label to write in the graph
  std::string writeInGraph() const override {
    return myfunc.getGraphInfo( getName() );
  }
/// This builds the task list for the action
  void calculate() override;
/// This ensures that we create some bookeeping stuff during the first step
  void setupStreamedComponents( const std::string& headstr, unsigned& nquants, unsigned& nmat, unsigned& maxcol, unsigned& nbookeeping ) override ;
/// Calculate the function
  void performTask( const unsigned& current, MultiValue& myvals ) const override ;
};

template <class T>
void FunctionOfVector<T>::registerKeywords(Keywords& keys ) {
  Action::registerKeywords(keys);
  ActionWithValue::registerKeywords(keys);
  ActionWithArguments::registerKeywords(keys);
  keys.use("ARG");
  std::string name = keys.getDisplayName();
  std::size_t und=name.find("_VECTOR");
  keys.setDisplayName( name.substr(0,und) );
  keys.reserve("compulsory","PERIODIC","if the output of your function is periodic then you should specify the periodicity of the function.  If the output is not periodic you must state this using PERIODIC=NO");
  keys.add("hidden","NO_ACTION_LOG","suppresses printing from action on the log");
  T tfunc;
  tfunc.registerKeywords( keys );
  if( keys.getDisplayName()=="SUM" ) {
    keys.setValueDescription("the sum of all the elements in the input vector");
  } else if( keys.getDisplayName()=="MEAN" ) {
    keys.setValueDescription("the mean of all the elements in the input vector");
  } else if( keys.getDisplayName()=="HIGHEST" ) {
    keys.setValueDescription("the largest element of the input vector");
  } else if( keys.getDisplayName()=="LOWEST" ) {
    keys.setValueDescription("the smallest element in the input vector");
  } else if( keys.getDisplayName()=="SORT" ) {
    keys.setValueDescription("a vector that has been sorted into ascending order");
  } else if( keys.outputComponentExists(".#!value") ) {
    keys.setValueDescription("the vector obtained by doing an element-wise application of " + keys.getOutputComponentDescription(".#!value") + " to the input vectors");
  }
}

template <class T>
FunctionOfVector<T>::FunctionOfVector(const ActionOptions&ao):
  Action(ao),
  ActionWithVector(ao),
  doAtEnd(true),
  firststep(true),
  nderivatives(0) {
  // Get the shape of the output
  std::vector<unsigned> shape(1);
  shape[0]=getNumberOfFinalTasks();
  // Read the input and do some checks
  myfunc.read( this );
  // Create the task list
  if( myfunc.doWithTasks() ) {
    doAtEnd=false;
    if( shape[0]>0 ) {
      done_in_chain=true;
    }
  } else {
    plumed_assert( getNumberOfArguments()==1 );
    done_in_chain=false;
    getPntrToArgument(0)->buildDataStore();
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
  // for(unsigned i=argstart; i<getNumberOfArguments();++i) {
  //   if( getPntrToArgument(i)->isTimeSeries() ) {
  //       for(unsigned i=0; i<getNumberOfComponents(); ++i) getPntrToOutput(i)->makeHistoryDependent();
  //       break;
  //   }
  // }
  // Set the periodicities of the output components
  myfunc.setPeriodicityForOutputs( this );
  // Check if we can put the function in a chain
  for(unsigned i=argstart; i<getNumberOfArguments(); ++i) {
    // CollectFrames* ab=dynamic_cast<CollectFrames*>( getPntrToArgument(i)->getPntrToAction() );
    // if( ab && ab->hasClear() ) { doNotChain=true; getPntrToArgument(i)->buildDataStore( getLabel() ); }
    // No chains if we are using a sum or a mean
    if( getPntrToArgument(i)->getRank()==0 ) {
      FunctionOfVector<Sum>* as = dynamic_cast<FunctionOfVector<Sum>*>( getPntrToArgument(i)->getPntrToAction() );
      if(as) {
        done_in_chain=false;
      }
    } else {
      ActionWithVector* av=dynamic_cast<ActionWithVector*>( getPntrToArgument(i)->getPntrToAction() );
      if( !av ) {
        done_in_chain=false;
      }
    }
  }
  // Don't need to do the calculation in a chain if the input is constant
  bool allconstant=true;
  for(unsigned i=argstart; i<getNumberOfArguments(); ++i) {
    if( !getPntrToArgument(i)->isConstant() ) {
      allconstant=false;
      break;
    }
  }
  if( allconstant ) {
    done_in_chain=false;
  }
  nderivatives = buildArgumentStore(myfunc.getArgStart());
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
  return nderivatives;
}

template <class T>
void FunctionOfVector<T>::prepare() {
  unsigned argstart = myfunc.getArgStart();
  std::vector<unsigned> shape(1);
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
void FunctionOfVector<T>::setupStreamedComponents( const std::string& headstr, unsigned& nquants, unsigned& nmat, unsigned& maxcol, unsigned& nbookeeping ) {
  if( firststep ) {
    stored_arguments.resize( getNumberOfArguments() );
    std::string control = getFirstActionInChain()->getLabel();
    for(unsigned i=0; i<stored_arguments.size(); ++i) {
      if( getPntrToArgument(i)->isConstant() ) {
        stored_arguments[i]=false;
      } else {
        stored_arguments[i] = !getPntrToArgument(i)->ignoreStoredValue( control );
      }
    }
    firststep=false;
  }
  ActionWithVector::setupStreamedComponents( headstr, nquants, nmat, maxcol, nbookeeping );
}

template <class T>
void FunctionOfVector<T>::performTask( const unsigned& current, MultiValue& myvals ) const {
  unsigned argstart=myfunc.getArgStart();
  std::vector<double> args( getNumberOfArguments()-argstart);
  if( actionInChain() ) {
    for(unsigned i=argstart; i<getNumberOfArguments(); ++i) {
      if(  getPntrToArgument(i)->getRank()==0 ) {
        args[i-argstart] = getPntrToArgument(i)->get();
      } else if( !getPntrToArgument(i)->valueHasBeenSet() ) {
        args[i-argstart] = myvals.get( getPntrToArgument(i)->getPositionInStream() );
      } else {
        args[i-argstart] = getPntrToArgument(i)->get( myvals.getTaskIndex() );
      }
    }
  } else {
    for(unsigned i=argstart; i<getNumberOfArguments(); ++i) {
      if( getPntrToArgument(i)->getRank()==1 ) {
        args[i-argstart]=getPntrToArgument(i)->get(current);
      } else {
        args[i-argstart] = getPntrToArgument(i)->get();
      }
    }
  }
  // Calculate the function and its derivatives
  std::vector<double> vals( getNumberOfComponents() );
  Matrix<double> derivatives( getNumberOfComponents(), args.size() );
  myfunc.calc( this, args, vals, derivatives );
  // And set the values
  for(unsigned i=0; i<vals.size(); ++i) {
    myvals.addValue( getConstPntrToComponent(i)->getPositionInStream(), vals[i] );
  }
  // Return if we are not computing derivatives
  if( doNotCalculateDerivatives() ) {
    return;
  }
  // And now compute the derivatives
  // Second condition here is only not true if actionInChain()==True if
  // input arguments the only non-constant objects in input are scalars.
  // In that case we can use the non chain version to calculate the derivatives
  // with respect to the scalar.
  if( actionInChain() ) {
    for(unsigned j=0; j<args.size(); ++j) {
      unsigned istrn = getPntrToArgument(argstart+j)->getPositionInStream();
      if( stored_arguments[argstart+j] ) {
        unsigned task_index = myvals.getTaskIndex();
        if( getPntrToArgument(argstart+j)->getRank()==0 ) {
          task_index=0;
        }
        myvals.addDerivative( istrn, task_index, 1.0 );
        myvals.updateIndex( istrn, task_index );
      }
      unsigned arg_deriv_s = arg_deriv_starts[argstart+j];
      for(unsigned k=0; k<myvals.getNumberActive(istrn); ++k) {
        unsigned kind=myvals.getActiveIndex(istrn,k);
        for(int i=0; i<getNumberOfComponents(); ++i) {
          unsigned ostrn=getConstPntrToComponent(i)->getPositionInStream();
          myvals.addDerivative( ostrn, arg_deriv_s + kind, derivatives(i,j)*myvals.getDerivative( istrn, kind ) );
        }
      }
      // Ensure we only store one lot of derivative indices
      bool found=false;
      ActionWithValue* aav=getPntrToArgument(argstart+j)->getPntrToAction();
      for(unsigned k=0; k<j; ++k) {
        if( arg_deriv_starts[argstart+k]==arg_deriv_s ) {
          if( getPntrToArgument(argstart+k)->getPntrToAction()!=aav ) {
            ActionWithVector* av = dynamic_cast<ActionWithVector*>( getPntrToArgument(argstart+j)->getPntrToAction() );
            if( av ) {
              for(int i=0; i<getNumberOfComponents(); ++i) {
                av->updateAdditionalIndices( getConstPntrToComponent(i)->getPositionInStream(), myvals );
              }
            }
          }
          found=true;
          break;
        }
      }
      if( found ) {
        continue;
      }
      for(unsigned k=0; k<myvals.getNumberActive(istrn); ++k) {
        unsigned kind=myvals.getActiveIndex(istrn,k);
        for(int i=0; i<getNumberOfComponents(); ++i) {
          unsigned ostrn=getConstPntrToComponent(i)->getPositionInStream();
          myvals.updateIndex( ostrn, arg_deriv_s + kind );
        }
      }
    }
  } else {
    unsigned base=0;
    for(unsigned j=0; j<args.size(); ++j) {
      if( getPntrToArgument(argstart+j)->getRank()==1 ) {
        for(int i=0; i<getNumberOfComponents(); ++i) {
          unsigned ostrn=getConstPntrToComponent(i)->getPositionInStream();
          myvals.addDerivative( ostrn, base+current, derivatives(i,j) );
          myvals.updateIndex( ostrn, base+current );
        }
      } else {
        for(int i=0; i<getNumberOfComponents(); ++i) {
          unsigned ostrn=getConstPntrToComponent(i)->getPositionInStream();
          myvals.addDerivative( ostrn, base, derivatives(i,j) );
          myvals.updateIndex( ostrn, base );
        }
      }
      base += getPntrToArgument(argstart+j)->getNumberOfValues();
    }
  }
}

template <class T>
unsigned FunctionOfVector<T>::getNumberOfFinalTasks() {
  unsigned nelements=0, argstart=myfunc.getArgStart();
  for(unsigned i=argstart; i<getNumberOfArguments(); ++i) {
    plumed_assert( getPntrToArgument(i)->getRank()<2 );
    if( getPntrToArgument(i)->getRank()==1 ) {
      if( nelements>0 ) {
        // if( getPntrToArgument(i)->isTimeSeries() && getPntrToArgument(i)->getShape()[0]<nelements ) nelements=getPntrToArgument(i)->isTimeSeries();
        // else
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
void FunctionOfVector<T>::areAllTasksRequired( std::vector<ActionWithVector*>& task_reducing_actions ) {
  if( task_reducing_actions.size()==0 ) {
    return;
  }
  if( !myfunc.allComponentsRequired( getArguments(), task_reducing_actions ) ) {
    task_reducing_actions.push_back(this);
  }
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
  // Everything is done elsewhere
  if( actionInChain() ) {
    return;
  }
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
