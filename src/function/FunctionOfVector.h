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

#include "core/ActionWithValue.h"
#include "core/ActionWithArguments.h"
#include "core/CollectFrames.h"
#include "core/ActionSetup.h"
#include "tools/Matrix.h"

namespace PLMD {
namespace function {

template <class T>
class FunctionOfVector : 
public ActionWithValue,
public ActionWithArguments
{
private:
/// The forces that we get from the values
  std::vector<double> forcesToApply;
/// The function that is being computed
  T myfunc;
/// The number of derivatives for this action
  unsigned nderivatives;
public:
  static void registerKeywords(Keywords&);
/// This method is used to run the calculation with functions such as highest/lowest and sort.  
/// It is static so we can reuse the functionality in FunctionOfMatrix
  static void runSingleTaskCalculation( const Value* arg, ActionWithValue* action, T& f );
  explicit FunctionOfVector(const ActionOptions&);
/// Get the size of the task list at the end of the run
  unsigned getNumberOfFinalTasks() override ;
/// Check if derivatives are available
  void turnOnDerivatives() override;
/// Get the number of derivatives for this action
  unsigned getNumberOfDerivatives() const override ;
/// Get the label to write in the graph
  std::string writeInGraph() const override { return myfunc.getGraphInfo( getName() ); }
/// This updates the number of tasks we need to do if there is a time seeries
  void buildCurrentTaskList( bool& forceAllTasks, std::vector<std::string>& actionsThatSelectTasks, std::vector<unsigned>& tflags ) override;
  void calculate() override;
  void update() override;
  void runFinalJobs() override;
/// Calculate the function
  void performTask( const unsigned& current, MultiValue& myvals ) const override ;
/// Add the forces 
  void apply() override;
};

template <class T>
void FunctionOfVector<T>::registerKeywords(Keywords& keys ) {
  Action::registerKeywords(keys); ActionWithValue::registerKeywords(keys); ActionWithArguments::registerKeywords(keys); keys.use("ARG");
  keys.reserve("compulsory","PERIODIC","if the output of your function is periodic then you should specify the periodicity of the function.  If the output is not periodic you must state this using PERIODIC=NO");
  keys.add("hidden","NO_ACTION_LOG","suppresses printing from action on the log");
  T tfunc; tfunc.registerKeywords( keys );
}

template <class T>
FunctionOfVector<T>::FunctionOfVector(const ActionOptions&ao):
Action(ao),
ActionWithValue(ao),
ActionWithArguments(ao),
nderivatives(getNumberOfScalarArguments())
{
  // Get the shape of the output
  std::vector<unsigned> shape(1); shape[0]=getNumberOfFinalTasks();
  // Read the input and do some checks
  myfunc.read( this );
  // Create the task list
  if( myfunc.doWithTasks() ) {
      for(unsigned i=0;i<shape[0];++i) addTaskToList(i);
  } else { plumed_assert( getNumberOfArguments()==1 ); done_over_stream=false; getPntrToArgument(0)->buildDataStore( getLabel() ); }
  // Check we are not calculating a sum
  if( myfunc.zeroRank() ) shape.resize(0);  
  // Get the names of the components
  std::vector<std::string> components( keywords.getAllOutputComponents() );
  // Create the values to hold the output
  if( components.size()==0 && myfunc.zeroRank() ) addValueWithDerivatives( shape );
  else if( components.size()==0 ) addValue( shape );
  else { 
    std::vector<std::string> str_ind( myfunc.getComponentsPerLabel() );
    for(unsigned i=0;i<components.size();++i) {
        if( str_ind.size()>0 ) {
            for(unsigned j=0;j<str_ind.size();++j) addComponent( components[i] + str_ind[j], shape );
        } else if( components[i].find_first_of("_")!=std::string::npos ) {
            if( getNumberOfArguments()==1 ) addValue( shape ); 
            else { for(unsigned i=0; i<getNumberOfArguments(); ++i) addComponent( getPntrToArgument(i)->getName() + components[i], shape ); }
        } else addComponent( components[i], shape );
    } 
  }
  // Set the periodicities of the output components
  myfunc.setPeriodicityForOutputs( this );
  // Check if we can put the function in a chain
  bool hasscalar=false;
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
      CollectFrames* ab=dynamic_cast<CollectFrames*>( getPntrToArgument(i)->getPntrToAction() );
      // No chains if we have a scalar input that is not calculated in setup 
      if( getPntrToArgument(i)->getRank()==0 ) {
          ActionSetup* as=dynamic_cast<ActionSetup*>( getPntrToArgument(i)->getPntrToAction() );
          if(!as) hasscalar=true; else getPntrToArgument(i)->buildDataStore( getLabel() );
      // No chains if we are doing some analysis of a stored time series
      } else if(ab) {
          if( ab->hasClear() ) { hasscalar=true; getPntrToArgument(i)->buildDataStore( getLabel() ); }
      }
  }
  if( myfunc.doWithTasks() && !hasscalar && distinct_arguments.size()>0 ) nderivatives = setupActionInChain(0); 
}

template <class T>
void FunctionOfVector<T>::turnOnDerivatives() {
  if( !myfunc.derivativesImplemented() ) error("derivatives have not been implemended for " + getName() );
  ActionWithValue::turnOnDerivatives(); 
}

template <class T>
unsigned FunctionOfVector<T>::getNumberOfDerivatives() const {
  return nderivatives;
}

template <class T>
void FunctionOfVector<T>::performTask( const unsigned& current, MultiValue& myvals ) const {
  std::vector<double> args( getNumberOfArguments() );
  if( actionInChain() ) {
      for(unsigned i=0;i<getNumberOfArguments();++i) {
          if(  getPntrToArgument(i)->getRank()==0 ) args[i] = getPntrToArgument(i)->get();
          else if( !getPntrToArgument(i)->valueHasBeenSet() ) args[i] = myvals.get( getPntrToArgument(i)->getPositionInStream() );
          else args[i] = getPntrToArgument(i)->get( myvals.getTaskIndex() );
      }
  } else {
      for(unsigned i=0;i<getNumberOfArguments();++i) {
          if( getPntrToArgument(i)->getRank()==1 ) args[i]=getPntrToArgument(i)->get(current);
          else args[i] = getPntrToArgument(i)->get();
      }
  } 
  // Calculate the function and its derivatives
  std::vector<double> vals( getNumberOfComponents() ); Matrix<double> derivatives( getNumberOfComponents(), getNumberOfArguments() );
  myfunc.calc( this, args, vals, derivatives );
  // And set the values
  for(unsigned i=0;i<vals.size();++i) myvals.addValue( getPntrToOutput(i)->getPositionInStream(), vals[i] );
  // Return if we are not computing derivatives
  if( doNotCalculateDerivatives() ) return;

  // And now compute the derivatives
  if( actionInChain() ) {
      for(unsigned j=0;j<getNumberOfArguments();++j) {
          unsigned istrn = getArgumentPositionInStream( j, myvals );
          for(unsigned k=0; k<myvals.getNumberActive(istrn); ++k) {
              unsigned kind=myvals.getActiveIndex(istrn,k);
              for(unsigned i=0;i<getNumberOfComponents();++i) {
                  unsigned ostrn=getPntrToOutput(i)->getPositionInStream();
                  myvals.addDerivative( ostrn, arg_deriv_starts[j] + kind, derivatives(i,j)*myvals.getDerivative( istrn, kind ) );
              }
          }
          // Ensure we only store one lot of derivative indices
          bool found=false;
          for(unsigned k=0; k<j; ++k) {
              if( arg_deriv_starts[k]==arg_deriv_starts[j] ) { found=true; break; }
          }
          if( found ) continue;
          for(unsigned k=0; k<myvals.getNumberActive(istrn); ++k) {
              unsigned kind=myvals.getActiveIndex(istrn,k);
              for(unsigned i=0;i<getNumberOfComponents();++i) {
                  unsigned ostrn=getPntrToOutput(i)->getPositionInStream();
                  myvals.updateIndex( ostrn, arg_deriv_starts[j] + kind );
              }
          } 
      }
  } else {
      unsigned base=0;
      for(unsigned j=0;j<getNumberOfArguments();++j) { 
          if( getPntrToArgument(j)->getRank()==1 ) {
              for(unsigned i=0;i<getNumberOfComponents();++i) {
                  unsigned ostrn=getPntrToOutput(i)->getPositionInStream(); 
                  myvals.addDerivative( ostrn, base+current, derivatives(i,j) ); 
                  myvals.updateIndex( ostrn, base+current );
              }
          } else {
              for(unsigned i=0;i<getNumberOfComponents();++i) {
                  unsigned ostrn=getPntrToOutput(i)->getPositionInStream();
                  myvals.addDerivative( ostrn, base, derivatives(i,j) );
                  myvals.updateIndex( ostrn, base );
              }
          }
          base += getPntrToArgument(j)->getNumberOfValues();
      }
  }
}

template <class T>
unsigned FunctionOfVector<T>::getNumberOfFinalTasks() {
  unsigned nelements=0;
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
      plumed_assert( getPntrToArgument(i)->getRank()<2 );
      if( getPntrToArgument(i)->getRank()==1 ) {
          if( nelements>0 ) { 
              if( getPntrToArgument(i)->isTimeSeries() && getPntrToArgument(i)->getShape()[0]<nelements ) nelements=getPntrToArgument(i)->isTimeSeries();
              else if(getPntrToArgument(i)->getShape()[0]!=nelements ) error("all vectors input should have the same length");
          } else if( nelements==0 ) nelements=getPntrToArgument(i)->getShape()[0];
          plumed_assert( !getPntrToArgument(i)->hasDerivatives() );
      }
  }
  // The prefactor for average and sum is set here so the number of input scalars is guaranteed to be correct
  myfunc.setPrefactor( this, 1.0 );
  return nelements;
}

template <class T>
void FunctionOfVector<T>::buildCurrentTaskList( bool& forceAllTasks, std::vector<std::string>& actionsThatSelectTasks, std::vector<unsigned>& tflags ) {
  // If we have a time series then resize the task list so this action does all the required tasks
  plumed_assert( myfunc.doWithTasks() ); unsigned nstart = getFullNumberOfTasks(), ndata = 0;
  for(unsigned i=0;i<getNumberOfArguments();++i) {
      if( getPntrToArgument(i)->getRank()<=1 && getPntrToArgument(i)->isTimeSeries() ) ndata = getPntrToArgument(i)->getNumberOfValues();
  }
  if( nstart<ndata ) {
      for(unsigned i=nstart;i<ndata;++i) addTaskToList(i);
      std::vector<unsigned> shape(1); shape[0]=ndata;
      for(unsigned i=0;i<getNumberOfComponents();++i) {
          if( getPntrToOutput(i)->getRank()==1 ) getPntrToOutput(i)->setShape( shape );
      }
  }
  // Check if there is a sepcial method for setting up the task list
  if( !myfunc.defaultTaskListBuilder() ) {
      myfunc.buildTaskList( this, actionsThatSelectTasks );
  } else {
      bool safeToChain=true, atLeastOneRank=false;
      for(unsigned i=0;i<getNumberOfArguments();++i) {
          if( getPntrToArgument(i)->getRank()>0 ) atLeastOneRank=true;
          Action* myact = getPntrToArgument(i)->getPntrToAction();
          if( myact ) {
              std::string argact = myact->getLabel(); bool found=false;
              for(unsigned j=0;j<actionsThatSelectTasks.size();++j) {
                  if( argact==actionsThatSelectTasks[j] ){ found=true; break; }
              }
              if( !found ) safeToChain=false;
          } else safeToChain=false;
      }
      plumed_assert( atLeastOneRank );
      if( safeToChain ) actionsThatSelectTasks.push_back( getLabel() );
  }
}

template <class T>
void FunctionOfVector<T>::runSingleTaskCalculation( const Value* arg, ActionWithValue* action, T& f ) {
  // This is used if we are doing sorting actions on a single vector
  unsigned nv = arg->getNumberOfValues(); std::vector<double> args( nv );
  for(unsigned i=0;i<nv;++i) args[i] = arg->get(i);
  std::vector<double> vals( action->getNumberOfComponents() ); Matrix<double> derivatives( action->getNumberOfComponents(), nv );
  f.calc( NULL, args, vals, derivatives );
  for(unsigned i=0;i<vals.size();++i) action->copyOutput(i)->set( vals[i] );
  // Return if we are not computing derivatives
  if( action->doNotCalculateDerivatives() ) return;
  // Now set the derivatives
  for(unsigned j=0; j<nv; ++j) {
      for(unsigned i=0;i<vals.size();++i) action->copyOutput(i)->setDerivative( j, derivatives(i,j) );
  }
} 

template <class T>
void FunctionOfVector<T>::calculate() {
  // Everything is done elsewhere
  if( actionInChain() ) return;
  // This is done if we are calculating a function of multiple cvs
  if( getFullNumberOfTasks()>0 ) runAllTasks();
  // This is used if we are doing sorting actions on a single vector
  else runSingleTaskCalculation( getPntrToArgument(0), this, myfunc );
}

template <class T>
void FunctionOfVector<T>::update() {
  if( skipUpdate() || actionInChain() ) return;
  plumed_dbg_assert( !actionInChain() ); plumed_assert( myfunc.doWithTasks() );
  if( getFullNumberOfTasks()>0 ) runAllTasks();
}
  
template <class T>
void FunctionOfVector<T>::runFinalJobs() {
  if( skipUpdate() || actionInChain() ) return;
  plumed_assert( myfunc.doWithTasks() );
  resizeForFinalTasks(); runAllTasks();
}

template <class T>
void FunctionOfVector<T>::apply() {
  if( doNotCalculateDerivatives() ) return;
          
  if( forcesToApply.size()!=getNumberOfDerivatives() ) forcesToApply.resize( getNumberOfDerivatives() );
  
  std::fill(forcesToApply.begin(),forcesToApply.end(),0); unsigned ss=0;
  if( getForcesFromValues( forcesToApply ) ) setForcesOnArguments( 0, forcesToApply, ss );
}

}
}
#endif
