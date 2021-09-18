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

#include "FunctionBase.h"
#include "tools/Matrix.h"

namespace PLMD {
namespace function {

template <class T>
class FunctionOfVector : public FunctionBase {
private:
/// The function that is being computed
  T myfunc;
/// The number of derivatives for this action
  unsigned nderivatives;
public:
  static void registerKeywords(Keywords&);
  explicit FunctionOfVector(const ActionOptions&);
/// Get the size of the task list at the end of the run
  unsigned getNumberOfFinalTasks() override ;
/// Check if derivatives are available
  void turnOnDerivatives() override;
/// Get the number of derivatives for this action
  unsigned getNumberOfDerivatives() const override ;
/// This updates the number of tasks we need to do if there is a time seeries
  void resizeTimeSeriesTaskList() override ;
/// Calculate the function
  void performTask( const unsigned& current, MultiValue& myvals ) const override ;
/// Add the forces 
  void apply() override;
};

template <class T>
void FunctionOfVector<T>::registerKeywords(Keywords& keys ) {
  FunctionBase::registerKeywords( keys );
  keys.add("hidden","NO_ACTION_LOG","suppresses printing from action on the log");
  T tfunc; tfunc.registerKeywords( keys );
}

template <class T>
FunctionOfVector<T>::FunctionOfVector(const ActionOptions&ao):
Action(ao),
FunctionBase(ao),
nderivatives(getNumberOfScalarArguments())
{
  // Get the shape of the output
  std::vector<unsigned> shape(1); shape[0]=getNumberOfFinalTasks();
  // Create the task list
  for(unsigned i=0;i<shape[0];++i) addTaskToList(i);
  // Read the input and do some checks
  myfunc.read( this );
  // Check we are not calculating a sum
  if( myfunc.getRank()==0 ) shape.resize(0);  
  // Get the names of the components
  std::vector<std::string> components( keywords.getAllOutputComponents() );
  // Create the values to hold the output
  if( components.size()==0 && myfunc.getRank()==0 ) addValueWithDerivatives( shape );
  else if( components.size()==0 ) addValue( shape );
  else { for(unsigned i=0;i<components.size();++i) addComponent( components[i], shape ); }
  // Set the periodicities of the output components
  myfunc.setPeriodicityForOutputs( this );
  // Now setup the action in the chain if we can
  if( distinct_arguments.size()>0 ) nderivatives = setupActionInChain(0); 
  // Fix any functions that are time series
  bool timeseries=false;
  for(unsigned i=0;i<getNumberOfArguments();++i) {
      if( getPntrToArgument(i)->isTimeSeries() ) { timeseries=true; break; }
  }
  if( timeseries ) {
      for(unsigned i=0;i<getNumberOfComponents();++i) getPntrToOutput(i)->makeTimeSeries();
  }
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
          if( getPntrToArgument(i)->getRank()==1 ) args[i] = myvals.get( getPntrToArgument(i)->getPositionInStream() );
          else args[i] = getPntrToArgument(i)->get();
      }
  } else {
      for(unsigned i=0;i<getNumberOfArguments();++i) {
          if( getPntrToArgument(i)->getRank()==1 ) args[i]=getPntrToArgument(i)->get(current);
          else args[i] = getPntrToArgument(i)->get();
      }
  } 
  // Calculate the function and its derivatives
  std::vector<double> vals( getNumberOfComponents() ); Matrix<double> derivatives( getNumberOfComponents(), getNumberOfArguments() );
  myfunc.calc( args, vals, derivatives );
  // And set the values
  for(unsigned i=0;i<vals.size();++i) myvals.addValue( getPntrToOutput(i)->getPositionInStream(), vals[i] );
  // Return if we are not computing derivatives
  if( doNotCalculateDerivatives() ) return;

  // And now compute the derivatives
  if( actionInChain() ) {
      for(unsigned j=0;j<getNumberOfArguments();++j) {
          if( getPntrToArgument(j)->getRank()==1 ) {
              unsigned istrn = getPntrToArgument(j)->getPositionInStream();
              for(unsigned k=0; k<myvals.getNumberActive(istrn); ++k) {
                  unsigned kind=myvals.getActiveIndex(istrn,k);
                  for(unsigned i=0;i<getNumberOfComponents();++i) {
                      unsigned ostrn=getPntrToOutput(i)->getPositionInStream();
                      myvals.addDerivative( ostrn, arg_deriv_starts[j] + kind, derivatives(i,j)*myvals.getDerivative( istrn, kind ) );
                      myvals.updateIndex( ostrn, arg_deriv_starts[j] + kind );
                  }
              } 
          } else plumed_merror("not implemented this yet");
      }
  } else {
      for(unsigned j=0;j<getNumberOfArguments();++j) { 
          if( getPntrToArgument(j)->getRank()==1 ) {
              for(unsigned i=0;i<getNumberOfComponents();++i) {
                  unsigned ostrn=getPntrToOutput(i)->getPositionInStream(); 
                  myvals.addDerivative( ostrn, current, derivatives(i,j) ); 
                  myvals.updateIndex( ostrn, current );
              }
          } else plumed_error();
      }
  }
}

template <class T>
unsigned FunctionOfVector<T>::getNumberOfFinalTasks() {
  unsigned nelements=0;
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
      plumed_assert( getPntrToArgument(i)->getRank()<2 );
      if( getPntrToArgument(i)->getRank()==1 ) {
          if( nelements>0 && getPntrToArgument(i)->getShape()[0]!=nelements ) error("all vectors input should have the same length");
          else if( nelements==0 ) nelements=getPntrToArgument(i)->getShape()[0];
          plumed_assert( !getPntrToArgument(i)->hasDerivatives() );
      }
  }
  // The prefactor for average and sum is set here so the number of input scalars is guaranteed to be correct
  myfunc.setPrefactor( this );
  return nelements;
}

template <class T>
void FunctionOfVector<T>::resizeTimeSeriesTaskList() {
  unsigned nstart = getFullNumberOfTasks(), ndata = 0;
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
}

template <class T>
void FunctionOfVector<T>::apply() {
  applyNonGrid();
}

}
}
#endif
