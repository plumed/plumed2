/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2020 The plumed team
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
#include "core/ActionWithArguments.h"
#include "core/ActionWithValue.h"
#include "core/ActionRegister.h"

namespace PLMD {
namespace matrixtools {

class Stack :
  public ActionWithArguments,
  public ActionWithValue
{
private: 
  std::vector<double> forcesToApply;
  std::vector<std::string> actionsLabelsInChain;
public:
  static void registerKeywords( Keywords& keys );
  explicit Stack(const ActionOptions&);
  unsigned getNumberOfDerivatives() const override;
  unsigned getNumberOfColumns() const override;
  void calculate() override; 
  void getTasksForParent( const std::string& parent, std::vector<std::string>& actionsThatSelectTasks, std::vector<unsigned>& tflags );
  void performTask( const unsigned& current, MultiValue& myvals ) const override {}
  void gatherStoredValue( const unsigned& valindex, const unsigned& code, const MultiValue& myvals,
                          const unsigned& bufstart, std::vector<double>& buffer ) const override;
  void apply() override;
};

PLUMED_REGISTER_ACTION(Stack,"HSTACK")
PLUMED_REGISTER_ACTION(Stack,"VSTACK")

void Stack::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys ); 
  ActionWithValue::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys ); keys.use("ARG"); 
}

Stack::Stack(const ActionOptions&ao):
Action(ao),
ActionWithArguments(ao),
ActionWithValue(ao)
{
   unsigned nvals = 0;
   for(unsigned i=arg_ends[0];i<arg_ends[1];++i) nvals += getPntrToArgument(i)->getNumberOfValues( getLabel() );
   // Check for consistent numbers of values in other actions
   for(unsigned j=0;j<arg_ends.size()-1;++j) {
       unsigned tvals = 0;
       for(unsigned i=arg_ends[j];i<arg_ends[j+1];++i) tvals += getPntrToArgument(i)->getNumberOfValues( getLabel() );
       if( tvals!=nvals ) error("mismatch between number of values in each vector that is to be combined");
   }
   std::vector<unsigned> shape(2); 
   if( getName()=="HSTACK" ){ shape[0]=arg_ends.size()-1; shape[1]=nvals; }
   else if( getName()=="VSTACK" ) { shape[0]=nvals; shape[1]=arg_ends.size()-1; }
   else error("unknown type of stack object");

   std::vector<Value*> args( getArguments() ); requestArguments( args, false );
   addValue( shape ); setNotPeriodic(); getPntrToOutput(0)->alwaysStoreValues(); bool istimeseries=false;
   for(unsigned i=0; i<getNumberOfArguments();++i) {
      if( getPntrToArgument(0)->isTimeSeries() ) { istimeseries=true; break; }
   }
   if( istimeseries  ) {
      for(unsigned i=0; i<getNumberOfArguments();++i) {
          if( !getPntrToArgument(0)->isTimeSeries() ) error( "on argument is time series but " + getPntrToArgument(i)->getName() + " is not a time series");
      }
      getPntrToOutput(0)->makeTimeSeries();
  }
  forcesToApply.resize( shape[0]*shape[1] );
  for(unsigned i=0;i<nvals;++i) addTaskToList(i);
}

unsigned Stack::getNumberOfDerivatives() const {
  return 0;
}

unsigned Stack::getNumberOfColumns() const {
  return getPntrToOutput(0)->getShape()[1];
}

void Stack::getTasksForParent( const std::string& parent, std::vector<std::string>& actionsThatSelectTasks, std::vector<unsigned>& tflags ) {
  if( actionsLabelsInChain.size()==0 ) getAllActionLabelsInChain( actionsLabelsInChain );
  bool ignore = checkUsedOutsideOfChain( actionsLabelsInChain, parent, actionsThatSelectTasks, tflags );
}

void Stack::calculate() {
  // Everything is done elsewhere
  if( actionInChain() ) return;
  // Run all the tasks
  runAllTasks();
}

void Stack::gatherStoredValue( const unsigned& valindex, const unsigned& code, const MultiValue& myvals,
                               const unsigned& bufstart, std::vector<double>& buffer ) const {
  std::vector<double> argsh( arg_ends.size()-1 ); retrieveArguments( myvals, argsh, 0 );
  if( getName()=="HSTACK" ) {
      unsigned sss=getPntrToOutput(0)->getShape()[1];
      for(unsigned i=0;i<argsh.size();++i) buffer[ bufstart + code + i*sss ] += argsh[i];  
  } else if( getName()=="VSTACK" ) {
      for(unsigned i=0;i<argsh.size();++i) buffer[ bufstart + argsh.size()*code + i ] += argsh[i];
  } else plumed_error();
}

void Stack::apply() {
  if( !getPntrToOutput(0)->forcesWereAdded() ) return;

  Value* mat=getPntrToOutput(0); std::vector<unsigned> shape( mat->getShape() );
  if( getName()=="HSTACK" ) { 
      for(unsigned i=0; i<shape[0]; ++i) {
          for(unsigned j=0; j<shape[1]; ++j) forcesToApply[i*shape[1] + j] = mat->getForce( i*shape[1] + j );
      }   
  } else if(  getName()=="VSTACK" ) { 
      for(unsigned i=0; i<shape[1]; ++i) {
          for(unsigned j=0; j<shape[0]; ++j) forcesToApply[i*shape[0] + j] = mat->getForce( j*shape[1] + i );
      } 
  } else plumed_error(); 
  unsigned mm=0; setForcesOnArguments( 0, forcesToApply, mm );
}

}
}
