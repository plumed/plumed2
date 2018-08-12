/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2018 The plumed team
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
#include "ParallelPlumedActions.h"
#include "ActionAtomistic.h"
#include "ActionWithArguments.h"
#include "ActionRegister.h"
#include "PlumedMain.h"
#include "ActionSet.h"

namespace PLMD {

PLUMED_REGISTER_ACTION(ParallelPlumedActions,"PLUMED_VECTOR")

void ParallelPlumedActions::registerKeywords( Keywords& keys ) {
  Action::registerKeywords(keys); ActionWithValue::registerKeywords(keys);
  keys.add("numbered","INPUT","the name of a PLUMED input file or alternatively the input for a single PLUMED action");
  keys.reset_style("INPUT","compulsory");
}

ParallelPlumedActions::ParallelPlumedActions(const ActionOptions&ao):
  Action(ao),
  ActionWithValue(ao),
  nderivatives(0)
{
  // Read the plumed input
  for(int i=1;; i++) {
      std::string input;
      if( !parseNumbered("INPUT",i,input) ) break;
      // Get the current size of the action set
      unsigned a_start = plumed.getActionSet().size();
      // Read the input to this command
      std::vector<std::string> input_words=Tools::getWords(input);
      if( input_words.size()==1 ) plumed.readInputFile(input);
      else {
        std::string remainder = input;
        while( remainder.find(";")!=std::string::npos ) {
            std::size_t semi = remainder.find_first_of(';');
            plumed.readInputLine( remainder.substr(0,semi) + " CALLER=" + getLabel() ); 
            remainder = remainder.substr(semi+1); 
        }
        plumed.readInputLine( remainder );
      } 
      // Now copy the actions that were just created to local storage
      action_lists.push_back( std::pair<unsigned,unsigned>( a_start, plumed.getActionSet().size() ) );
      // Make sure these actions are not run by the PlumedMain object
      for(unsigned j=a_start;j<plumed.getActionSet().size();++j) {
          // Add in any dependencies on other actions
          for(const auto & p : plumed.getActionSet()[j].get()->getDependencies() ) {
              if( p->getCaller()=="plumedmain" ) addDependency(p); 
          } 
      }
      // Check that final command in input is just a scalar
      const ActionWithValue* av=dynamic_cast<const ActionWithValue*>( plumed.getActionSet()[action_lists[i-1].second-1].get() );
      if( !av ) error("final action in each created set of action should have a value");
      if( av->getNumberOfComponents()!=1 ) error("final action in each created set of actions should calculate only one scalar");
      if( av->copyOutput(0)->isPeriodic() ) error("final action should not have periodic domain");
      if( av->copyOutput(0)->getRank()!=0 ) error("final action should calculate a scalar");
      // Now store the value that we need
      valuesToGet.push_back( av->copyOutput(0) );
  } 
  // Create a value to hold the output from the plumed inputs
  std::vector<unsigned> shape(1); shape[0]=action_lists.size();
  addValue( shape ); setNotPeriodic(); forcesToApply.resize( action_lists.size() );
  // Now create some tasks
  for(unsigned i=0;i<action_lists.size();++i) addTaskToList(i);
  // Now work out the number of derivatives we have in this action
  der_starts.resize( action_lists.size() );
  for(unsigned i=0;i<action_lists.size();++i) {
      ActionWithValue* av=dynamic_cast<ActionWithValue*>( plumed.getActionSet()[action_lists[i].second-1].get() );
      der_starts[i]=nderivatives; nderivatives += av->getNumberOfDerivatives();
  }
}

void ParallelPlumedActions::turnOnDerivatives() {
  ActionWithValue::turnOnDerivatives();
  // Check all action lists have the same size (ultimately what to check for single chains in each action)
  for(unsigned i=0;i<action_lists.size();++i) {
      // Turn on derivatives for all associated actions
      for(unsigned j=action_lists[i].first;j<action_lists[i].second;++j) {
          ActionWithValue* av=dynamic_cast<ActionWithValue*>( plumed.getActionSet()[j].get() );
          if(av) {
             av->turnOnDerivatives();
             if( j>action_lists[i].first && !av->actionInChain() ) {
                 error("cannot use derivatives with multiple chains in input");
             }
          } 
      }
  }
}

unsigned ParallelPlumedActions::getNumberOfDerivatives() const {
  return nderivatives;
}

void ParallelPlumedActions::clearDerivatives( const bool& force ) {
  ActionWithValue::clearDerivatives( force );
  for(unsigned i=0;i<action_lists.size();++i) {
      for(unsigned j=action_lists[i].first;j<action_lists[i].second;++j) {
          Action*p=plumed.getActionSet()[j].get();
          ActionWithValue*av=dynamic_cast<ActionWithValue*>(p);
          ActionAtomistic*aa=dynamic_cast<ActionAtomistic*>(p);
          {
            if(av) av->clearInputForces();
            if(av) av->clearDerivatives();
          }
          {
            if(aa) aa->clearOutputForces();
            if(aa) aa->retrieveAtoms();
          }
      }
  }
}

void ParallelPlumedActions::calculate() {
  plumed_assert( !actionInChain() );
  runAllTasks();
}

void ParallelPlumedActions::performTask( const unsigned& task_index, MultiValue& myvals ) const {
  // Run calculate for these actions
  for(unsigned i=action_lists[task_index].first;i<action_lists[task_index].second;++i) plumed.getActionSet()[i].get()->activate();
  for(unsigned i=action_lists[task_index].first;i<action_lists[task_index].second;++i) {
      Action* aa = plumed.getActionSet()[i].get(); aa->calculate(); aa->deactivate();
      
  }
  // Retrieve the value
  myvals.setValue( getPntrToOutput(0)->getPositionInStream(), valuesToGet[task_index]->get() );
  // Set the derivatives
  if( !doNotCalculateDerivatives() ) {
      unsigned nstart=der_starts[task_index]; unsigned jval = getPntrToOutput(0)->getPositionInStream();
      for(unsigned i=0;i<valuesToGet[task_index]->getNumberOfDerivatives();++i){
         myvals.addDerivative( jval, nstart + i, valuesToGet[task_index]->getDerivative(i) ); myvals.updateIndex( jval, nstart + i );
      }
  }
}

void ParallelPlumedActions::setForcesOnPlumedActions( const std::vector<double>& forces, unsigned& start ) {
  // Pass these forces to the underlying PLUMED objects
  for(unsigned i=0;i<action_lists.size();++i) {
      ActionAtomistic* aa = dynamic_cast<ActionAtomistic*>( plumed.getActionSet()[action_lists[i].first].get() );
      if( aa ) aa->setForcesOnAtoms( forces, start );
      ActionWithArguments* av = dynamic_cast<ActionWithArguments*>( plumed.getActionSet()[action_lists[i].first].get() );
      if( av ) av->setForcesOnArguments( forces, start );   
  }
}

void ParallelPlumedActions::apply() {
  if( doNotCalculateDerivatives() ) return;
  std::fill(forcesToApply.begin(),forcesToApply.end(),0);
  if( getForcesFromValues( forcesToApply ) ) { 
      error("This is untested and I am pretty sure it doesn't work in parallel - try at your own peril");
      // Add forces to final values
      for(unsigned i=0;i<action_lists.size();++i) {
          const ActionWithValue* av=dynamic_cast<const ActionWithValue*>( plumed.getActionSet()[action_lists[i].second-1].get() );
          (av->copyOutput(0))->addForce( 0, forcesToApply[i] );
      }
      for(unsigned i=0;i<action_lists.size();++i) {
          for(unsigned j=action_lists[i].second-1;j<=action_lists[i].first;--j) plumed.getActionSet()[j].get()->apply();
      }
  }
  // Now do apply for all actions
  for(unsigned i=0;i<action_lists.size();++i) {
      for(unsigned j=action_lists[i].first;j<action_lists[i].second;++j) {
          ActionAtomistic* a=dynamic_cast<ActionAtomistic*>( plumed.getActionSet()[j].get() );
          if( a ) a->applyForces();
      }
  }
}


}
