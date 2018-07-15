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
  ActionWithValue(ao)
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
      else plumed.readInputLine( input );
      // Now copy the actions that were just created to local storage
      action_lists.push_back( std::pair<unsigned,unsigned>( a_start, plumed.getActionSet().size() ) );
      // Make sure these actions are not run by the PlumedMain object
      for(unsigned j=a_start;j<plumed.getActionSet().size();++j) {
          plumed.getActionSet()[j].get()->setCallingAction( getLabel() );
      }
      // Check that final command in input is just a scalar
      const ActionWithValue* av=dynamic_cast<const ActionWithValue*>( plumed.getActionSet()[action_lists[i-1].second-1].get() );
      if( !av ) error("final action in each created set of action should have a value");
      if( av->getNumberOfComponents()!=1 ) error("final action in each created set of actions should calculate only one scalar");
      if( av->copyOutput(0)->isPeriodic() ) error("final action should not have periodic domain");
      if( av->copyOutput(0)->getRank()!=0 ) error("final action should calculate a scalar");
  } 
  // Create a value to hold the output from the plumed inputs
  std::vector<unsigned> shape(1); shape[0]=action_lists.size();
  addValue( shape ); setNotPeriodic();
  // Now create some tasks
  for(unsigned i=0;i<action_lists.size();++i) addTaskToList(i);
}

unsigned ParallelPlumedActions::getNumberOfDerivatives() const {
  return 0;
}

void ParallelPlumedActions::calculate() {
  plumed_assert( !actionInChain() );
  runAllTasks();
}

void ParallelPlumedActions::performTask( const unsigned& task_index, MultiValue& myvals ) const {
  for(unsigned i=action_lists[task_index].first;i<action_lists[task_index].second;++i) {
      Action*p=plumed.getActionSet()[i].get();
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
      p->calculate();
  }
  const ActionWithValue* av=dynamic_cast<const ActionWithValue*>( plumed.getActionSet()[action_lists[task_index].second-1].get() );
  myvals.setValue( getPntrToOutput(0)->getPositionInStream(), av->copyOutput(0)->get() );
}

void ParallelPlumedActions::apply() {

}


}
