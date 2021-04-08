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
  forcesWereSet(false),
  nderivatives(0)
{
  // Read the plumed input
  unsigned ncols=0; bool periodic=false; std::string pmin, pmax;
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
        plumed.readInputLine( remainder + " CALLER=" + getLabel() );
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
      if( av->getNumberOfComponents()!=1 && av->getName()!="RMSD" ) error("final action in each created set of actions should calculate only one scalar");
      if( i==1 && av->copyOutput(0)->isPeriodic() ) { periodic=true; av->copyOutput(0)->getDomain( pmin, pmax ); }
      else if( !periodic && av->copyOutput(0)->isPeriodic() ) error("mismatched periodicities in inputs");
      else if( periodic && !av->copyOutput(0)->isPeriodic() ) error("mismatched periodicities in inputs");
      else if( periodic ) { 
         std::string cmin, cmax; av->copyOutput(0)->getDomain( cmin, cmax ); 
         if( cmin!=pmin || cmax!=pmax ) error("mismatched domains in inputs");
      }
      if( av->copyOutput(0)->getRank()!=0 ) {
          av->copyOutput(0)->buildDataStore( getLabel() );
          if( i==1 ) ncols = av->copyOutput(0)->getNumberOfValues( getLabel() );
          else if( ncols!=av->copyOutput(0)->getNumberOfValues( getLabel() ) ) {
             error("mismatched sizes of values");
          }
      }
      // Now store the value that we need
      valuesToGet.push_back( av->copyOutput(0) );
  } 
  // Create a value to hold the output from the plumed inputs
  std::vector<unsigned> shape; 
  if( ncols==0 ){ shape.resize(1); shape[0]=action_lists.size(); }
  else { shape.resize(2); shape[0]=action_lists.size(); shape[1]=ncols; }
  addValue( shape ); if( periodic ) setPeriodic( pmin, pmax ); else setNotPeriodic(); 
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

unsigned ParallelPlumedActions::getNumberOfColumns() const {
  plumed_assert( getPntrToOutput(0)->getRank()==2 );
  return getPntrToOutput(0)->getShape()[1];
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
  forcesWereSet=false;
}

void ParallelPlumedActions::getInfoForGridHeader( std::string& gtype, std::vector<std::string>& argn, std::vector<std::string>& min,
                                                  std::vector<std::string>& max, std::vector<unsigned>& nbin,
                                                  std::vector<double>& spacing, std::vector<bool>& pbc, const bool& dumpcube ) const {
  gtype="flat"; nbin[0] = getPntrToOutput(0)->getNumberOfValues( getLabel() ); spacing[0] = 1;
}

void ParallelPlumedActions::getGridPointIndicesAndCoordinates( const unsigned& ind, std::vector<unsigned>& indices, std::vector<double>& coords ) const {
  coords[0] = ind;
}

void ParallelPlumedActions::activate() {
  // Activate all actions here so that atoms are collected for all actions
  for(unsigned i=0;i<action_lists.size();++i) {
      for(unsigned j=action_lists[i].first;j<action_lists[i].second;++j) plumed.getActionSet()[j].get()->activate();
  }
  Action::activate();
}

void ParallelPlumedActions::calculate() {
  plumed_assert( !actionInChain() );
  runAllTasks();
}

void ParallelPlumedActions::prepareForTasks( const unsigned& nactive, const std::vector<unsigned>& pTaskList ) {
  // Deactivate all actions that are called by this thing -- currently all are active on all nodes (was done in activate)
  for(unsigned i=0;i<action_lists.size();++i) {
      for(unsigned j=action_lists[i].first;j<action_lists[i].second;++j) plumed.getActionSet()[j].get()->deactivate();
  }
  // Now activate tasks that are actually on this node
  unsigned stride=comm.Get_size(); unsigned rank=comm.Get_rank();
  if( runInSerial() ) { stride=1; rank=0; }
  for(unsigned i=rank;i<nactive;i+=stride) {
      unsigned itask = pTaskList[i];
      for(unsigned j=action_lists[itask].first;j<action_lists[itask].second;++j) plumed.getActionSet()[j].get()->activate();
  }
}

void ParallelPlumedActions::performTask( const unsigned& task_index, MultiValue& myvals ) const {
  // Run calculate for these actions
  for(unsigned i=action_lists[task_index].first;i<action_lists[task_index].second;++i) {
      Action* aa = plumed.getActionSet()[i].get(); aa->calculate(); aa->deactivate();
  }
  // Retrieve the value
  if( getPntrToOutput(0)->getRank()==1 ) { 
      myvals.setValue( getPntrToOutput(0)->getPositionInStream(), valuesToGet[task_index]->get() );
      // Set the derivatives
      if( !doNotCalculateDerivatives() ) {
          unsigned nstart=der_starts[task_index]; unsigned jval = getPntrToOutput(0)->getPositionInStream();
          for(unsigned i=0;i<valuesToGet[task_index]->getNumberOfDerivatives();++i){
             myvals.addDerivative( jval, nstart + i, valuesToGet[task_index]->getDerivative(i) ); myvals.updateIndex( jval, nstart + i );
          }
      }
  } else {
      unsigned nvals = getPntrToOutput(0)->getShape()[1], matind = getPntrToOutput(0)->getPositionInMatrixStash();
      for(unsigned i=0;i<nvals;++i) { myvals.stashMatrixElement( matind, i, valuesToGet[task_index]->get(i) ); }
  }
}

void ParallelPlumedActions::setForcesOnPlumedActions( const std::vector<double>& forces, unsigned& start ) {
  // Pass these forces to the underlying PLUMED objects
  for(unsigned i=0;i<action_lists.size();++i) {
      ActionAtomistic* aa = dynamic_cast<ActionAtomistic*>( plumed.getActionSet()[action_lists[i].first].get() );
      if( aa ) aa->setForcesOnAtoms( forces, start );
      ActionWithArguments* av = dynamic_cast<ActionWithArguments*>( plumed.getActionSet()[action_lists[i].first].get() );
      if( av ) av->setForcesOnArguments( 0, forces, start );   
  }
  forcesWereSet=true;
}

void ParallelPlumedActions::apply() {
  if( doNotCalculateDerivatives() ) return;

  // Check what tasks need to be redone on all nodes
  if( !forcesWereSet ) {
      unsigned ncols = 1; if( getPntrToOutput(0)->getRank()==2 ) ncols = getPntrToOutput(0)->getShape()[1];
      for(unsigned i=0;i<action_lists.size();++i) {
          bool hasforce=false; 
          for(unsigned j=0;j<ncols;++j){ if( fabs( getPntrToOutput(0)->getForce(i*ncols+j) )>epsilon) hasforce=true; }
          if( hasforce ) {
              for(unsigned j=action_lists[i].first;j<action_lists[i].second;++j) {
                  Action* p=plumed.getActionSet()[j].get(); p->activate(); 
                  ActionWithValue* av=dynamic_cast<ActionWithValue*>(p);
                  if(av) av->clearDerivatives(); 
                  p->calculate();
              }
          }
      }
      // Run tasks and apply
      for(unsigned i=0;i<action_lists.size();++i) {
          if( !plumed.getActionSet()[action_lists[i].first].get()->isActive() ) continue;
          // Add forces on this value to underlying values
          for(unsigned j=0;j<ncols;++j) valuesToGet[i]->addForce( j, getPntrToOutput(0)->getForce(i*ncols+j) );
          // Do apply
          for(unsigned j=action_lists[i].second-1;j>=action_lists[i].first;--j) plumed.getActionSet()[j].get()->apply();
          // Deactivate everything we have activated
          for(unsigned j=action_lists[i].first;j<action_lists[i].second;++j) plumed.getActionSet()[j].get()->deactivate(); 
      }
  }
  // Now do apply for all actions
  for(unsigned i=0;i<action_lists.size();++i) {
      for(unsigned j=action_lists[i].first;j<action_lists[i].second;++j) {
          // Apply forces on atoms
          ActionAtomistic* a=dynamic_cast<ActionAtomistic*>( plumed.getActionSet()[j].get() );
          if( a ) a->applyForces();
      }
  }
}


}
