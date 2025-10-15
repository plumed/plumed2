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
#include "ActionWithVector.h"
#include "ActionWithMatrix.h"
#include "PlumedMain.h"
#include "ActionSet.h"
#include "tools/OpenMP.h"
#include "tools/Communicator.h"

namespace PLMD {

void ActionWithVector::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  keys.remove("NUMERICAL_DERIVATIVES");
  ActionWithArguments::registerKeywords( keys );
}

ActionWithVector::ActionWithVector(const ActionOptions&ao):
  Action(ao),
  ActionAtomistic(ao),
  ActionWithValue(ao),
  ActionWithArguments(ao) {
  if( !keywords.exists("MASKED_INPUT_ALLOWED") ) {
    for(unsigned i=0; i<getNumberOfArguments(); ++i) {
      ActionWithVector* av = dynamic_cast<ActionWithVector*>( getPntrToArgument(i)->getPntrToAction() );
      if( av && av->getNumberOfMasks()>=0 ) {
        nmask=0;
      }
    }
  }

  if( keywords.exists("MASK") ) {
    std::vector<Value*> mask;
    parseArgumentList("MASK",mask);
    if( mask.size()>0 ) {
      if( nmask>=0 && getNumberOfArguments()==1 ) {
        ActionWithVector* av=dynamic_cast<ActionWithVector*>( getPntrToArgument(0)->getPntrToAction() );
        plumed_massert( av, "input should be a vector from ActionWithVector" );
        unsigned j=0, nargs = av->getNumberOfArguments();
        for(unsigned i=nargs-av->nmask; i<nargs; ++i) {
          if( av->getPntrToArgument(i)!=mask[j] ) {
            error("the masks in subsequent actions do not match");
          }
          j++;
        }
      }
      if( getNumberOfArguments()>0 && getName().find("EVALUATE_FUNCTION_FROM_GRID")==std::string::npos && getPntrToArgument(0)->hasDerivatives() ) {
        error("input for mask should be vector or matrix");
      } else if( mask[0]->getRank()==2 ) {
        if( mask.size()>1 ) {
          error("MASK should only have one argument");
        }
        log.printf("  only computing elements of matrix that correspond to non-zero elements of matrix %s \n", mask[0]->getName().c_str() );
      } else if( mask[0]->getRank()==1 ) {
        log.printf("  only computing elements of vector that correspond to non-zero elements of vectors %s", mask[0]->getName().c_str() );
        for(unsigned i=1; i<mask.size(); ++i) {
          if( mask[i]->getRank()!=1 ) {
            log.printf("\n");
            error("input to mask should be vector");
          }
          log.printf(", %s", mask[i]->getName().c_str() );
        }
        log.printf("\n");
      }
      std::vector<Value*> allargs( getArguments() );
      nmask=mask.size();
      for(unsigned i=0; i<mask.size(); ++i) {
        allargs.push_back( mask[i] );
      }
      requestArguments( allargs );
    }
  }
}

void ActionWithVector::lockRequests() {
  ActionAtomistic::lockRequests();
  ActionWithArguments::lockRequests();
}

void ActionWithVector::unlockRequests() {
  ActionAtomistic::unlockRequests();
  ActionWithArguments::unlockRequests();
}

void ActionWithVector::calculateNumericalDerivatives(ActionWithValue* av) {
  plumed_merror("cannot calculate numerical derivative for action " + getName() + " with label " + getLabel() );
}

void ActionWithVector::prepare() {
  active_tasks.resize(0);
}

int ActionWithVector::checkTaskIsActive( const unsigned& itask ) const {
  unsigned nargs = getNumberOfArguments();
  if( nargs==0 ) {
    return 1;
  } else if( nmask>0 ) {
    for(unsigned j=nargs-nmask; j<nargs; ++j) {
      Value* myarg = getPntrToArgument(j);
      if( myarg->getRank()==1 && !myarg->hasDerivatives() ) {
        if( fabs(myarg->get(itask))>0.0 ) {
          return 1;
        }
      } else if( myarg->getRank()==2 && !myarg->hasDerivatives() ) {
        if( myarg->getRowLength(itask)>0 ) {
          return 1;
        }
      } else {
        plumed_merror("only matrices and vectors should be used as masks");
      }
    }
  } else {
    for(unsigned i=0; i<nargs; ++i) {
      Value* myarg = getPntrToArgument(i);
      if( !myarg->isDerivativeZeroWhenValueIsZero() ) {
        return 1;
      }

      if( myarg->getRank()==0 ) {
        return 1;
      } else if( myarg->getRank()==1 && !myarg->hasDerivatives() ) {
        if( fabs(myarg->get(itask))>0.0 ) {
          return 1;
        }
      } else if( myarg->getRank()==2 && !myarg->hasDerivatives() ) {
        const unsigned ncol = myarg->getRowLength(itask);
        const unsigned base = itask*myarg->getNumberOfColumns();
        for(unsigned k=0; k<ncol; ++k) {
          if( fabs(myarg->get(base+k,false))>0.0 ) {
            return 1;
          }
        }
      } else if( myarg->getRank()>0 ) {
        return 1;
      } else {
        plumed_merror("should not be in action " + getName() );
      }
    }
  }
  return -1;
}

std::vector<unsigned>& ActionWithVector::getListOfActiveTasks( ActionWithVector* action ) {
  if( active_tasks.size()>0 ) {
    return active_tasks;
  }
  unsigned ntasks=0;
  getNumberOfTasks( ntasks );

  active_tasks.resize(0);
  active_tasks.reserve(ntasks);
  for(unsigned i=0; i<ntasks; ++i) {
    if( checkTaskIsActive(i)>0 ) {
//no resize are triggered, since we have reserved the number of tasks
      active_tasks.push_back(i);
    }
  }
  return active_tasks;
}

void ActionWithVector::getInputData( std::vector<double>& inputdata ) const {
  plumed_dbg_assert( getNumberOfAtoms()==0 );
  unsigned nargs = getNumberOfArguments();
  unsigned nmasks=getNumberOfMasks();
  // getNumberOfMasks(); returns nmask, that it is an int
  // nmasks cant be <0 (it is unsigned), so I check nmask for that
  if( nargs>=nmasks && nmask>0 ) {
    nargs = nargs - nmasks;
  }

  std::size_t total_args = 0;
  for(unsigned i=0; i<nargs; ++i) {
    total_args += getPntrToArgument(i)->getNumberOfStoredValues();
  }

  if( inputdata.size()!=total_args ) {
    inputdata.resize( total_args );
  }

  total_args = 0;
  for(unsigned i=0; i<nargs; ++i) {
    Value* myarg = getPntrToArgument(i);
    total_args+= myarg->assignValues(View{&inputdata[total_args],inputdata.size()-total_args});
  }
}

void ActionWithVector::transferStashToValues( const std::vector<unsigned>& partialTaskList, const std::vector<double>& stash ) {
  unsigned ntask = partialTaskList.size();
  unsigned ncomponents = getNumberOfComponents();
  for(unsigned i=0; i<ncomponents; ++i) {
    Value* myval = copyOutput(i);
    for(unsigned j=0; j<ntask; ++j) {
      myval->set( partialTaskList[j], stash[partialTaskList[j]*ncomponents+i] );
    }
  }
}

void ActionWithVector::transferForcesToStash( const std::vector<unsigned>& partialTaskList, std::vector<double>& stash ) const {
  unsigned ntask = partialTaskList.size();
  unsigned ncomponents = getNumberOfComponents();
  for(unsigned i=0; i<ncomponents; ++i) {
    auto myval = getConstPntrToComponent(i);
    for(unsigned j=0; j<ntask; ++j) {
      stash[partialTaskList[j]*ncomponents+i] = myval->getForce( partialTaskList[j] );
    }
  }
}

void ActionWithVector::getNumberOfTasks( unsigned& ntasks ) {
  if( ntasks==0 ) {
    if( getNumberOfArguments()==1 && getNumberOfComponents()==1 && getPntrToComponent(0)->getRank()==0 ) {
      if( !getPntrToArgument(0)->hasDerivatives() && getPntrToArgument(0)->getRank()==2 ) {
        ntasks = getPntrToArgument(0)->getShape()[0];
      } else {
        ntasks = getPntrToArgument(0)->getNumberOfValues();
      }
    } else {
      plumed_assert( getNumberOfComponents()>0 && getPntrToComponent(0)->getRank()>0 );
      if( getPntrToComponent(0)->hasDerivatives() ) {
        ntasks = getPntrToComponent(0)->getNumberOfValues();
      } else {
        ntasks = getPntrToComponent(0)->getShape()[0];
      }
    }
  }
  for(unsigned i=0; i<getNumberOfComponents(); ++i) {
    if( getPntrToComponent(i)->getRank()==0 ) {
      if( getNumberOfArguments()!=1 ) {
        error("mismatched numbers of tasks in streamed quantities");
      }
      if( getPntrToArgument(0)->hasDerivatives() && ntasks!=getPntrToArgument(0)->getNumberOfValues() ) {
        error("mismatched numbers of tasks in streamed quantities");
      } else if ( !getPntrToArgument(0)->hasDerivatives() && ntasks!=getPntrToArgument(0)->getShape()[0] ) {
        error("mismatched numbers of tasks in streamed quantities");
      }
    } else if( getPntrToComponent(i)->hasDerivatives() && ntasks!=getPntrToComponent(i)->getNumberOfValues() ) {
      error("mismatched numbers of tasks in streamed quantities");
    } else if( !getPntrToComponent(i)->hasDerivatives() && ntasks!=getPntrToComponent(i)->getShape()[0] ) {
      error("mismatched numbers of tasks in streamed quantities");
    }
  }
}

unsigned ActionWithVector::getNumberOfForceDerivatives() const {
  unsigned nforces=0;
  unsigned nargs = getNumberOfArguments();
  unsigned  nmasks = getNumberOfMasks();
  if( nargs>=nmasks && nmasks>0 ) {
    nargs = nargs - nmasks;
  }
  if( getNumberOfAtoms()>0 ) {
    nforces += 3*getNumberOfAtoms() + 9;
  }
  for(unsigned i=0; i<nargs; ++i) {
    nforces += getPntrToArgument(i)->getNumberOfStoredValues();
  }
  return nforces;
}

bool ActionWithVector::checkForForces() {
  if( getPntrToComponent(0)->getRank()==0 ) {
    return ActionWithValue::checkForForces();
  }

  // Check if there are any forces
  bool hasforce=false;
  for(unsigned i=0; i<getNumberOfComponents(); ++i) {
    if( getConstPntrToComponent(i)->getRank()>0 && getConstPntrToComponent(i)->forcesWereAdded() ) {
      hasforce=true;
      break;
    }
  }
  if( !hasforce ) {
    return false;
  }
  applyNonZeroRankForces( forcesForApply );
  return true;
}

void ActionWithVector::apply() {
  unsigned nf =  getNumberOfForceDerivatives();
  if( forcesForApply.size()!=nf ) {
    forcesForApply.resize( nf, 0 );
  }

  if( !checkForForces() ) {
    return;
  }
  // Find the top of the chain and add forces
  unsigned ind=0;
  addForcesOnArguments( 0, forcesForApply, ind );
  setForcesOnAtoms( forcesForApply, ind );
}

}
