/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2017 The plumed team
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
#include "MultiBias.h"


namespace PLMD {
namespace bias {

void MultiBias::registerKeywords( Keywords& keys ) {
  Action::registerKeywords(keys);
  ActionPilot::registerKeywords(keys);
  ActionWithValue::registerKeywords(keys);
  ActionWithArguments::registerKeywords(keys);
  keys.add("hidden","STRIDE","the frequency with which the forces due to the bias should be calculated.  This can be used to correctly set up multistep algorithms");
  componentsAreNotOptional(keys);
  keys.addOutputComponent("bias","default","the instantaneous value of the bias potential");
}

MultiBias::MultiBias(const ActionOptions&ao):
  Action(ao),
  ActionPilot(ao),
  ActionWithValue(ao),
  ActionWithArguments(ao),
  forcesToApply(getNumberOfScalarArguments(),0.0)
{
  createTasksFromArguments(); nderivatives = getNumberOfScalarArguments();
  if( distinct_arguments.size()>0 ) {
    // Create the chain of actions that will calculate the function
    nderivatives = setupActionInChain(0);
    // Set forces to apply to correct size
    forcesToApply.resize( nderivatives );
  }

  // Notice bias is always an object with zero rank
  addComponentWithDerivatives("bias"); componentIsNotPeriodic("bias");

  if(getStride()>1) {
    log<<"  multiple time step "<<getStride()<<" ";
    log<<cite("Ferrarotti, Bottaro, Perez-Villa, and Bussi, J. Chem. Theory Comput. 11, 139 (2015)")<<"\n";
  }
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    (getPntrToArgument(i)->getPntrToAction())->turnOnDerivatives();
  }

  turnOnDerivatives();
}

void MultiBias::calculate() {
  // Everything is done elsewhere
  if( actionInChain() ) return;
  // This is done if we are calculating a function of multiple cvs
  plumed_dbg_assert( getFullNumberOfTasks()>0 ); runAllTasks();
}

void MultiBias::performTask( const unsigned& current, MultiValue& myvals ) const {
  // Get the values of all the arguments
  std::vector<double> args( getNumberOfArgumentsPerTask() ); retrieveArguments( myvals, args, 0 );
  // Calculate whatever we are calculating
  calculateBias( args, myvals );
  if( actionInChain() ) {
    bool matout=false, matinp=getPntrToArgument(0)->getRank()==2;
#ifdef DNDEBUG
    if( matinp ) {
      for(unsigned i=1; i<getNumberOfArguments(); ++i) plumed_dbg_assert( getPntrToArgument(i)->getRank()==2 );
    }
#endif
    if( matinp ) {
      matout=getPntrToOutput(0)->getRank()==2;
#ifdef DNDEBUG
      if( matout ) {
        for(unsigned i=1; i<getNumberOfComponents(); ++i) plumed_dbg_assert( getPntrToOutput(i)->getRank()==2 );
      }
#endif
    }
    if( (matinp && matout) || !matinp ) {
      unsigned ostrn = getPntrToOutput(0)->getPositionInStream();
      for(unsigned i=0; i<distinct_arguments.size(); ++i) {
        unsigned istrn = (distinct_arguments[i].first->copyOutput(0))->getPositionInStream();
        for(unsigned k=0; k<myvals.getNumberActive(istrn); ++k) {
          unsigned kind = myvals.getActiveIndex(istrn,k);
          myvals.updateIndex( ostrn, arg_deriv_starts[i] + kind );
        }
      }
    } else if( myvals.inVectorCall() ) {
      // This requires further thought to make it future proof
      std::vector<unsigned> & indices( myvals.getIndices() );
      unsigned ostrn = getPntrToOutput(0)->getPositionInStream();
      myvals.clearActiveMembers( ostrn );
      for(unsigned i=0; i<myvals.getNumberOfIndices(); ++i) {
        myvals.updateIndex( ostrn, 3*indices[i]+0 ); myvals.updateIndex( ostrn, 3*indices[i]+1 ); myvals.updateIndex( ostrn, 3*indices[i]+2 );
      }
      unsigned nbase = nderivatives - 9;
      for(unsigned i=0; i<9; ++i) myvals.updateIndex( ostrn, nbase + i );
    }
  } else {
    unsigned ostrn = getPntrToOutput(0)->getPositionInStream();
    for(unsigned i=0; i<getNumberOfScalarArguments(); ++i) myvals.updateIndex( ostrn, i );
  }
}

void MultiBias::apply() {
  // Add forces due to this bias
  if( onStep()) {
    double gstr = static_cast<double>(getStride());
    getPntrToComponent(0)->addForce(0, -1.0*gstr );
  }
  // And add all the forces from here and elsewhere
  std::fill(forcesToApply.begin(),forcesToApply.end(),0); unsigned ss=0;
  if( getForcesFromValues( forcesToApply ) ) setForcesOnArguments( 0, forcesToApply, ss );
}

}
}


