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
#include "FunctionBase.h"

namespace PLMD {
namespace function {

void FunctionBase::registerKeywords(Keywords& keys) {
  Action::registerKeywords(keys); ActionWithValue::registerKeywords(keys);
  ActionWithArguments::registerKeywords(keys); keys.use("ARG");
  keys.reserve("compulsory","PERIODIC","if the output of your function is periodic then you should specify the periodicity of the function.  If the output is not periodic you must state this using PERIODIC=NO");
}

FunctionBase::FunctionBase(const ActionOptions&ao):
Action(ao),
ActionWithValue(ao),
ActionWithArguments(ao),
firststep(true)
{
}

void FunctionBase::buildCurrentTaskList( bool& forceAllTasks, std::vector<std::string>& actionsThatSelectTasks, std::vector<unsigned>& tflags ) {
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

void FunctionBase::evaluateAllFunctions() {
/// NEEDS FIGURING OUT GARETH
  if( firststep ) { setupOnFirstStep(); firststep=false; }
  //   std::vector<unsigned> shape( getShape() );
  //   unsigned ival = getPntrToOutput(0)->getNumberOfValues();
  //   getPntrToOutput(0)->setShape( shape ); firststep=false;
  //   if( ival<getPntrToOutput(0)->getNumberOfValues() ) {
  //     for(unsigned j=ival; j<getPntrToOutput(0)->getNumberOfValues(); ++j) addTaskToList(j);
  //   }
  runAllTasks();
}

void FunctionBase::calculate() {
  // Everything is done elsewhere
  if( actionInChain() ) return;
  // This is done if we are calculating a function of multiple cvs
  evaluateAllFunctions();
}

void FunctionBase::update() {
  if( skipUpdate() || actionInChain() ) return;
  plumed_dbg_assert( !actionInChain() );
  if( getFullNumberOfTasks()>0 ) evaluateAllFunctions();
}

void FunctionBase::runFinalJobs() {
  if( skipUpdate() || actionInChain() ) return;
  resizeForFinalTasks(); evaluateAllFunctions();
}

void FunctionBase::applyNonGrid() {
  if( doNotCalculateDerivatives() ) return;

  if( forcesToApply.size()!=getNumberOfDerivatives() ) forcesToApply.resize( getNumberOfDerivatives() );

  std::fill(forcesToApply.begin(),forcesToApply.end(),0); unsigned ss=0;
  if( getForcesFromValues( forcesToApply ) ) setForcesOnArguments( 0, forcesToApply, ss );
}

}
}
