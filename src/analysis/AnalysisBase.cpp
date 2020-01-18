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
#include "AnalysisBase.h"
#include "ReadAnalysisFrames.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"

namespace PLMD {
namespace analysis {

void AnalysisBase::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys ); ActionPilot::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys ); ActionAtomistic::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys ); keys.remove("NUMERICAL_DERIVATIVES");
  ActionWithVessel::registerKeywords( keys ); keys.remove("TOL"); keys.reset_style("TIMINGS","hidden"); keys.isAnalysis();
  keys.add("atoms-2","USE_OUTPUT_DATA_FROM","use the output of the analysis performed by this object as input to your new analysis object");
}

AnalysisBase::AnalysisBase(const ActionOptions&ao):
  Action(ao),
  ActionPilot(ao),
  ActionWithValue(ao),
  ActionAtomistic(ao),
  ActionWithArguments(ao),
  ActionWithVessel(ao),
  my_input_data(NULL)
{
  // We have an if statement here so that this doesn't break with READ_DISSIMILARITIES
  if( keywords.exists("USE_OUTPUT_DATA_FROM") ) {
    std::string datastr; parse("USE_OUTPUT_DATA_FROM",datastr);
    if( keywords.style("USE_OUTPUT_DATA_FROM","atoms") && datastr.length()==0 ) error("input analysis action was not specified use USE_OUTPUT_DATA_FROM");
    if( datastr.length()>0 ) {
      my_input_data=plumed.getActionSet().selectWithLabel<AnalysisBase*>( datastr );
      log.printf("  performing analysis on output from %s \n",datastr.c_str() );
      if( !my_input_data ) error("could not find analysis action named " + datastr );
      addDependency( my_input_data );
    }
  }
}

std::vector<std::string> AnalysisBase::getArgumentNames() {
  std::vector<Value*> arg_p( getArgumentList() );
  std::vector<std::string> argn( arg_p.size() );
  for(unsigned i=0; i<arg_p.size(); ++i) {
    plumed_assert( i<argn.size() && i<arg_p.size() );
    argn[i]=arg_p[i]->getName();
  }
  return argn;
}

void AnalysisBase::update() {
  if( getStep()==0 || ( getStride()>0 && !onStep() ) ) return;
  // And do the analysis
  performAnalysis();
}

void AnalysisBase::runFinalJobs() {
  if( getStride()>0 ) return;
  performAnalysis();
}

}
}
