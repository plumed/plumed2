/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2015 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)
  
   See http://www.plumed-code.org for more information.

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

void AnalysisBase::registerKeywords( Keywords& keys ){
  Action::registerKeywords( keys );
  ActionPilot::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys );
  ActionWithVessel::registerKeywords( keys ); keys.remove("TOL"); keys.reset_style("TIMINGS","hidden"); keys.isAnalysis();  
  keys.add("atoms-2","USE_OUTPUT_DATA_FROM","use the ouput of the analysis performed by this object as input to your new analysis object");
}

AnalysisBase::AnalysisBase(const ActionOptions&ao):
Action(ao),
ActionPilot(ao),
ActionAtomistic(ao),        
ActionWithArguments(ao),
ActionWithVessel(ao),
mydata(NULL)
{
  // We have an if statement here so that this doesn't break with READ_DISSIMILARITIES
  std::string datastr; if( keywords.exists("USE_OUTPUT_DATA_FROM") ) parse("USE_OUTPUT_DATA_FROM",datastr);
  if( keywords.exists("USE_OUTPUT_DATA_FROM") && 
      datastr.length()==0 && 
      !keywords.exists("REUSE_INPUT_DATA_FROM") ) error("input analysis action was not specified use USE_OUTPUT_DATA_FROM");
  if( datastr.length()>0 ){
      mydata=plumed.getActionSet().selectWithLabel<AnalysisBase*>( datastr );
      ReadAnalysisFrames* checkt = dynamic_cast<ReadAnalysisFrames*>( mydata );
      if( checkt ) error("READ_ANALYSIS_FRAMES should only be used in input to the FRAMES keyword"); 
      log.printf("  performing analysis on output from %s \n",datastr.c_str() );
      if( !mydata ) error("could not find analysis action named " + datastr );
      freq=mydata->freq; use_all_data=mydata->use_all_data;
      if( !use_all_data ) setStride( freq );
  }
}

void AnalysisBase::update(){
  // Do nothing if we are just analysing at the end of the calculation
  if( use_all_data ) return ;
  // Check that we are on step
  plumed_dbg_assert( getStep()%freq==0 );
  // And do the analysis
  if( getStep()>0 ) performAnalysis();
}

void AnalysisBase::runFinalJobs(){
  // Nothing to do if we are not analysing all the data in the trajectory
  if( !use_all_data ) return;
  // Erm ... user has done something weird
  if( getNumberOfDataPoints()==0 ) error("no data is available for analysis");
  // And do the analysis
  if( use_all_data ) performAnalysis();
}

}
}
