/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013,2014 The plumed team
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
#include "ReadAnalysisFrames.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/ActionRegister.h"

//+PLUMEDOC ANALYSIS READ_ANALYSIS_FRAMES
/* 
This allows you to convert a trajectory and a dissimilarity matrix into a dissimilarity object

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace analysis {

PLUMED_REGISTER_ACTION(ReadAnalysisFrames,"READ_ANALYSIS_FRAMES")

void ReadAnalysisFrames::registerKeywords( Keywords& keys ){
  AnalysisWithDataCollection::registerKeywords( keys ); 
  keys.reset_style("USE_ALL_DATA","hidden"); keys.remove("RUN"); keys.remove("REWEIGHT_BIAS"); keys.remove("REWEIGHT_TEMP"); 
  keys.remove("TEMP"); keys.remove("WRITE_CHECKPOINT"); keys.remove("NOMEMORY"); keys.remove("RESTART"); 
  keys.remove("UPDATE_FROM"); keys.remove("UPDATE_UNTIL"); keys.remove("USE_OUTPUT_DATA_FROM");
  keys.remove("ARG"); keys.remove("SERIAL"); keys.remove("REUSE_INPUT_DATA_FROM");
}

ReadAnalysisFrames::ReadAnalysisFrames( const ActionOptions& ao ):
Action(ao),
AnalysisWithDataCollection(ao)
{
  if( plumed.getActionSet().size()!=0 ) error("read analysis frames command must be at top of input file");
}

}
}
