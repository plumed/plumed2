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
#include "AnalysisWithAnalysableOutput.h"
#include "reference/ReferenceConfiguration.h"
#include "reference/MetricRegister.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"

namespace PLMD {
namespace analysis {

void AnalysisWithAnalysableOutput::registerKeywords( Keywords& keys ){
  Analysis::registerKeywords( keys );
}

AnalysisWithAnalysableOutput::AnalysisWithAnalysableOutput( const ActionOptions& ao ):
Action(ao),
Analysis(ao)
{
}

void AnalysisWithAnalysableOutput::setOutputWeights( const std::vector<double>& wwwin ){
  plumed_dbg_assert( wwwin.size()==noutput_points );
  if( oweights.size()!=wwwin.size() ) oweights.resize( wwwin.size() );
  for(unsigned i=0;i<wwwin.size();++i) oweights[i]=wwwin[i];
}

void AnalysisWithAnalysableOutput::getOutputForPoint( const unsigned& idata, std::vector<double>& point ){
  error("cannot construct histograms using this kind of object");
}

}
}
