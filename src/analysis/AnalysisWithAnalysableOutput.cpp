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
#include "DissimilarityMatrixBase.h"
#include "reference/ReferenceConfiguration.h"
#include "reference/MetricRegister.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"

namespace PLMD {
namespace analysis {

void AnalysisWithAnalysableOutput::registerKeywords( Keywords& keys ){
  Analysis::registerKeywords( keys );
  keys.add("optional","DISSIMILARITIES","the label of the action that calculates the dissimilarities");
}

AnalysisWithAnalysableOutput::AnalysisWithAnalysableOutput( const ActionOptions& ao ):
Action(ao),
Analysis(ao),
mydissims(NULL),
myinput(NULL)
{
  std::string dlab; parse("DISSIMILARITIES",dlab);
  if( dlab.length()>0 ){
    log.printf("  calculating dissimilarities using Action with label %s \n",dlab.c_str() );
    mydissims = plumed.getActionSet().selectWithLabel<DissimilarityMatrixBase*>(dlab);
    if(!mydissims) error("action labelled " +  dlab + " does not exist or is not of DissimilarityMatrixBase type");
    addDependency( mydissims );
  } else if( dimredstash ){
    myinput = dynamic_cast<AnalysisWithAnalysableOutput*>( dimredstash ); 
  }
}

void AnalysisWithAnalysableOutput::setNumberOfOutputPoints( const unsigned& n ){
  noutput_points=n;
}

void AnalysisWithAnalysableOutput::setOutputWeights( const std::vector<double>& wwwin ){
  plumed_dbg_assert( wwwin.size()==noutput_points );
  if( oweights.size()!=wwwin.size() ) oweights.resize( wwwin.size() );
  for(unsigned i=0;i<wwwin.size();++i) oweights[i]=wwwin[i];
}

void AnalysisWithAnalysableOutput::getOutputForPoint( const unsigned& idata, std::vector<double>& point ){
  error("cannot construct histograms using this kind of object");
}

double AnalysisWithAnalysableOutput::getDissimilarity( const unsigned& idata, const unsigned& jdata ){
  if( mydissims ) return mydissims->getDissimilarity( idata, jdata );
  if( dimredstash ) return myinput->getOutputDissimilarity( idata, jdata );
  plumed_error(); return 1.0;  
}

bool AnalysisWithAnalysableOutput::dissimilaritiesWereSet(){
  return ( mydissims || dimredstash );
}

}
}
