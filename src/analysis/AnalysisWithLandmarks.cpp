/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2016 The plumed team
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
#include "AnalysisWithLandmarks.h"
#include "LandmarkRegister.h"
#include "LandmarkSelectionBase.h"

//+PLUMEDOC INTERNAL landmarkselection
/*
This is currently a filler page.  

Just use LANDMARKS=ALL.  More complex versions will appear in later versions.

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace analysis {

void AnalysisWithLandmarks::registerKeywords( Keywords& keys ){
  Analysis::registerKeywords( keys );
  keys.add("compulsory","LANDMARKS","ALL","only use a subset of the data that was collected. "
                                          "For more information on the landmark selection algorithms that are available in "
                                          "plumed see \\ref landmarkselection.");
}

AnalysisWithLandmarks::AnalysisWithLandmarks( const ActionOptions& ao):
Action(ao),
Analysis(ao),
data_to_analyze(NULL)
{
   std::string linput; parse("LANDMARKS",linput);
   std::vector<std::string> words=Tools::getWords(linput); 
   landmarkSelector=landmarkRegister().create( LandmarkSelectionOptions(words,this) );
   log.printf("  %s\n", landmarkSelector->description().c_str() );
}

AnalysisWithLandmarks::~AnalysisWithLandmarks(){
   delete landmarkSelector;
}

void AnalysisWithLandmarks::setDataToAnalyze( MultiReferenceBase* mydata ){
   data_to_analyze=mydata;
}

unsigned AnalysisWithLandmarks::getNumberOfLandmarks() const {
  return landmarkSelector->getNumberOfLandmarks();
}

void AnalysisWithLandmarks::performAnalysis(){
  plumed_assert( data_to_analyze );
  landmarkSelector->selectLandmarks( data_to_analyze );
  analyzeLandmarks();
}

void AnalysisWithLandmarks::performTask( const unsigned& taskIndex, const unsigned& current, MultiValue& myvals ) const {
  plumed_merror("Should not be here");
}

}
}
