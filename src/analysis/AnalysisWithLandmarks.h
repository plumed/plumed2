/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2017 The plumed team
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
#ifndef __PLUMED_analysis_AnalysisWithLandmarks_h
#define __PLUMED_analysis_AnalysisWithLandmarks_h

#include "Analysis.h"

namespace PLMD {

class MultiReferenceBase;

namespace analysis {

class LandmarkSelectionBase;

class AnalysisWithLandmarks : public Analysis {
  friend class LandmarkSelectionBase;
  friend class CopyAllFrames;
private:
/// This object selects landmarks from the data
  LandmarkSelectionBase* landmarkSelector;
/// A pointer to the data we are analyzing
  MultiReferenceBase* data_to_analyze;
protected:
/// Set the data that needs to be analyzed
  void setDataToAnalyze( MultiReferenceBase* mydata );
/// Return the number of landmarks we are selecting
  unsigned getNumberOfLandmarks() const ;
public:
  static void registerKeywords( Keywords& keys );
  explicit AnalysisWithLandmarks( const ActionOptions& );
  ~AnalysisWithLandmarks();
/// Do the analysis
  void performAnalysis();
  virtual void analyzeLandmarks()=0;
/// This does nothing
  void performTask( const unsigned&, const unsigned&, MultiValue& ) const ;
};

}
}
#endif
