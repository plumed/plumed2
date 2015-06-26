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
#ifndef __PLUMED_analysis_LandmarkSelectionBase_h
#define __PLUMED_analysis_LandmarkSelectionBase_h

#include "AnalysisWithAnalysableOutput.h"

namespace PLMD {
namespace analysis {

class LandmarkSelectionBase : public AnalysisWithAnalysableOutput {
private:
/// The number of landmarks we are selecting
  unsigned nlandmarks;
/// The indices of the landmarks in the original data set
  std::vector<unsigned> landmark_indices;
/// How do we treat weights
  bool novoronoi, noweights;
protected:
/// Return the number of landmarks
  unsigned getNumberOfLandmarks() const ;
/// Transfer frame i in the underlying action to the object we are going to analyze
  void selectFrame( const unsigned& );
public:
  static void registerKeywords( Keywords& keys );  
  LandmarkSelectionBase( const ActionOptions& ao );
  void performAnalysis();
  virtual void selectLandmarks()=0;
  ReferenceConfiguration* getOutputConfiguration( const unsigned& idata );
  double getDistanceBetweenFrames( const unsigned& iframe, const unsigned& jframe, const bool& sq );
  void performTask(){ plumed_error(); }
  double getOutputDissimilarity( const unsigned& idata, const unsigned& jdata );
};

inline
unsigned LandmarkSelectionBase::getNumberOfLandmarks() const {
  return nlandmarks;
}

inline
double LandmarkSelectionBase::getDistanceBetweenFrames( const unsigned& iframe, const unsigned& jframe, const bool& sq ){
  return Analysis::getDistanceBetweenFrames( landmark_indices[iframe], landmark_indices[jframe], true );
}

}
}
#endif
