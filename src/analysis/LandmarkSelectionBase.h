/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2019 The plumed team
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
#ifndef __PLUMED_analysis_LandmarkSelectionBase_h
#define __PLUMED_analysis_LandmarkSelectionBase_h

#include "AnalysisBase.h"

namespace PLMD {
namespace analysis {

class LandmarkSelectionBase : public AnalysisBase {
  friend class ReselectLandmarks;
private:
/// The number of landmarks we are selecting
  unsigned nlandmarks;
/// The weights of the landmark points
  std::vector<double> lweights;
/// The indices of the landmarks in the original data set
  std::vector<unsigned> landmark_indices;
/// How do we treat weights
  bool novoronoi, noweights;
protected:
/// Transfer frame i in the underlying action to the object we are going to analyze
  void selectFrame( const unsigned& );
/// Do a voronoi analysis
  void voronoiAnalysis( const std::vector<unsigned>& myindices, std::vector<double>& lweights, std::vector<unsigned>& assignments ) const ;
public:
  static void registerKeywords( Keywords& keys );
  explicit LandmarkSelectionBase( const ActionOptions& ao );
/// Return the number of data points
  unsigned getNumberOfDataPoints() const override;
/// Return the index of the data point in the base class
  unsigned getDataPointIndexInBase( const unsigned& idata ) const override;
/// Get the weight
  double getWeight( const unsigned& idata ) override;
/// Get a reference configuration
  DataCollectionObject& getStoredData( const unsigned& idat, const bool& calcdist ) override;
/// Select landmark configurations
  void performAnalysis() override;
  virtual void selectLandmarks()=0;
/// Get the squared dissimilarity between two reference configurations
  double getDissimilarity( const unsigned& i, const unsigned& j ) override;
/// This does nothing - it just ensures the final class is not abstract
  void performTask( const unsigned&, const unsigned&, MultiValue& ) const override { plumed_error(); }
};

inline
unsigned LandmarkSelectionBase::getNumberOfDataPoints() const {
  return nlandmarks;
}

inline
unsigned LandmarkSelectionBase::getDataPointIndexInBase( const unsigned& idata ) const {
  return AnalysisBase::getDataPointIndexInBase( landmark_indices[idata] );
}

inline
double LandmarkSelectionBase::getWeight( const unsigned& idata ) {
  return lweights[idata];
}

inline
DataCollectionObject& LandmarkSelectionBase::getStoredData( const unsigned& idat, const bool& calcdist ) {
  return AnalysisBase::getStoredData( landmark_indices[idat], calcdist );
}

inline
double LandmarkSelectionBase::getDissimilarity( const unsigned& i, const unsigned& j ) {
  return AnalysisBase::getDissimilarity( landmark_indices[i], landmark_indices[j] );
}

}
}
#endif
