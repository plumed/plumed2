/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2014 The plumed team
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
#ifndef __PLUMED_analysis_DimensionalityReductionBase_h
#define __PLUMED_analysis_DimensionalityReductionBase_h

#include "AnalysisWithLandmarks.h"

namespace PLMD {

class PointWiseMapping;

namespace analysis {

class DimensionalityReductionBase : public AnalysisWithLandmarks {
private:
  unsigned nlow;
  double cgtol;
  std::string ofilename;
  std::string efilename;
  PointWiseMapping* myembedding;
  std::vector<double> fframes;
  Matrix<double> targetDisimilarities;
protected:
  const Matrix<double>& getTargets() const ;
public:
  static void registerKeywords( Keywords& keys );
  DimensionalityReductionBase( const ActionOptions& );
  ~DimensionalityReductionBase();
  void analyzeLandmarks();
  virtual std::string getAlgorithmName() const=0;
  virtual void calculateAllDistances( PointWiseMapping* mymap, Matrix<double>& targets )=0;
  virtual void generateProjections( PointWiseMapping* mymap )=0;
  virtual double transformHD( const double& val, double& df ) const=0;
  virtual double transformLD( const double& val, double& df ) const=0;
  void findClosestPoint( const int& ii, ReferenceConfiguration* myref, std::vector<double>& pp );
  double calculateStress( const std::vector<double>& pp, std::vector<double>& der );
  void analyzeAllData();
  unsigned getLowDimensionSize() const ;
  void getPropertyNames( std::vector<std::string>& dimnames );
  void getProjectedPoint( const unsigned& idata, std::vector<double>& pp );
};

inline
const Matrix<double>& DimensionalityReductionBase::getTargets() const {
  return targetDisimilarities;
}

inline
unsigned DimensionalityReductionBase::getLowDimensionSize() const {
  return nlow;
}

}
}
#endif
