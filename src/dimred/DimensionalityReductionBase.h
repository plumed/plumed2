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
#ifndef __PLUMED_dimred_DimensionalityReductionBase_h
#define __PLUMED_dimred_DimensionalityReductionBase_h

#include "analysis/AnalysisWithAnalysableOutput.h"

namespace PLMD {
namespace dimred {

class DimensionalityReductionBase : public analysis::AnalysisWithAnalysableOutput {
friend class ProjectNonLandmarkPoints;
friend class SketchMapBase;
private:
  bool use_dimred_dissims;
  std::vector<double> dtargets;
  ReferenceConfiguration* mydata;
  Matrix<double> projections;
protected:
  unsigned nlow;
  DimensionalityReductionBase* dimredbase;
  double getInputDissimilarity( const unsigned& idata, const unsigned& jdata );
public:
  static void registerKeywords( Keywords& keys );
  DimensionalityReductionBase( const ActionOptions& );
  ~DimensionalityReductionBase();
  ReferenceConfiguration* getOutputConfiguration( const unsigned& idata );
  unsigned getNumberOfOutputPoints() const ;
  void getOutputForPoint( const unsigned& idata, std::vector<double>& point );
  double getOutputDissimilarity( const unsigned& idata, const unsigned& jdata );
  void performTask(){}
  void performAnalysis();
  virtual void calculateProjections( const Matrix<double>& , Matrix<double>& )=0;
  virtual void setTargetDistance( const unsigned& , const double& );
  virtual double calculateStress( const std::vector<double>& pp, std::vector<double>& der );
};

inline
unsigned DimensionalityReductionBase::getNumberOfOutputPoints() const {
  return getNumberOfDataPoints();
}

inline
void DimensionalityReductionBase::setTargetDistance( const unsigned& idata, const double& dist ){
  dtargets[idata]=dist;
}

}
}
#endif
