/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2020 The plumed team
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
#ifndef __PLUMED_dimred_SketchMapBase_h
#define __PLUMED_dimred_SketchMapBase_h

#include "tools/SwitchingFunction.h"
#include "DimensionalityReductionBase.h"

namespace PLMD {
namespace dimred {

class SketchMapBase : public DimensionalityReductionBase {
private:
/// To save us retyping switching functions many times the code will reuse the
/// ones from previous sketch-map objects
  bool reuse_hd, reuse_ld;
/// Was previous action in chain a sketch-map action
  SketchMapBase* smapbase;
/// Switching functions for low and high dimensional space
  SwitchingFunction lowdf, highdf;
/// This is used within calculate stress to hold the target distances and the
/// target values for the high dimensional switching function
  std::vector<double> dtargets, ftargets, pweights;
/// Stress normalization (sum_ij w_i w_j)
  double normw;
protected:
/// This holds the target distances and target transformed distances
  Matrix<double> distances, transformed;
/// The fraction of pure distances to mix in when optimising
  double mixparam;
public:
  static void registerKeywords( Keywords& keys );
  explicit SketchMapBase( const ActionOptions& );
/// This starts the process of calculating the projections
  void calculateProjections( const Matrix<double>&, Matrix<double>& ) override;
/// This finishes the process of calculating the prjections
  virtual void minimise( Matrix<double>& )=0;
/// Apply the low dimensional switching function to the value val
  double transformLowDimensionalDistance( const double& val, double& df ) const ;
/// Apply the high dimensional switching function to the value val
  double transformHighDimensionalDistance( const double& val, double& df ) const ;
/// Set the target distance and from it calculate the target value for the switching function
/// This target vector is used when we use calculateStress when finding the projections of individual points.
/// For example this function is used in PLMD::dimred::ProjectOutOfSample
  void setTargetDistance( const unsigned& idata, const double& dist ) override;
/// Calculate the pointwise stress on one point when it is located at p.
/// This function makes use of the distance data in dtargets and ftargets
/// It is used in PLMD::dimred::ProjectOutOfSample and in pointwise optimisation
  double calculateStress( const std::vector<double>& p, std::vector<double>& d ) override;
/// Calculate the total stress when the projections are placed at point p.  Notice
/// this is a vectorized version of the matrix of projections
  double calculateFullStress( const std::vector<double>& p, std::vector<double>& d );
};

inline
double SketchMapBase::transformLowDimensionalDistance( const double& val, double& df ) const {
  if( reuse_ld ) return smapbase->transformLowDimensionalDistance( val, df );
  double ans=1.0 - lowdf.calculate( val, df ); df*=-val; return ans;
}

inline
double SketchMapBase::transformHighDimensionalDistance( const double& val, double& df ) const {
  if( reuse_hd ) return smapbase->transformHighDimensionalDistance( val, df );
  double ans=1.0 - highdf.calculate( val, df ); df*=-val; return ans;
}

inline
void SketchMapBase::setTargetDistance( const unsigned& idata, const double& dist ) {
  double df; dtargets[idata]=dist; ftargets[idata]=transformHighDimensionalDistance( dist, df );
}

}
}
#endif
