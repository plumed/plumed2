/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2019 The plumed team
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
#include "core/ActionRegister.h"
#include "SketchMapBase.h"
#include "SMACOF.h"

//+PLUMEDOC DIMRED SKETCHMAP_SMACOF
/*
Optimize the sketch-map stress function using the SMACOF algorithm.

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace dimred {

class SketchMapSmacof : public SketchMapBase {
private:
  unsigned max_smap, maxiter;
  double smap_tol, iter_tol, regulariser;
  double recalculateWeights( const Matrix<double>& projections, Matrix<double>& weights );
public:
  static void registerKeywords( Keywords& keys );
  explicit SketchMapSmacof( const ActionOptions& ao );
  void minimise( Matrix<double>& ) override;
};

PLUMED_REGISTER_ACTION(SketchMapSmacof,"SKETCHMAP_SMACOF")

void SketchMapSmacof::registerKeywords( Keywords& keys ) {
  SketchMapBase::registerKeywords( keys );
  keys.add("compulsory","SMACOF_TOL","1E-4","the tolerance for each SMACOF cycle");
  keys.add("compulsory","SMACOF_MAXCYC","1000","maximum number of optimization cycles for SMACOF algorithm");
  keys.add("compulsory","SMAP_TOL","1E-4","the tolerance for sketch-map");
  keys.add("compulsory","SMAP_MAXCYC","100","maximum number of optimization cycles for iterative sketch-map algorithm");
  keys.add("compulsory","REGULARISE_PARAM","0.001","this is used to ensure that we don't divide by zero when updating weights");
}

SketchMapSmacof::SketchMapSmacof( const ActionOptions& ao ):
  Action(ao),
  SketchMapBase(ao)
{
  parse("REGULARISE_PARAM",regulariser);
  parse("SMACOF_MAXCYC",max_smap); parse("SMAP_MAXCYC",maxiter);
  parse("SMACOF_TOL",smap_tol); parse("SMAP_TOL",iter_tol);
}

void SketchMapSmacof::minimise( Matrix<double>& projections ) {
  Matrix<double> weights( distances.nrows(), distances.ncols() ); weights=0.;
  double filt = recalculateWeights( projections, weights );

  for(unsigned i=0; i<maxiter; ++i) {
    SMACOF::run( weights, distances, smap_tol, max_smap, projections );
    // Recalculate weights matrix and sigma
    double newsig = recalculateWeights( projections, weights );
    // Test whether or not the algorithm has converged
    if( fabs( newsig - filt )<iter_tol ) break;
    // Make initial sigma into new sigma so that the value of new sigma is used every time so that the error can be reduced
    filt=newsig;
  }
}

double SketchMapSmacof::recalculateWeights( const Matrix<double>& projections, Matrix<double>& weights ) {
  double filt=0, totalWeight=0.;; double dr;
  for(unsigned i=1; i<weights.nrows(); ++i) {
    for(unsigned j=0; j<i; ++j) {
      double ninj=getWeight(i)*getWeight(j); totalWeight += ninj;

      double tempd=0;
      for(unsigned k=0; k<projections.ncols(); ++k) {
        double tmp = projections(i,k) - projections(j,k);
        tempd += tmp*tmp;
      }
      double dij=sqrt(tempd);

      double fij = transformLowDimensionalDistance( dij, dr );
      double filter=transformed(i,j)-fij;
      double diff=distances(i,j) - dij;

      if( fabs(diff)<regulariser ) weights(i,j)=weights(j,i)=0.0;
      else weights(i,j)=weights(j,i) = ninj*( (1-mixparam)*( filter*dr )/diff + mixparam );
      filt += ninj*( (1-mixparam)*filter*filter + mixparam*diff*diff );
    }
  }
  return filt / totalWeight;
}

}
}
