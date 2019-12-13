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
#include "DimensionalityReductionBase.h"
#include "core/ActionRegister.h"
#include "SMACOF.h"

//+PLUMEDOC DIMRED SMACOF_MDS
/*
Optimize the multidimensional scaling stress function using the SMACOF algorithm.

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace dimred {

class SmacofMDS : public DimensionalityReductionBase {
private:
  unsigned maxloops;
  double tol;
public:
  static void registerKeywords( Keywords& keys );
  explicit SmacofMDS( const ActionOptions& );
  void calculateProjections( const Matrix<double>&, Matrix<double>& ) override;
};

PLUMED_REGISTER_ACTION(SmacofMDS,"SMACOF_MDS")

void SmacofMDS::registerKeywords( Keywords& keys ) {
  DimensionalityReductionBase::registerKeywords( keys );
  keys.remove("NLOW_DIM");
  keys.add("compulsory","SMACOF_TOL","1E-4","tolerance for the SMACOF optimization algorithm");
  keys.add("compulsory","SMACOF_MAXCYC","1000","maximum number of optimization cycles for SMACOF algorithm");
}

SmacofMDS::SmacofMDS( const ActionOptions& ao):
  Action(ao),
  DimensionalityReductionBase(ao)
{
  if( !dimredbase ) error("SMACOF must be initialized using output from dimensionality reduction object");

  parse("SMACOF_TOL",tol); parse("SMACOF_MAXCYC",maxloops);
  log.printf("  running smacof to convergence at %f or for a maximum of %u steps \n",tol,maxloops);
}

void SmacofMDS::calculateProjections( const Matrix<double>& targets, Matrix<double>& projections ) {
  // Take the square root of all the distances and the weights
  Matrix<double> weights( targets.nrows(), targets.ncols() );
  Matrix<double> distances( targets.nrows(), targets.ncols() );
  for(unsigned i=1; i<distances.ncols(); ++i) {
    for(unsigned j=0; j<i; ++j) {
      distances(i,j)=distances(j,i)=sqrt( targets(i,j) );
      weights(i,j)=weights(j,i)=getWeight(i)*getWeight(j);
    }
  }
  // And run SMACOF
  SMACOF::run( weights, targets, tol, maxloops, projections );
}

}
}
