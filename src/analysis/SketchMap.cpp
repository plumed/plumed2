/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

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
#include "ClassicalScaling.h"
#include "SMACOF.h"
#include "reference/PointWiseMapping.h"
#include "core/ActionRegister.h"
#include "tools/SwitchingFunction.h"

//+PLUMEDOC ANALYSIS SKETCHMAP
/*
Perform a dimensionality reduction using the sketch-map algorithm

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace analysis {

class SketchMap : public DimensionalityReductionBase {
private:
  bool nosmacof;
  double smactol, smaptol, regulariser;
  SwitchingFunction lowdf, highdf;
  double recalculateWeights( const Matrix<double>& Distances, const Matrix<double>& F, PointWiseMapping* mymap, Matrix<double>& Weights );
public:
  static void registerKeywords( Keywords& keys );
  SketchMap( const ActionOptions& ao );
  std::string getAlgorithmName() const { return "sketch-map"; }
  void calculateAllDistances( PointWiseMapping* mymap, Matrix<double>& targets );
  void generateProjections( PointWiseMapping* mymap );
  double transformHD( const double& val, double& df ) const ;
  double transformLD( const double& val, double& df ) const ;
};

PLUMED_REGISTER_ACTION(SketchMap,"SKETCHMAP")

void SketchMap::registerKeywords( Keywords& keys ){
  DimensionalityReductionBase::registerKeywords( keys );
  keys.add("compulsory","HIGH_DIM_FUNCTION","the parameters of the switching function in the high dimensional space");
  keys.add("compulsory","LOW_DIM_FUNCTION","the parameters of the switching function in the low dimensional space");
  keys.add("compulsory","SMACOF_TOL","1E-4","the tolerance for each SMACOF cycle");
  keys.add("compulsory","SMAP_TOL","1E-4","the tolerance for sketch-map");
  keys.add("compulsory","REGULARISE_PARAM","0.001","this is used to ensure that we don't divide by zero when updating weights");
}

SketchMap::SketchMap( const ActionOptions& ao ):
Action(ao),
DimensionalityReductionBase(ao)
{
  // Read in the switching functions
  std::string linput,hinput, errors;
  parse("HIGH_DIM_FUNCTION",hinput);
  highdf.set(hinput,errors);
  if(errors.length()>0) error(errors);
  parse("LOW_DIM_FUNCTION",linput);
  lowdf.set(hinput,errors);
  if(errors.length()>0) error(errors);

  // Read tolerances
  parse("SMACOF_TOL",smactol);
  parse("SMAP_TOL",smaptol);
  parse("REGULARISE_PARAM",regulariser);
}

void SketchMap::calculateAllDistances( PointWiseMapping* mymap, Matrix<double>& targets ){
  // Calculate matrix of dissimilarities (High dimensional space) 
  mymap->calculateAllDistances( getPbc(), getArguments(), comm, mymap->modifyDmat(), false );
  double dr; unsigned M = mymap->getNumberOfReferenceFrames(); 
  for(unsigned i=1; i<M; ++i){
      for(unsigned j=0; j<i; ++j) {
          targets(i,j) = targets(j,i) = 1.0 - highdf.calculate( mymap->modifyDmat()(i,j), dr ); // high dim space
      }
  }
}


void SketchMap::generateProjections( PointWiseMapping* mymap ){
  Matrix<double> Distances( mymap->modifyDmat() ); //making a copy
  // Calculates the first guess of projections in LD space 
  ClassicalScaling::run( mymap );
  
  // Calculate the value of sigma and the weights
  unsigned M = mymap->getNumberOfReferenceFrames();
  Matrix<double> Weights(M,M); double filt = recalculateWeights( Distances, getTargets(), mymap, Weights );

  unsigned MAXSTEPS=100; double newsig;
  for(unsigned i=0;i<MAXSTEPS;++i){
      // Run the smacof algorithm
      SMACOF::run( Weights, mymap, smactol );
      // Recalculate weights matrix and sigma
      newsig = recalculateWeights( Distances, getTargets(), mymap, Weights );
      printf("HELLO GARETH AND RACHEL %d %f %f %f \n",i, newsig, filt, fabs( newsig - filt ) );
      // Test whether or not the algorithm has converged
      if( fabs( newsig - filt )<smaptol ) break;
      // Make initial sigma into new sigma so that the value of new sigma is used every time so that the error can be reduced
      filt=newsig;
  } 
}

double SketchMap::recalculateWeights( const Matrix<double>& Distances, const Matrix<double>& F, PointWiseMapping* mymap, Matrix<double>& Weights ){
  double filt=0, totalWeight=0.;; double dr;
  for(unsigned i=1; i<Weights.nrows(); ++i){
      for(unsigned j=0; j<i; ++j){
          double tempd=0;
          for(unsigned k=0;k<mymap->getNumberOfProperties();++k){
             double tmp = mymap->getProjectionCoordinate( i, k ) - mymap->getProjectionCoordinate( j, k );
             tempd += tmp*tmp;
          }
          double ninj=mymap->getWeight(i)*mymap->getWeight(j);
          totalWeight += ninj;

          double dij=sqrt(tempd);
          double fij = 1.0 - lowdf.calculate( dij, dr );
          double filter=F(i,j)-fij;
          double diff=Distances(i,j) - dij;
          if( fabs(diff)<regulariser ) Weights(i,j)=Weights(j,i)=0.0;
          else Weights(i,j)=Weights(j,i) = ( -ninj*filter*dij*dr ) / diff;
          filt += ninj*filter*filter;
      }
  }
  return filt / totalWeight;
}

double SketchMap::transformHD( const double& val, double& df ) const { 
  double vv=1.0 - highdf.calculate( val, df );
  df=-df; return vv;
}
double SketchMap::transformLD( const double& val, double& df ) const {
  double vv=1.0 - lowdf.calculate( val, df );
  df=-df; return vv;
}

}
}
