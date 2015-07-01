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
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "tools/Random.h"
#include "reference/MetricRegister.h"
#include "tools/ConjugateGradient.h"
#include "analysis/AnalysisWithAnalysableOutput.h"
#include "DimensionalityReductionBase.h"

namespace PLMD {
namespace dimred {

class ProjectNonLandmarkPoints : public analysis::AnalysisWithAnalysableOutput {
private:
  double cgtol;
  unsigned nlow;
  ReferenceConfiguration* mydata; 
  DimensionalityReductionBase* mybase;
  void generateProjection( const unsigned& idata, std::vector<double>& point );
public:
  static void registerKeywords( Keywords& keys );
  ProjectNonLandmarkPoints( const ActionOptions& ao );
  ~ProjectNonLandmarkPoints();
  ReferenceConfiguration* getOutputConfiguration( const unsigned& idata );
  unsigned getNumberOfOutputPoints() const ;
  void getOutputForPoint( const unsigned& idata, std::vector<double>& point );
  double getOutputDissimilarity( const unsigned& idata, const unsigned& jdata );
  void performTask(){}
  void performAnalysis();
  double calculateStress( const std::vector<double>& pp, std::vector<double>& der );
};

PLUMED_REGISTER_ACTION(ProjectNonLandmarkPoints,"PROJECT_ALL_ANALYSIS_DATA")

void ProjectNonLandmarkPoints::registerKeywords( Keywords& keys ){
  analysis::AnalysisWithAnalysableOutput::registerKeywords( keys );
  keys.add("compulsory","PROJECTION","the projection that you wish to generate out-of-sample projections with");
  keys.add("compulsory","CGTOL","1E-6","the tolerance for the conjugate gradient optimisation");
}

ProjectNonLandmarkPoints::ProjectNonLandmarkPoints( const ActionOptions& ao ):
Action(ao),
analysis::AnalysisWithAnalysableOutput(ao),
mybase(NULL)
{
  std::string myproj; parse("PROJECTION",myproj);
  mybase = plumed.getActionSet().selectWithLabel<DimensionalityReductionBase*>( myproj );
  if( !mybase ) error("could not find projection of data named " + myproj ); 
  nlow = mybase->nlow;

  if( mybase->getBaseDataLabel()!=getBaseDataLabel() ) error("mismatch between base data labels"); 

  log.printf("  generating out-of-sample projections using projection with label %s \n",myproj.c_str() );
  parse("CGTOL",cgtol);

  ReferenceConfigurationOptions("EUCLIDEAN");
  mydata=metricRegister().create<ReferenceConfiguration>("EUCLIDEAN");
  std::vector<std::string> dimnames(nlow); std::string num;
  for(unsigned i=0;i<nlow;++i){ Tools::convert(i+1,num); dimnames[i] = getLabel() + "." + num; }
  mydata->setNamesAndAtomNumbers( std::vector<AtomNumber>(), dimnames );
}

ProjectNonLandmarkPoints::~ProjectNonLandmarkPoints(){
  delete mydata;
}

void ProjectNonLandmarkPoints::generateProjection( const unsigned& idata, std::vector<double>& point ){
  ConjugateGradient<ProjectNonLandmarkPoints> myminimiser( this );
  unsigned closest=0; double mindist = sqrt( getDissimilarity( idata, mybase->getDataPointIndexInBase(0) ) );
  mybase->setTargetDistance( 0, mindist );
  for(unsigned i=1;i<mybase->getNumberOfDataPoints();++i){
      double dist = sqrt( getDissimilarity( idata, mybase->getDataPointIndexInBase(i) ) );
      mybase->setTargetDistance( i, dist );
      if( dist<mindist ){ mindist=dist; closest=i; }
  }
  // Put the initial guess near to the closest landmark  -- may wish to use grid here again Sandip??
  Random random; random.setSeed(-1234);
  for(unsigned j=0;j<nlow;++j) point[j]=mybase->projections(closest,j) + (random.RandU01() - 0.5)*0.01;
  myminimiser.minimise( cgtol, point, &ProjectNonLandmarkPoints::calculateStress );
}

ReferenceConfiguration* ProjectNonLandmarkPoints::getOutputConfiguration( const unsigned& idata ){
  std::vector<double> pp(nlow); std::vector<double> empty( pp.size() ); generateProjection( idata, pp );
  mydata->setReferenceConfig( std::vector<Vector>(), pp, empty );
  return mydata;
}

///void ProjectNonLandmarkPoints::printAdditionalDataForFrameToPDB( const unsigned& idata, OFile& afile, const std::string& fmt ){
///  std::size_t psign=fmt.find("%"); plumed_assert( psign!=std::string::npos ); std::string num, descr2="%s=%-" + fmt.substr(psign+1);
///  for(unsigned j=0;j<nlow;++j){ Tools::convert(j+1,num); afile.printf(descr2.c_str(), (getLabel()+num).c_str(), projections(idata,j) ); }
///  afile.printf("\n"); 
///} 

unsigned ProjectNonLandmarkPoints::getNumberOfOutputPoints() const {
  return getNumberOfDataPoints();
}

void ProjectNonLandmarkPoints::getOutputForPoint( const unsigned& idata, std::vector<double>& point ){
  if( point.size()!=nlow ) point.resize( nlow );
  generateProjection( idata, point );
}

double ProjectNonLandmarkPoints::getOutputDissimilarity( const unsigned& idata, const unsigned& jdata ){
  std::vector<double> proj1( nlow ), proj2( nlow ); generateProjection( idata, proj1 ); generateProjection( jdata, proj2 );
  double dissim=0; for(unsigned i=0;i<nlow;++i){ double tmp=proj1[i]-proj2[i]; dissim+=tmp*tmp; }
  return dissim;
}

void ProjectNonLandmarkPoints::performAnalysis(){
  // Retrieve the weights from the previous calculation
  std::vector<double> lweights( getNumberOfDataPoints() );
  for(unsigned i=0;i<getNumberOfDataPoints();++i) lweights[i]=getWeight(i);
  setOutputWeights( lweights );

}

double ProjectNonLandmarkPoints::calculateStress( const std::vector<double>& pp, std::vector<double>& der ){
  return mybase->calculateStress( pp, der );
} 

}
}

