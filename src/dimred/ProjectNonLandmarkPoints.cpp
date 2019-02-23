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
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "tools/Random.h"
#include "tools/ConjugateGradient.h"
#include "analysis/AnalysisBase.h"
#include "reference/ReferenceConfiguration.h"
#include "DimensionalityReductionBase.h"
#include "PCA.h"

//+PLUMEDOC DIMRED PROJECT_ALL_ANALYSIS_DATA
/*
Find projections of all non-landmark points using the embedding calculated by a dimensionality reduction optimization calculation.

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace dimred {

class ProjectNonLandmarkPoints : public analysis::AnalysisBase {
private:
/// Tolerance for conjugate gradient algorithm
  double cgtol;
/// Number of diemsions in low dimensional space
  unsigned nlow;
/// The class that calculates the projection of the data that is required
  DimensionalityReductionBase* mybase;
/// Generate a projection of the ith data point - this is called in two routine
  void generateProjection( const unsigned& idat, std::vector<double>& point );
public:
  static void registerKeywords( Keywords& keys );
  ProjectNonLandmarkPoints( const ActionOptions& ao );
/// Get a reference configuration (this returns the projection)
  analysis::DataCollectionObject& getStoredData( const unsigned& idat, const bool& calcdist );
/// Overwrite getArguments so we get arguments from underlying class
  std::vector<Value*> getArgumentList();
/// This does nothing -- projections are calculated when getDataPoint and getReferenceConfiguration are called
  void performAnalysis() {}
/// This just calls calculate stress in the underlying projection object
  double calculateStress( const std::vector<double>& pp, std::vector<double>& der );
/// Overwrite virtual function in ActionWithVessel
  void performTask( const unsigned&, const unsigned&, MultiValue& ) const { plumed_error(); }
};

PLUMED_REGISTER_ACTION(ProjectNonLandmarkPoints,"PROJECT_ALL_ANALYSIS_DATA")

void ProjectNonLandmarkPoints::registerKeywords( Keywords& keys ) {
  analysis::AnalysisBase::registerKeywords( keys );
  keys.add("compulsory","PROJECTION","the projection that you wish to generate out-of-sample projections with");
  keys.add("compulsory","CGTOL","1E-6","the tolerance for the conjugate gradient optimization");
  keys.addOutputComponent("coord","default","the low-dimensional projections of the various input configurations");
}

ProjectNonLandmarkPoints::ProjectNonLandmarkPoints( const ActionOptions& ao ):
  Action(ao),
  analysis::AnalysisBase(ao),
  mybase(NULL)
{
  std::string myproj; parse("PROJECTION",myproj);
  mybase = plumed.getActionSet().selectWithLabel<DimensionalityReductionBase*>( myproj );
  if( !mybase ) error("could not find projection of data named " + myproj );
  // Add the dependency and set the dimensionality
  addDependency( mybase ); nlow = mybase->nlow;
  // Add fake components to the underlying ActionWithValue for the arguments
  std::string num;
  for(unsigned i=0; i<nlow; ++i) {
    Tools::convert(i+1,num); addComponent( "coord-" + num ); componentIsNotPeriodic( "coord-" + num );
  }

  log.printf("  generating out-of-sample projections using projection with label %s \n",myproj.c_str() );
  parse("CGTOL",cgtol);
}

std::vector<Value*> ProjectNonLandmarkPoints::getArgumentList() {
  std::vector<Value*> arglist( analysis::AnalysisBase::getArgumentList() );
  for(unsigned i=0; i<nlow; ++i) arglist.push_back( getPntrToComponent(i) );
  return arglist;
}

void ProjectNonLandmarkPoints::generateProjection( const unsigned& idat, std::vector<double>& point ) {
  PCA* ispca = dynamic_cast<PCA*>( mybase );
  if( ispca ) {
    ispca->getProjection( my_input_data->getStoredData(idat,false), point );
  } else {
    ConjugateGradient<ProjectNonLandmarkPoints> myminimiser( this );
    unsigned closest=0; double mindist = sqrt( getDissimilarity( mybase->getDataPointIndexInBase(0), idat ) );
    mybase->setTargetDistance( 0, mindist );
    for(unsigned i=1; i<mybase->getNumberOfDataPoints(); ++i) {
      double dist = sqrt( getDissimilarity( mybase->getDataPointIndexInBase(i), idat ) );
      mybase->setTargetDistance( i, dist );
      if( dist<mindist ) { mindist=dist; closest=i; }
    }
    // Put the initial guess near to the closest landmark  -- may wish to use grid here again Sandip??
    Random random; random.setSeed(-1234);
    for(unsigned j=0; j<nlow; ++j) point[j]=mybase->projections(closest,j) + (random.RandU01() - 0.5)*0.01;
    myminimiser.minimise( cgtol, point, &ProjectNonLandmarkPoints::calculateStress );
  }
}

analysis::DataCollectionObject& ProjectNonLandmarkPoints::getStoredData( const unsigned& idat, const bool& calcdist ) {
  std::vector<double> pp(nlow); generateProjection( idat, pp ); std::string num;
  analysis::DataCollectionObject& myref=AnalysisBase::getStoredData(idat,calcdist);
  for(unsigned i=0; i<nlow; ++i) { Tools::convert(i+1,num); myref.setArgument( getLabel() + ".coord-" + num, pp[i] ); }
  return myref;
}

double ProjectNonLandmarkPoints::calculateStress( const std::vector<double>& pp, std::vector<double>& der ) {
  return mybase->calculateStress( pp, der );
}

}
}

