/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013 The plumed team
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
#include "core/ActionWithValue.h"
#include "core/ActionPilot.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/ActionRegister.h"
#include "ActionWithVessel.h"
#include "ActionWithInputVessel.h"
#include "FunctionOnGrid.h"
#include "GridVesselBase.h"
#include "InterpolationBase.h"
#include "NearestNeighborInterpolation.h"

namespace PLMD {
namespace gridtools {

class GridStats : 
  public ActionWithValue,
  public vesselbase::ActionWithInputVessel
{
private:
  bool noskew, serial;
  unsigned dimension, npoints;
  double i2sigma2, cellvolume;
  std::vector<unsigned> indices, ngrid;
  std::vector<double> min, pos, delx;
  GridVesselBase* myf;
  InterpolationBase* myfield;
  void getCoordinates( const unsigned& );
public:
  static void registerKeywords( Keywords& keys );
  GridStats(const ActionOptions& ao);
  ~GridStats();
  void calculate();
  void apply(){}
  unsigned getNumberOfDerivatives();
};

PLUMED_REGISTER_ACTION(GridStats,"GRIDSTATS")

void GridStats::registerKeywords( Keywords& keys ){
  Action::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  vesselbase::ActionWithInputVessel::registerKeywords( keys );
  keys.remove("DATA"); keys.use("FUNC");
  keys.add("compulsory","SIGMA","The sigma parameter");
  keys.add("compulsory","NGRIDPOINTS","the number of gridpoints to use for the integration");
  keys.add("compulsory","INTERPOLATION","cubic","what algorithm should be used for interpolation");
  keys.addFlag("NOSKEW",false,"do not bother calcualting the skewness and kurtosis of the field");
  keys.addFlag("SERIAL",false,"do the calculation in serial.  Do not parallelize");
  keys.addOutputComponent("_mean","default","the global minimum calculated as described above");
  keys.addOutputComponent("_var","default","the eigenvalues of covariance matrix for the distribution described above");
  keys.addOutputComponent("varthet","default","the angle between the eigenvector of the covariance with the largest eigenvalue and the x axis (2d distribution only)");
  keys.addOutputComponent("_skew","NOSKEW","the skewness projected along the eigenvectors of the covariance matrix");
  keys.addOutputComponent("_kurt","NOSKEW","the kurtosis projected along the eigenvectors of the covariance matrix");
}

GridStats::GridStats(const ActionOptions& ao):
Action(ao),
ActionWithValue(ao),
ActionWithInputVessel(ao),
myf(NULL),
myfield(NULL)
{
  readArgument( "func" );
  myf = dynamic_cast<GridVesselBase*>( getPntrToArgument() );   

  // Create interpolators for fields
  std::string interpols; parse("INTERPOLATION",interpols);
  FunctionOnGrid* testf=dynamic_cast<FunctionOnGrid*>( myf );
  if( testf && interpols!="nearest" ) error("cannot interpolate function if there are no derivatives");

  parseVector("NGRIDPOINTS",ngrid); dimension=myf->getDimension();
  if( ngrid.size()!=dimension ) error("mismatched dimensionality between field and grid points");

  if( interpols=="cubic" ){
     log.printf("  using cubically interpolated grid \n");
     myfield = InterpolationBase::createCubicInterpolator( myf, 0 );
  } else if ( interpols=="nearest" ){
     log.printf("  no interpolation of grid \n");
     std::vector<unsigned> nbin( myf->getNbin() );
     for(unsigned i=0;i<dimension;++i){
         if( nbin[i]!=ngrid[i] ){
             ngrid[i]=nbin[i];
             warning("mismatch between number of calculated points and number of integration points.  Using number of calculated points");
         }
     }
     myfield = new NearestNeighborInterpolation( myf, 0 );
  } else {
     error(interpols + " is not a valid interpolation algorithm");
  }

  // Create the input for the bias grid 
  npoints=1; cellvolume=1;
  double max; pos.resize( dimension ); indices.resize( dimension );
  min.resize( dimension ); delx.resize( dimension );
  std::vector<std::string> gmin( myf->getMin() ), gmax( myf->getMax() );
  for(unsigned i=0;i<dimension;++i){
     Tools::convert( gmin[i], min[i] );
     Tools::convert( gmax[i], max );
     delx[i] = ( max - min[i] ) / static_cast<double>( ngrid[i] );
     ngrid[i] += 1; npoints *= ngrid[i]; cellvolume*=delx[i];
  }

  // Read in sigma parameter
  double sigma; parse("SIGMA",sigma); 
  i2sigma2= 1. / (2.*sigma*sigma); 

  parseFlag("NOSKEW",noskew);
  parseFlag("SERIAL",serial);

  for(unsigned i=0;i<dimension;++i){
     addComponent(myf->getQuantityDescription(i)+"_mean"); 
     componentIsNotPeriodic(myf->getQuantityDescription(i)+"_mean");
  }
  for(unsigned i=0;i<dimension;++i){
     addComponent(myf->getQuantityDescription(i)+"_var"); 
     componentIsNotPeriodic(myf->getQuantityDescription(i)+"_var");
  }
  // variance/skewness/kurtosis are always given in directions parallel to the 
  // eigenvectors of the covariance matrix.  If dimension=2 then it is 
  // straightforward to calculate the angle of the eigenvectors to the x and
  // y axis so this is displayed.  For higher dimensionalities this is harder 
  // to calculate so we don't bother.
  if( dimension==2 ){
     addComponent("var-thet"); componentIsNotPeriodic("var-thet");
  }
  if(!noskew){
     for(unsigned i=0;i<dimension;++i){
       addComponent(myf->getQuantityDescription(i)+"_skew"); 
       componentIsNotPeriodic(myf->getQuantityDescription(i)+"_skew");
     }
     for(unsigned i=0;i<dimension;++i){
       addComponent(myf->getQuantityDescription(i)+"_kurt"); 
       componentIsNotPeriodic(myf->getQuantityDescription(i)+"_kurt");
     }
  }
}

GridStats::~GridStats(){
  delete myfield; 
}

unsigned GridStats::getNumberOfDerivatives(){
  return 0;
}

void GridStats::getCoordinates( const unsigned& index ){
 unsigned kk=index;
 indices[0]=index%ngrid[0];
 for(unsigned i=1;i<dimension-1;++i){
   kk=(kk-indices[i-1])/ngrid[i-1];
   indices[i]=kk%ngrid[i];
 }
 if(dimension>=2){  // I think this is wrong
    indices[dimension-1]=(kk-indices[dimension-2])/ngrid[dimension-2];
 }
 for(unsigned i=0;i<dimension;++i) pos[i] = min[i] + indices[i]*delx[i];
}

void GridStats::calculate(){
  unsigned rank, stride;

  if( serial ){
     rank=0; stride=1;
  } else {
     rank=comm.Get_rank();
     stride=comm.Get_size();
  }

  // Setup the interpolator for the fields
  myfield->set_table();

  // Calculate the position of the mean and the covariance
  double norm=0.0; Matrix<double> covar( dimension, dimension ); 
  covar=0.0; std::vector<double> mean( dimension, 0.0 );
  for(unsigned i=rank;i<npoints;i+=stride){
      getCoordinates( i );
      double myspot = exp( -i2sigma2*myfield->getFunctionValue( pos ) );
      norm += myspot;  
      for(unsigned j=0;j<dimension;++j) mean[j] += myspot*pos[j]; 
      for(unsigned j=0;j<dimension;++j){
          for(unsigned k=0;k<dimension;++k) covar(j,k) += myspot*pos[j]*pos[k];
      }
  }

  // Finish calculation of mean and covariance
  norm *= cellvolume; comm.Sum( norm ); 
  for(unsigned i=0;i<dimension;++i){ 
      mean[i] *= cellvolume / norm;
      for(unsigned j=0;j<dimension;++j) covar(i,j) *= cellvolume / norm;
  }
  comm.Sum( mean ); comm.Sum( covar );
  for(unsigned i=0;i<dimension;++i){
      for(unsigned j=0;j<dimension;++j) covar(i,j) -= mean[i]*mean[j];
  }

  // Now calculate the eigenvalues and eigenvectors of the covariance matrix
  std::vector<double> eigvals( dimension ); 
  Matrix<double> eigvecs( dimension, dimension );
  diagMat( covar, eigvals, eigvecs );

  // Now set the mean and covariance values
  for(unsigned i=0;i<dimension;++i){
     getPntrToComponent(myf->getQuantityDescription(i)+"_mean")->set( mean[i] );  
     getPntrToComponent(myf->getQuantityDescription(i)+"_var")->set( eigvals[i] );
  }
  if( dimension==2 ) getPntrToComponent("var-thet")->set( atan2( eigvecs(0,1), eigvecs(0,0) ) );

  if(noskew) return;

  // Calculate the skewness and kurtosis
  std::vector<double> proj( dimension, 0.0 );
  std::vector<double> skew( dimension, 0.0 ), kurt( dimension, 0.0 );
  for(unsigned i=rank;i<npoints;i+=stride){
      getCoordinates( i );
      double myspot = exp( -i2sigma2*myfield->getFunctionValue( pos ) );
      for(unsigned j=0;j<dimension;++j) pos[j]-=mean[j];
      for(unsigned j=0;j<dimension;++j){
         proj[j]=0.;
         for(unsigned k=0;k<dimension;++k) proj[j] += eigvecs(j,k)*pos[k];
         double temp2=proj[j]*proj[j];
         skew[j] += temp2*proj[j]*myspot; kurt[j] += temp2*temp2*myspot;
      }
  }

  // Finish calculation of skewness and kurtosis 
  for(unsigned i=0;i<dimension;++i){
      skew[i] *= cellvolume / norm;
      kurt[i] *= cellvolume / norm;
  }
  comm.Sum( skew ); comm.Sum( kurt );

  for(unsigned i=0;i<dimension;++i){
     getPntrToComponent(myf->getQuantityDescription(i)+"_skew")->set( skew[i]/sqrt(eigvals[i]*eigvals[i]*eigvals[i]) );
     getPntrToComponent(myf->getQuantityDescription(i)+"_kurt")->set( kurt[i]/(eigvals[i]*eigvals[i]) - 3.0 );
  }
}

}
}
