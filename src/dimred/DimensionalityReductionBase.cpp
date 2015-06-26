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
#include "tools/ConjugateGradient.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "reference/PointWiseMapping.h"

namespace PLMD {
namespace analysis {

void DimensionalityReductionBase::registerKeywords( Keywords& keys ){
  AnalysisWithLandmarks::registerKeywords( keys );
  keys.add("compulsory","NLOW_DIM","number of low-dimensional coordinates required");
  keys.add("compulsory","OUTPUT_FILE","file on which to output the final embedding coordinates");
  keys.add("compulsory","EMBEDDING_OFILE","dont output","file on which to output the embedding in plumed input format");
  keys.add("compulsory","OSAMPLE_CGTOL","1E-6","the tolerance on the out of sample conjugate gradient minimisation");
  keys.add("compulsory","OSAMPLE_FILE","dont output","the file on which the projections of all data points are output");
}

DimensionalityReductionBase::DimensionalityReductionBase( const ActionOptions& ao ):
Action(ao),
AnalysisWithLandmarks(ao)
{
  myembedding = new PointWiseMapping( getMetricName(), false );
  setDataToAnalyze( dynamic_cast<MultiReferenceBase*>(myembedding) );

  parse("NLOW_DIM",nlow); parse("OSAMPLE_CGTOL",cgtol);
  if( nlow<1 ) error("dimensionality of low dimensional space must be at least one");

  // Setup the property names
  std::vector<std::string> propnames( nlow ); std::string num;
  for(unsigned i=0;i<propnames.size();++i){
     Tools::convert(i+1,num); std::string lab=getLabel();
     if(lab.find("@")!=std::string::npos) propnames[i]=getName() + "." + num;
     else propnames[i]=getLabel() + "." + num;
  }
  myembedding->setPropertyNames( propnames, false );

  parseOutputFile("EMBEDDING_OFILE",efilename);
  parseOutputFile("OUTPUT_FILE",ofilename);
}

DimensionalityReductionBase::~DimensionalityReductionBase(){
  delete myembedding;
}

void DimensionalityReductionBase::analyzeLandmarks(){
  // These hold data on the distances from high dimensional points
  fframes.resize( getNumberOfLandmarks() );
  targetDisimilarities.resize( getNumberOfLandmarks(), getNumberOfLandmarks() );
  // This calculates all the distanaces between the high dimensional points
  calculateAllDistances( myembedding, targetDisimilarities );
  // This generates the projections of the points
  generateProjections( myembedding );

  // Output the embedding as long lists of data
  OFile gfile; gfile.link(*this); 
  gfile.setBackupString("analysis");
  gfile.fmtField(getOutputFormat()+" ");
  gfile.open( ofilename.c_str() );
  
  // Print embedding coordinates
  for(unsigned i=0;i<getNumberOfLandmarks();++i){
      for(unsigned j=0;j<nlow;++j){
          std::string num; Tools::convert(j+1,num);
          gfile.printField( getLabel() + "." + num , myembedding->getProjectionCoordinate(i,j) );
      }
      gfile.printField();
  }  
  gfile.close();

  // Output the embedding in plumed format
  if( efilename!="dont output"){
     // std::string ifname=saveResultsFromPreviousAnalyses( efilename );
     OFile afile; afile.link(*this); afile.setBackupString("analysis");
     afile.open( efilename.c_str() );
     if( plumed.getAtoms().usingNaturalUnits() ) myembedding->print( 1.0, getAlgorithmName(), getTime(), afile, getOutputFormat() );
     else myembedding->print( plumed.getAtoms().getUnits().getLength()/0.1, getAlgorithmName(), getTime(), afile, getOutputFormat() );
     afile.close();
  }
}

void DimensionalityReductionBase::findClosestPoint( const unsigned& ii, std::vector<double>& pp ){
  plumed_dbg_assert( ii<getNumberOfDataPoints() &&  pp.size()==myembedding->getNumberOfProperties() );
  unsigned pnum; double mindist, df;

  if( ii==getLandmarkIndex(0) ){
      pnum=1; fframes[0]=0.0; mindist=fframes[1]=transformHD( getDistanceBetweenFrames( getLandmarkIndex(1), ii, false ), df );
      for(unsigned i=2;i<myembedding->getNumberOfReferenceFrames();++i){
          fframes[i]=transformHD( getDistanceBetweenFrames( getLandmarkIndex(i), ii, false ), df );
          if( fframes[i]<mindist ){ mindist=fframes[i]; pnum=i; }
      } 
  } else {
      pnum=0; mindist=fframes[0]=transformHD( getDistanceBetweenFrames( getLandmarkIndex(0), ii, false ), df );
      for(unsigned i=1;i<myembedding->getNumberOfReferenceFrames();++i){
          fframes[i]=transformHD( getDistanceBetweenFrames( getLandmarkIndex(i), ii, false ), df );
          if( ii==getLandmarkIndex(i) ) continue;
          if( fframes[i]<mindist ){ mindist=fframes[i]; pnum=i; }
      }
  } 

  for(unsigned i=0;i<myembedding->getNumberOfProperties();++i) pp[i]=myembedding->getProjectionCoordinate( pnum, i );
}

double DimensionalityReductionBase::calculateStress( const std::vector<double>& pp, std::vector<double>& der ){
  plumed_dbg_assert( pp.size()==myembedding->getNumberOfProperties() && der.size()==myembedding->getNumberOfProperties() );

  std::vector<double> tmpder( myembedding->getNumberOfProperties() );
  double df, chi2=0.0; der.assign(der.size(),0.0);
  for(unsigned i=0;i<myembedding->getNumberOfReferenceFrames();++i){
     // Ensure that the identical point is skipped if we are doing point wise global optimisation
     if( fframes[i]<epsilon ) continue ;

     double dist=0.0;
     for(unsigned j=0;j<myembedding->getNumberOfProperties();++j){
         tmpder[j] = pp[j] - myembedding->getProjectionCoordinate( i, j );
         dist += tmpder[j]*tmpder[j];
     }
     dist=transformLD( sqrt(dist), df );
     double tmp = fframes[i] - dist;
     // Accumulate the stress 
     chi2 += myembedding->getWeight(i)*tmp*tmp;
     // Accumulate the derivatives
     for(unsigned j=0;j<myembedding->getNumberOfProperties();++j) der[j] += 2*myembedding->getWeight(i)*tmp*df*tmpder[j];
  }
  return chi2;
}

void DimensionalityReductionBase::getPropertyNames( std::vector<std::string>& dimnames ){
  plumed_dbg_assert( dimnames.size()==nlow );
  for(unsigned i=0;i<nlow;++i) dimnames[i]=myembedding->getPropertyName(i);
}

void DimensionalityReductionBase::getProjectedPoint( const unsigned& idata, std::vector<double>& pp ){
  plumed_dbg_assert( pp.size()==nlow );
  findClosestPoint( idata, pp );
  ConjugateGradient<DimensionalityReductionBase> myminimiser( this );
  myminimiser.minimise( cgtol, pp, &DimensionalityReductionBase::calculateStress );
}

}
}
