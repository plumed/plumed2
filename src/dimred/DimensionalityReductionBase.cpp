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
#include "DimensionalityReductionBase.h"
#include "reference/ReferenceConfiguration.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"

namespace PLMD {
namespace dimred {

void DimensionalityReductionBase::registerKeywords( Keywords& keys ) {
  analysis::AnalysisBase::registerKeywords( keys );
  keys.add("compulsory","NLOW_DIM","number of low-dimensional coordinates required");
  keys.addOutputComponent("coord","default","the low-dimensional projections of the various input configurations");
}

DimensionalityReductionBase::DimensionalityReductionBase( const ActionOptions& ao ):
  Action(ao),
  analysis::AnalysisBase(ao),
  dimredbase(NULL)
{
  // Check that some dissimilarity information is available
  if( my_input_data ) {
    if( getName()!="PCA" && !dissimilaritiesWereSet() ) error("dissimilarities have not been calcualted in input actions");
    // Now we check if the input was a dimensionality reduction object
    dimredbase = dynamic_cast<DimensionalityReductionBase*>( my_input_data );
  }

  // Retrieve the dimension in the low dimensionality space
  nlow=0;
  if( dimredbase ) {
    nlow=dimredbase->nlow;
    log.printf("  projecting in %u dimensional space \n",nlow);
  } else if( keywords.exists("NLOW_DIM") ) {
    parse("NLOW_DIM",nlow);
    if( nlow<1 ) error("dimensionality of low dimensional space must be at least one");
    log.printf("  projecting in %u dimensional space \n",nlow);
  }
  // Now add fake components to the underlying ActionWithValue for the arguments
  std::string num;
  for(unsigned i=0; i<nlow; ++i) {
    Tools::convert(i+1,num); addComponent( "coord-" + num ); componentIsNotPeriodic( "coord-" + num );
  }
}

std::vector<Value*> DimensionalityReductionBase::getArgumentList() {
  std::vector<Value*> arglist( analysis::AnalysisBase::getArgumentList() );
  for(unsigned i=0; i<nlow; ++i) arglist.push_back( getPntrToComponent(i) );
  return arglist;
}

void DimensionalityReductionBase::getProjection( const unsigned& idata, std::vector<double>& point, double& weight ) {
  if( point.size()!=nlow ) point.resize( nlow );
  weight = getWeight(idata); for(unsigned i=0; i<nlow; ++i) point[i]=projections(idata,i);
}

void DimensionalityReductionBase::performAnalysis() {
  log.printf("Generating projections required by action %s \n",getLabel().c_str() );
  // Resize the tempory array (this is used for out of sample)
  dtargets.resize( getNumberOfDataPoints() );
  // Resize the projections array
  projections.resize( getNumberOfDataPoints(), nlow );
  // Retreive the projections from the previous calculation
  if( dimredbase ) {
    std::vector<double> newp( nlow ); double w;
    for(unsigned i=0; i<getNumberOfDataPoints(); ++i) {
      dimredbase->getProjection( i, newp, w ); plumed_dbg_assert( newp.size()==nlow );
      for(unsigned j=0; j<nlow; ++j) projections(i,j)=newp[j];
    }
  }
  // Calculate matrix of dissimilarities
  Matrix<double> targets( getNumberOfDataPoints(), getNumberOfDataPoints() ); targets=0;
  for(unsigned i=1; i<getNumberOfDataPoints(); ++i) {
    for(unsigned j=0; j<i; ++j) targets(i,j)=targets(j,i)=getDissimilarity( i, j );
  }
  // This calculates the projections of the points
  calculateProjections( targets, projections );
  // Now set the projection values in the underlying object
  if( my_input_data ) {
    for(unsigned idat=0; idat<getNumberOfDataPoints(); ++idat) {
      analysis::DataCollectionObject& myref=AnalysisBase::getStoredData(idat,false); std::string num;
      for(unsigned i=0; i<nlow; ++i) { Tools::convert(i+1,num); myref.setArgument( getLabel() + ".coord-" + num, projections(idat,i) ); }
    }
  }
  log.printf("Generated projections required by action %s \n",getLabel().c_str() );
}

double DimensionalityReductionBase::calculateStress( const std::vector<double>& p, std::vector<double>& d ) {

  // Zero derivative and stress accumulators
  for(unsigned i=0; i<p.size(); ++i) d[i]=0.0;
  double stress=0; std::vector<double> dtmp( p.size() );

  // Now accumulate total stress on system
  for(unsigned i=0; i<dtargets.size(); ++i) {
    if( dtargets[i]<epsilon ) continue ;

    // Calculate distance in low dimensional space
    double dd=0;
    for(unsigned j=0; j<p.size(); ++j) { dtmp[j]=p[j]-projections(i,j); dd+=dtmp[j]*dtmp[j]; }
    dd = sqrt(dd);

    // Now do transformations and calculate differences
    double ddiff = dd - dtargets[i];

    // Calculate derivatives
    double pref = 2.*getWeight(i) / dd;
    for(unsigned j=0; j<p.size(); ++j) d[j] += pref*ddiff*dtmp[j];

    // Accumulate the total stress
    stress += getWeight(i)*ddiff*ddiff;
  }
  return stress;
}

}
}
