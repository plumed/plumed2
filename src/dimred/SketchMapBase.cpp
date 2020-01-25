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
#include "SketchMapBase.h"

namespace PLMD {
namespace dimred {

void SketchMapBase::registerKeywords( Keywords& keys ) {
  DimensionalityReductionBase::registerKeywords( keys );
  keys.remove("NLOW_DIM");
  keys.add("compulsory","HIGH_DIM_FUNCTION","as in input action","the parameters of the switching function in the high dimensional space");
  keys.add("compulsory","LOW_DIM_FUNCTION","as in input action","the parameters of the switching function in the low dimensional space");
  keys.add("compulsory","MIXPARAM","0.0","the amount of the pure distances to mix into the stress function");
}

SketchMapBase::SketchMapBase( const ActionOptions& ao ):
  Action(ao),
  DimensionalityReductionBase(ao),
  smapbase(NULL),
  normw(0.0)
{
  // Check if we have data from a input sketch-map object - we can reuse switching functions wahoo!!
  if( dimredbase ) smapbase = dynamic_cast<SketchMapBase*>( dimredbase );

  // Read in the switching functions
  std::string linput,hinput, errors;
  parse("HIGH_DIM_FUNCTION",hinput);
  if( hinput=="as in input action" ) {
    if( !smapbase ) error("high dimensional switching funciton has not been set - use HIGH_DIM_FUNCTION");
    reuse_hd=true;
    log.printf("  reusing high dimensional filter function defined in previous sketch-map action\n");
  } else {
    reuse_hd=false;
    highdf.set(hinput,errors);
    if(errors.length()>0) error(errors);
    log.printf("  filter function for dissimilarities in high dimensional space has cutoff %s \n",highdf.description().c_str() );
  }

  parse("LOW_DIM_FUNCTION",linput);
  if( linput=="as in input action" ) {
    if( !smapbase ) error("low dimensional switching funciton has not been set - use LOW_DIM_FUNCTION");
    reuse_ld=true;
    log.printf("  reusing low dimensional filter function defined in previous sketch-map action\n");
  } else {
    reuse_ld=false;
    lowdf.set(linput,errors);
    if(errors.length()>0) error(errors);
    log.printf("  filter function for distances in low dimensionality space has cutoff %s \n",lowdf.description().c_str() );
  }

  // Read the mixing parameter
  parse("MIXPARAM",mixparam);
  if( mixparam<0 || mixparam>1 ) error("mixing parameter must be between 0 and 1");
  log.printf("  mixing %f of pure distances with %f of filtered distances \n",mixparam,1.-mixparam);
}

void SketchMapBase::calculateProjections( const Matrix<double>& targets, Matrix<double>& projections ) {
  if( dtargets.size()!=targets.nrows() ) {
    // These hold data so that we can do stress calculations
    dtargets.resize( targets.nrows() ); ftargets.resize( targets.nrows() ); pweights.resize( targets.nrows() );
    // Matrices for storing input data
    transformed.resize( targets.nrows(), targets.ncols() );
    distances.resize( targets.nrows(), targets.ncols() );
  }

  // Stores the weights in an array for faster access, as well as the normalization
  normw=0;
  for(unsigned i=0; i<targets.nrows() ; ++i) { pweights[i] = getWeight(i); normw+=pweights[i]; }
  normw*=normw;

  // Transform the high dimensional distances
  double df; distances=0.; transformed=0.;
  for(unsigned i=1; i<distances.ncols(); ++i) {
    for(unsigned j=0; j<i; ++j) {
      distances(i,j)=distances(j,i)=sqrt( targets(i,j) );
      transformed(i,j)=transformed(j,i)=transformHighDimensionalDistance( distances(i,j), df );
    }
  }
  // And minimse
  minimise( projections );
}

double SketchMapBase::calculateStress( const std::vector<double>& p, std::vector<double>& d ) {
  // Zero derivative and stress accumulators
  for(unsigned i=0; i<p.size(); ++i) d[i]=0.0;
  double stress=0; std::vector<double> dtmp( p.size() );
  // Now accumulate total stress on system
  for(unsigned i=0; i<ftargets.size(); ++i) {
    if( dtargets[i]<epsilon ) continue ;

    // Calculate distance in low dimensional space
    double dd=0;
    for(unsigned j=0; j<p.size(); ++j) { dtmp[j]=p[j]-projections(i,j); dd+=dtmp[j]*dtmp[j]; }
    dd = sqrt(dd);

    // Now do transformations and calculate differences
    double df, fd = transformLowDimensionalDistance( dd, df );
    double ddiff = dd - dtargets[i];
    double fdiff = fd - ftargets[i];

    // Calculate derivatives
    double pref = 2.*getWeight(i) / dd ;
    for(unsigned j=0; j<p.size(); ++j) d[j] += pref*( (1-mixparam)*fdiff*df + mixparam*ddiff )*dtmp[j];

    // Accumulate the total stress
    stress += getWeight(i)*( (1-mixparam)*fdiff*fdiff + mixparam*ddiff*ddiff );
  }
  return stress;
}

double SketchMapBase::calculateFullStress( const std::vector<double>& p, std::vector<double>& d ) {
  // Zero derivative and stress accumulators
  for(unsigned i=0; i<p.size(); ++i) d[i]=0.0;
  double stress=0; std::vector<double> dtmp( nlow );

  for(unsigned i=1; i<distances.nrows(); ++i) {
    double iweight = pweights[i];
    for(unsigned j=0; j<i; ++j) {
      double jweight =  pweights[j];
      // Calculate distance in low dimensional space
      double dd=0;
      for(unsigned k=0; k<nlow; ++k) { dtmp[k]=p[nlow*i+k] - p[nlow*j+k]; dd+=dtmp[k]*dtmp[k]; }
      dd = sqrt(dd);

      // Now do transformations and calculate differences
      double df, fd = transformLowDimensionalDistance( dd, df );
      double ddiff = dd - distances(i,j);
      double fdiff = fd - transformed(i,j);;

      // Calculate derivatives
      double pref = 2.*iweight*jweight*( (1-mixparam)*fdiff*df + mixparam*ddiff ) / dd;
      for(unsigned k=0; k<nlow; ++k) {
        double dterm=pref*dtmp[k]; d[nlow*i+k]+=dterm; d[nlow*j+k]-=dterm;
      }

      // Accumulate the total stress
      stress += iweight*jweight*( (1-mixparam)*fdiff*fdiff + mixparam*ddiff*ddiff );
    }
  }
  stress /= normw; for (unsigned k=0; k < d.size(); ++k) d[k] /= normw;
  return stress;
}

}
}


