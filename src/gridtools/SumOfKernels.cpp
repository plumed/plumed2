/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2023 The plumed team
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
#include "SumOfKernels.h"

namespace PLMD {
namespace gridtools {

double DiscreteKernel::calc( const DiscreteKernel& params, const DiagonalKernelParams& kp, View<const double> x, View<double> der, View<double> paramderivs ) {
  return kp.height;
}

void HistogramBeadKernel::registerKeywords( Keywords& keys ) {
  keys.add("compulsory","KERNEL","GAUSSIAN","the kernel function you are using.");
}

void HistogramBeadKernel::read( HistogramBeadKernel& p, ActionWithArguments* action, const std::vector<Value*>& args ) {
  std::string kerneltype;
  action->parse("KERNEL",kerneltype);
  p.gridspacing.resize( args.size() );
  p.beads.resize( args.size(), HistogramBead(HistogramBead::getKernelType(kerneltype), 0.0, 1.0, 0.5 ) );
}

void HistogramBeadKernel::setArgumentDomain( const unsigned& i, HistogramBeadKernel& params, const double& spacing, const bool isp, const std::string& min1, const std::string& max1 ) {
  params.gridspacing[i] = spacing;
  if( isp ) {
    double lcoord,  ucoord;
    Tools::convert( min1, lcoord );
    Tools::convert( max1, ucoord );
    params.beads[i].isPeriodic( lcoord, ucoord );
  } else {
    params.beads[i].isNotPeriodic();
  }
}

void HistogramBeadKernel::getSupport( HistogramBeadKernel& params, const DiagonalKernelParams& kp, double dp2cutoff, std::vector<double>& support ) {
  for(unsigned i=0; i<support.size(); ++i) {
    params.beads[i].set( 0, params.gridspacing[i], kp.sigma[i] );
    support[i] = params.beads[i].getCutoff();
  }
}

double HistogramBeadKernel::calc( const HistogramBeadKernel& params, const DiagonalKernelParams& kp, View<const double> x, View<double> der, View<double> paramderivs ) {
  double val = kp.height;
  for(unsigned i=0; i<x.size(); ++i) {
    paramderivs[i] = params.beads[i].calculateWithCutoff( kp.at[i], x[i], x[i]+params.gridspacing[i], kp.sigma[i], der[i] );
    val = val*paramderivs[i];
  }
  for(unsigned i=0; i<x.size(); ++i) {
    if( fabs(paramderivs[i])>epsilon ) {
      paramderivs[i] = der[i]*val / paramderivs[i];
      der[i] = 0.0;   // This derivative is set equal to zero because I am not sure what its proper value should be.
    }
  }
  paramderivs[2*kp.at.size()] = val / kp.height;
  return val;
}

bool DiagonalKernelParams::setKernelAndCheckHeight( DiagonalKernelParams& kp, std::size_t ndim, const std::vector<double>& argval ) {
  if( kp.at.size()==ndim && fabs( argval[argval.size()-1])<epsilon ) {
    return false;
  }
  if( kp.at.size()!=ndim ) {
    kp.at.resize( ndim );
    kp.sigma.resize( ndim );
  }
  if( argval.size()==ndim+1 ) {
    for(unsigned i=0; i<ndim; ++i) {
      kp.at[i] = argval[i];
    }
  } else {
    plumed_assert( argval.size()==2*ndim+1 );
    for(unsigned i=0; i<ndim; ++i) {
      kp.at[i] = argval[i];
      kp.sigma[i] = argval[ndim+i];
    }
  }
  kp.height = argval[argval.size()-1];
  return fabs( kp.height )>epsilon;
}

bool DiagonalKernelParams::bandwidthIsConstant( std::size_t ndim, const std::vector<Value*>& args ) {
  for(unsigned i=ndim; i<args.size()-1; ++i) {
    if( !args[i]->isConstant() ) {
      return false;
    }
  }
  return true;
}

bool DiagonalKernelParams::bandwidthsAllSame( std::size_t ndim, const std::vector<Value*>& args ) {
  for(unsigned i=ndim; i<args.size()-1; ++i) {
    if( !args[i]->allElementsEqual() ) {
      return false;
    }
  }
  return true;
}

std::size_t DiagonalKernelParams::getNumberOfParameters( const DiagonalKernelParams& kp ) {
  return kp.at.size() + kp.sigma.size() + 1;
}

void DiagonalKernelParams::getSigmaProjections( const DiagonalKernelParams& kp, std::vector<double>& support ) {
  for(unsigned i=0; i<support.size(); ++i) {
    support[i] = kp.sigma[i];
  }
}

double DiagonalKernelParams::evaluateR2( const RegularKernel<DiagonalKernelParams>& p, const DiagonalKernelParams& kp, View<const double> x, View<double> paramderivs ) {
  double r2 = 0;
  for(unsigned i=0; i<x.size(); ++i) {
    double tmp = RegularKernel<DiagonalKernelParams>::difference( p, i, x[i], kp.at[i] );
    if( fabs(tmp)<epsilon ) {
      paramderivs[i]=0;
    } else {
      paramderivs[i] = tmp/(kp.sigma[i]*kp.sigma[i]);
    }
    r2 += tmp*paramderivs[i];
  }
  return r2;
}

bool NonDiagonalKernelParams::setKernelAndCheckHeight( NonDiagonalKernelParams& kp, std::size_t ndim, const std::vector<double>& argval ) {
  plumed_assert( argval.size()==2*ndim+1 );
  if( kp.at.size()==ndim && fabs( argval[argval.size()-1])<epsilon ) {
    return false;
  }
  if( kp.at.size()!=ndim ) {
    kp.at.resize( ndim );
    kp.sigma.resize( ndim, ndim );
    kp.metric.resize( ndim, ndim );
  }
  for(unsigned i=0; i<ndim; ++i) {
    kp.at[i] = argval[i];
    for(unsigned j=0; j<ndim; ++j) {
      kp.sigma[i][j] = argval[ndim + i*ndim + j];
    }
  }
  Invert( kp.sigma, kp.metric );
  kp.height = argval[argval.size()-1];
  return fabs( kp.height )>epsilon;
}

bool NonDiagonalKernelParams::bandwidthIsConstant( std::size_t ndim, const std::vector<Value*>& args ) {
  return DiagonalKernelParams::bandwidthIsConstant( ndim, args );
}

bool NonDiagonalKernelParams::bandwidthsAllSame( std::size_t ndim, const std::vector<Value*>& args ) {
  return DiagonalKernelParams::bandwidthsAllSame( ndim, args );
}

std::size_t NonDiagonalKernelParams::getNumberOfParameters( const NonDiagonalKernelParams& kp ) {
  return kp.at.size() + kp.at.size()*kp.at.size() + 1;
}

void NonDiagonalKernelParams::getSigmaProjections( const NonDiagonalKernelParams& kp, std::vector<double>& support ) {
  std::vector<double> myautoval(support.size());
  Matrix<double> myautovec(support.size(),support.size());
  diagMat( kp.sigma, myautoval, myautovec);
  double maxautoval=myautoval[0];
  unsigned ind_maxautoval=0;
  for(unsigned i=1; i<support.size(); i++) {
    double neweig=myautoval[i];
    if(neweig>maxautoval) {
      maxautoval=neweig;
      ind_maxautoval=i;
    }
  }
  for(unsigned i=0; i<support.size(); i++) {
    support[i] = fabs(sqrt(maxautoval)*myautovec(i,ind_maxautoval));
  }
}

double NonDiagonalKernelParams::evaluateR2( const RegularKernel<NonDiagonalKernelParams>& p, const NonDiagonalKernelParams& kp, View<const double> x, View<double> paramderivs ) {
  double r2 = 0;
  double dp_j, dp_k;
  for(unsigned j=0; j<x.size(); ++j) {
    dp_j = RegularKernel<NonDiagonalKernelParams>::difference( p, j, x[j], kp.at[j] );
    for(unsigned k=0; k<x.size(); ++k) {
      if( j==k ) {
        dp_k = dp_j;
      } else {
        dp_k = RegularKernel<NonDiagonalKernelParams>::difference( p, k, x[k], kp.at[k] );
      }
      paramderivs[j] += kp.metric[j][k]*dp_k;
      r2 += dp_j*dp_k*kp.metric[j][k];
    }
  }
  return r2;
}

double UniversalVonMisses::calc( const UniversalVonMisses& params, const VonMissesKernelParams& kp, View<const double> x, View<double> der, View<double> paramderivs ) {
  double dot=x[0]*kp.at[0] + x[1]*kp.at[1] + x[2]*kp.at[2];
  double newval = kp.height*exp( kp.concentration*dot );
  der[0] += newval*kp.concentration*kp.at[0];
  der[1] += newval*kp.concentration*kp.at[1];
  der[2] += newval*kp.concentration*kp.at[2];
  paramderivs[0] = newval*kp.concentration*x[0];
  paramderivs[1] = newval*kp.concentration*x[1];
  paramderivs[2] = newval*kp.concentration*x[2];
  if( fabs(kp.height)>epsilon ) {
    paramderivs[4] = newval / kp.height;
  }
  return newval;
}

bool VonMissesKernelParams::bandwidthIsConstant( std::size_t ndim, const std::vector<Value*>& args ) {
  return args[3]->isConstant();
}

bool VonMissesKernelParams::bandwidthsAllSame( std::size_t ndim, const std::vector<Value*>& args ) {
  return args[3]->allElementsEqual();
}

bool VonMissesKernelParams::setKernelAndCheckHeight( VonMissesKernelParams& kp, std::size_t ndim, const std::vector<double>& argval ) {
  plumed_dbg_assert( argval.size()==5 );
  if( kp.at.size()!=3 ) {
    kp.at.resize(3);
  }

  kp.at[0] = argval[0];
  kp.at[1] = argval[1];
  kp.at[2] = argval[2];
  kp.concentration = argval[3];
  kp.height = argval[4];
  return true;
}

std::size_t VonMissesKernelParams::getNumberOfParameters( const VonMissesKernelParams& kp ) {
  return 5;
}

}
}
