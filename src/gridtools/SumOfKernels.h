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
#ifndef __PLUMED_gridtools_SumOfKernels_h
#define __PLUMED_gridtools_SumOfKernels_h

#include "function/FunctionSetup.h"
#include "tools/SwitchingFunction.h"
#include "tools/HistogramBead.h"
#include "tools/Matrix.h"

namespace PLMD {
namespace gridtools {

template <class K>
class RegularKernel;

class DiagonalKernelParams {
public:
  View<const double> at;
  View<const double> sigma;
  View<const double,1> height;
  DiagonalKernelParams( const View<const double>& params, const std::size_t ndim );
  static bool bandwidthIsConstant( std::size_t ndim, const std::vector<Value*>& args );
  static bool bandwidthsAllSame( std::size_t ndim, const std::vector<Value*>& args );
  static bool setKernelAndCheckHeight( View<double>& params, std::size_t ndim, const std::vector<double>& args );
  static std::size_t getNumberOfParameters( const std::size_t& nargs );
  static std::size_t getNumberOfShapeParameters( const std::size_t& nargs );
  static void getSigmaProjections( const View<const double>& shape, std::vector<double>& support );
  static double evaluateR2( const RegularKernel<DiagonalKernelParams>& p, const DiagonalKernelParams& kp, View<const double> x, View<double> paramderivs );
};

class NonDiagonalKernelParams {
public:
  View<const double> at;
  View2D<const double> metric;
  View<const double,1> height;
  NonDiagonalKernelParams( const View<const double>& params, const std::size_t ndim );
  static bool bandwidthIsConstant( std::size_t ndim, const std::vector<Value*>& args );
  static bool bandwidthsAllSame( std::size_t ndim, const std::vector<Value*>& args );
  static bool setKernelAndCheckHeight( View<double>& params, std::size_t ndim, const std::vector<double>& args );
  static std::size_t getNumberOfParameters( const std::size_t& nargs );
  static std::size_t getNumberOfShapeParameters( const std::size_t& nargs );
  static void getSigmaProjections( const View<const double>& shape, std::vector<double>& support );
  static double evaluateR2( const RegularKernel<NonDiagonalKernelParams>& p, const NonDiagonalKernelParams& kp, View<const double> x, View<double> paramderivs );
};

class DiscreteKernel {
public:
  static void registerKeywords( Keywords& keys ) {}
  static void read( DiscreteKernel& p, ActionWithArguments* action, const std::vector<Value*>& args ) {}
  static void setArgumentDomain( const unsigned& i, DiscreteKernel& params, const double& spacing, const bool isp, const std::string& min1, const std::string& max1 ) {}
  static void getSupport( DiscreteKernel& params, const View<const double>& shape, double dp2cutoff, std::vector<double>& support ) {}
  static double calc( const DiscreteKernel& params, const DiagonalKernelParams& kp, View<const double> x, View<double> der, View<double> paramderivs );
  void update() {}
  void toACCDevice() const {}
  void removeFromACCDevice() const  {}
};

class HistogramBeadKernel {
public:
  std::vector<HistogramBead> beads;
  std::vector<double> gridspacing;
  static void registerKeywords( Keywords& keys );
  static void read( HistogramBeadKernel& p, ActionWithArguments* action, const std::vector<Value*>& args );
  static void setArgumentDomain( const unsigned& i, HistogramBeadKernel& params, const double& spacing, const bool isp, const std::string& min1, const std::string& max1 );
  static void getSupport( HistogramBeadKernel& params, const View<const double>& shape, double dp2cutoff, std::vector<double>& support );
  static double calc( const HistogramBeadKernel& params, const DiagonalKernelParams& kp, View<const double> x, View<double> der, View<double> paramderivs );
  void update() {}
  void toACCDevice() const {}
  void removeFromACCDevice() const  {}
};

template <class K>
class RegularKernel {
public:
  bool canusevol;
#ifdef __PLUMED_HAS_OPENACC
  SwitchingFunctionAccelerable switchingFunction;
#else
  SwitchingFunction switchingFunction;
#endif //__PLUMED_HAS_OPENACC
  std::vector<int> periodic_v;
  int* periodic{nullptr};
  std::vector<double> max_minus_min_v;
  double* max_minus_min{nullptr};
  std::vector<double> inv_max_minus_min_v;
  double* inv_max_minus_min{nullptr};
  static void registerKeywords( Keywords& keys );
  static void read( RegularKernel& p, ActionWithArguments* action, const std::vector<Value*>& args );
  static void setArgumentDomain( const unsigned& i, RegularKernel& params, const double& spacing, const bool isp, const std::string& min1, const std::string& max1 );
  static double difference( const RegularKernel& params, unsigned i, const double& val1, const double& val2 );
  static void getSupport( const RegularKernel& params, const View<const double>& shape, double dp2cutoff, std::vector<double>& support );
  static double calc( const RegularKernel& params, const K& kp, View<const double> x, View<double> der, View<double> paramderivs );
  void update() {
    periodic = periodic_v.data();
    max_minus_min = max_minus_min_v.data();
    inv_max_minus_min = inv_max_minus_min_v.data();
  }
#ifdef __PLUMED_HAS_OPENACC
  void toACCDevice() const {
#pragma acc enter data copyin(this[0:1],periodic[0:periodic_v.size()],max_minus_min[0:max_minus_min_v.size()],inv_max_minus_min[0:inv_max_minus_min.size()])
    switchingFunction.toACCDevice();
  }
  void removeFromACCDevice() const  {
    switchingFunction.removeFromACCDevice();
#pragma acc exit data delete(inv_max_minus_min[0:inv_max_minus_min.size()],max_minus_min[0:max_minus_min_v.size()],periodic[0:periodic_v.size()],this[0:1])
  }
#endif //__PLUMED_HAS_OPENACC
};

template <class K>
void RegularKernel<K>::registerKeywords( Keywords& keys ) {
  keys.add("compulsory","KERNEL","GAUSSIAN","the kernel function you are using.");
}

template <class K>
void RegularKernel<K>::read( RegularKernel& p, ActionWithArguments* action, const std::vector<Value*>& args ) {
  std::string kerneltype;
  action->parse("KERNEL",kerneltype);
  std::string errors;
  for(auto & c: kerneltype) {
    c = std::toupper(c);
  }
  p.canusevol = (kerneltype=="GAUSSIAN");
  p.switchingFunction.set( kerneltype + " R_0=1.0 NOSTRETCH", errors );
  if( errors.length()!=0 ) {
    action->error("problem reading switching function description " + errors);
  }
  p.periodic_v.resize( args.size(), 0 );
  p.max_minus_min_v.resize( args.size() );
  p.inv_max_minus_min_v.resize( args.size() );
}

template <class K>
void RegularKernel<K>::setArgumentDomain( const unsigned& i, RegularKernel& params, const double& spacing, const bool isp, const std::string& min1, const std::string& max1 ) {
  if( isp ) {
    double min, max;
    Tools::convert( min1, min );
    Tools::convert( max1, max );
    params.periodic_v[i] = 1;
    params.max_minus_min_v[i]=max-min;
    params.inv_max_minus_min_v[i]=1.0/params.max_minus_min_v[i];
  }
}

template <class K>
double RegularKernel<K>::difference( const RegularKernel<K>& params, unsigned i, const double& val1, const double& val2 ) {
  if( params.periodic[i]==0 ) {
    return val1 - val2;
  }
  return params.max_minus_min[i]*Tools::pbc( params.inv_max_minus_min[i]*( val1 - val2 ) );
}

template <class K>
void RegularKernel<K>::getSupport( const RegularKernel<K>& params, const View<const double>& shape, double dp2cutoff, std::vector<double>& support ) {
  K::getSigmaProjections( shape, support );
  for(unsigned i=0; i<support.size(); ++i) {
    support[i] = sqrt(2.0*dp2cutoff)*support[i];
  }
}

template <class K>
double RegularKernel<K>::calc( const RegularKernel<K>& params, const K& kp, View<const double> x, View<double> der, View<double> paramderivs ) {
  double r2 = K::evaluateR2( params, kp, x, paramderivs );
  double dval, val = kp.height[0]*params.switchingFunction.calculateSqr( r2, dval );
  dval *= kp.height[0];
  for(unsigned i=0; i<der.size(); ++i) {
    der[i] += dval*paramderivs[i];
    paramderivs[i] = -dval*paramderivs[i];
  }
  paramderivs[2*kp.at.size()] = val / kp.height[0];
  return val;
}

class VonMissesKernelParams {
public:
  View<const double,3> at;
  View<const double,1> concentration;
  View<const double,1> height;
  VonMissesKernelParams( const View<const double>& params, const std::size_t ndim );
  static bool bandwidthIsConstant( std::size_t ndim, const std::vector<Value*>& args );
  static bool bandwidthsAllSame( std::size_t ndim, const std::vector<Value*>& args );
  static bool setKernelAndCheckHeight( View<double>& params, std::size_t ndim, const std::vector<double>& argval );
  static std::size_t getNumberOfParameters( const std::size_t& nargs );
  static std::size_t getNumberOfShapeParameters( const std::size_t& nargs );
};

class UniversalVonMisses {
public:
  static void registerKeywords( Keywords& keys ) {}
  static void read( UniversalVonMisses& p, ActionWithArguments* action, const std::vector<Value*>& args ) {}
  static void setArgumentDomain( const unsigned& i, UniversalVonMisses& params, const double& spacing, const bool isp, const std::string& min1, const std::string& max1 ) {}
  static double calc( const UniversalVonMisses& params, const VonMissesKernelParams& kp, View<const double> x, View<double> der, View<double> paramderivs );
  void update() {}
  void toACCDevice() const {}
  void removeFromACCDevice() const  {}
};

template <class K, class P>
class SumOfKernels {
public:
  P params;
  std::vector<double> kernelParams_v;
  double* kernelParams{nullptr};
/// This is used to setup the input gridobject's bounds with the grid data from values
  static void registerKeywords( Keywords& keys );
  static void read( SumOfKernels<K,P>& func, ActionWithArguments* action, const std::vector<Value*>& args, function::FunctionOptions& options );
  static void calc( View<const std::size_t> klist, const SumOfKernels<K,P>& func, View<const double> args, View<double> values, View<double> der, View<double> paramderivs );
  void update() {
    kernelParams=kernelParams_v.data();
  }
  void toACCDevice() const {
#pragma acc enter data copyin(this[0:1],kernelParams[0:kernelParams_v.size()])
    params.toACCDevice();
  }
  void removeFromACCDevice() const  {
    params.removeFromACCDevice();
#pragma acc exit data delete(kernelParams[0:kernelParams_v.size()],this[0:1])
  }
};

template <class K, class P>
void SumOfKernels<K,P>::registerKeywords( Keywords& keys ) {
  P::registerKeywords( keys );
}

template <class K, class P>
void SumOfKernels<K,P>::read( SumOfKernels<K,P>& func, ActionWithArguments* action, const std::vector<Value*>& args, function::FunctionOptions& options ) {
// Read the universal parameters for the kernel
  P::read( func.params, action, args );
}

template <class K, class P>
void SumOfKernels<K,P>::calc( View<const std::size_t> klist, const SumOfKernels<K,P>& func, View<const double> args, View<double> values, View<double> der, View<double> paramderivs ) {
  values[0] = 0;
  for(unsigned i=0; i<der.size(); ++i) {
    der[i] = 0;
  }
  std::size_t ndim = args.size();
  std::size_t nparams = K::getNumberOfParameters( ndim );
  for(unsigned i=0; i<klist.size(); ++i) {
    values[0] += P::calc( func.params, K(View<const double>(func.kernelParams + nparams*klist[i],nparams),ndim), args, der, View<double>( paramderivs.data() + i*nparams, nparams ) );
  }
}

}
}
#endif
