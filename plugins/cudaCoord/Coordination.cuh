/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2024 Daniele Rapetti, The plumed team

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   cudaOnPlumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   cudaOnPlumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

// cfloat for DLB_EPSILON and FLT_EPSILON
#include <cfloat>

#include <iostream>
#include <limits>
#include <numeric>

using std::cerr;

// #define vdbg(...) std::cerr << __LINE__ << ":" << #__VA_ARGS__ << " " <<
// (__VA_ARGS__) << '\n'
#define vdbg(...)

namespace PLMD {
namespace GPU {

// these constant will be used within the kernels
template <typename calculateFloat> struct rationalSwitchParameters {
  calculateFloat dmaxSQ = std::numeric_limits<calculateFloat>::max();
  calculateFloat invr0_2 = 1.0; // r0=1
  calculateFloat stretch = 1.0;
  calculateFloat shift = 0.0;
  int nn = 6;
  int mm = 12;
};

template <typename T> struct invData {
  T val = 1.0;
  T inv = 1.0;
  // this makes the `X = x;` work like "X.val=x;X.inv=1/x;"
  // and the compiler will do some inline magic for you
  invData (T const v) : val{v}, inv{T (1.0) / v} {}
  invData &operator= (T const v) {
    val = v;
    inv = T (1.0) / v;
    return *this;
  }
};
template <typename calculateFloat> struct ortoPBCs {
  invData<calculateFloat> X{1.0};
  invData<calculateFloat> Y{1.0};
  invData<calculateFloat> Z{1.0};
};

template <typename calculateFloat>
__device__ calculateFloat pbcClamp (calculateFloat x) {
  return 0.0;
}

template <> __device__ __forceinline__ double pbcClamp<double> (double x) {
  // convert a double to a signed int in round-to-nearest-even mode.
  // return __double2int_rn (x) - x;
  // return x - floor (x + 0.5);
  // Round argument x to an integer value in single precision floating-point
  // format.
  // Uses round to nearest rounding, with ties rounding to even.
  return x - nearbyint (x);
}

template <> __device__ __forceinline__ float pbcClamp<float> (float x) {
  // convert a double to a signed int in round-to-nearest-even mode.
  // return __float2int_rn (x) - x;
  // return x - floorf (x + 0.5f);
  return x - nearbyintf (x);
}

template <typename calculateFloat>
__device__ __forceinline__ calculateFloat pcuda_fastpow (calculateFloat base,
    int expo) {
  if (expo < 0) {
    expo = -expo;
    base = 1.0 / base;
  }
  calculateFloat result = 1.0;
  while (expo) {
    if (expo & 1) {
      result *= base;
    }
    expo >>= 1;
    base *= base;
  }
  return result;
}

template <typename calculateFloat> __device__ calculateFloat pcuda_eps() {
  return 0;
}

template <> constexpr __device__ float pcuda_eps<float>() {
  return FLT_EPSILON * 10.0f;
}
template <> constexpr __device__ double pcuda_eps<double>() {
  return DBL_EPSILON * 10.0;
}

template <typename calculateFloat>
__device__ __forceinline__ calculateFloat
pcuda_Rational (const calculateFloat rdist,
                const int NN,
                const int MM,
                calculateFloat &dfunc) {
  calculateFloat result;
  if (2 * NN == MM) {
    // if 2*N==M, then (1.0-rdist^N)/(1.0-rdist^M) = 1.0/(1.0+rdist^N)
    calculateFloat rNdist = pcuda_fastpow (rdist, NN - 1);
    result = 1.0 / (1 + rNdist * rdist);
    dfunc = -NN * rNdist * result * result;
  } else {
    if (rdist > (1. - pcuda_eps<calculateFloat>()) &&
        rdist < (1 + pcuda_eps<calculateFloat>())) {

      result = NN / MM;
      dfunc = 0.5 * NN * (NN - MM) / MM;
    } else {
      calculateFloat rNdist = pcuda_fastpow (rdist, NN - 1);
      calculateFloat rMdist = pcuda_fastpow (rdist, MM - 1);
      calculateFloat num = 1. - rNdist * rdist;
      calculateFloat iden = 1.0 / (1.0 - rMdist * rdist);
      result = num * iden;
      dfunc = ((-NN * rNdist * iden) + (result * (iden * MM) * rMdist));
    }
  }
  return result;
}

template <typename calculateFloat>
__global__ void getpcuda_Rational (const calculateFloat *rdists,
                                   const int NN,
                                   const int MM,
                                   calculateFloat *dfunc,
                                   calculateFloat *res) {
  const int i = threadIdx.x + blockIdx.x * blockDim.x;
  if (rdists[i] <= 0.) {
    res[i] = 1.;
    dfunc[i] = 0.0;
  } else {
    res[i] = pcuda_Rational (rdists[i], NN, MM, dfunc[i]);
  }
  // printf("stretch: %i: %f -> %f\n",i,rdists[i],res[i]);
}

template <typename calculateFloat>
__device__ __forceinline__ calculateFloat calculateSqr (
  const calculateFloat distancesq,
  const rationalSwitchParameters<calculateFloat> switchingParameters,
  calculateFloat &dfunc) {
  calculateFloat result = 0.0;
  dfunc = 0.0;
  if (distancesq < switchingParameters.dmaxSQ) {
    const calculateFloat rdist_2 = distancesq * switchingParameters.invr0_2;
    result = pcuda_Rational (
               rdist_2, switchingParameters.nn / 2, switchingParameters.mm / 2, dfunc);
    // chain rule:
    dfunc *= 2 * switchingParameters.invr0_2;
    // cu_stretch:
    result = result * switchingParameters.stretch + switchingParameters.shift;
    dfunc *= switchingParameters.stretch;
  }
  return result;
}

} // namespace GPU
} // namespace PLMD
