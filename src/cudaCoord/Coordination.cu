/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2023 The plumed team
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
#include "plumed/colvar/CoordinationBase.h"
#include "plumed/core/ActionRegister.h"
#include "plumed/tools/NeighborList.h"
#include "plumed/tools/SwitchingFunction.h"

#include "cudaHelpers.cuh"
#include "ndReduction.h"

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <iostream>
#include <limits>
#include <numeric>

using std::cerr;

#define vdbg(...)                                                              \
  std::cerr << __LINE__ << ":" << #__VA_ARGS__ << " " << (__VA_ARGS__) << '\n'
// #define vdbg(...)

namespace PLMD {
namespace colvar {
//+PLUMEDOC COLVAR CUDACOORDINATION
/*
Calculate coordination numbers. Like coordination, but on nvdia gpu and with no
swithcing function

This keyword can be used to calculate the number of contacts between two groups
of atoms and is defined as \f[ \sum_{i\in A} \sum_{i\in B} s_{ij} \f] where
\f$s_{ij}\f$ is 1 if the contact between atoms \f$i\f$ and \f$j\f$ is formed,
zero otherwise.
In actuality, \f$s_{ij}\f$ is replaced with a switching function so as to ensure
that the calculated CV has continuous derivatives. The default switching
function is: \f[ s_{ij} = \frac{ 1 - \left(\frac{{\bf r}_{ij}-d_0}{r_0}\right)^n
} { 1 - \left(\frac{{\bf r}_{ij}-d_0}{r_0}\right)^m } \f] but it can be changed
using the optional SWITCH option.

To make your calculation faster you can use a neighbor list, which makes it that
only a relevant subset of the pairwise distance are calculated at every step.

If GROUPB is empty, it will sum the \f$\frac{N(N-1)}{2}\f$ pairs in GROUPA. This
avoids computing twice permuted indexes (e.g. pair (i,j) and (j,i)) thus running
at twice the speed.

Notice that if there are common atoms between GROUPA and GROUPB the switching
function should be equal to one. These "self contacts" are discarded by plumed
(since version 2.1), so that they actually count as "zero".


\par Examples

The following example instructs plumed to calculate the total coordination
number of the atoms in group 1-10 with the atoms in group 20-100.  For atoms
1-10 coordination numbers are calculated that count the number of atoms from the
second group that are within 0.3 nm of the central atom.  A neighbor list is
used to make this calculation faster, this neighbor list is updated every 100
steps. \plumedfile COORDINATION GROUPA=1-10 GROUPB=20-100 R_0=0.3 NLIST
NL_CUTOFF=0.5 NL_STRIDE=100 \endplumedfile

The following is a dummy example which should compute the value 0 because the
self interaction of atom 1 is skipped. Notice that in plumed 2.0 "self
interactions" were not skipped, and the same calculation should return 1.
\plumedfile
c: COORDINATION GROUPA=1 GROUPB=1 R_0=0.3
PRINT ARG=c STRIDE=10
\endplumedfile

Here's an example that shows what happens when providing COORDINATION with
a single group:
\plumedfile
# define some huge group:
group: GROUP ATOMS=1-1000
# Here's coordination of a group against itself:
c1: COORDINATION GROUPA=group GROUPB=group R_0=0.3
# Here's coordination within a single group:
x: COORDINATION GROUPA=group R_0=0.3
# This is just multiplying times 2 the variable x:
c2: COMBINE ARG=x COEFFICIENTS=2 PERIODIC=NO

# the two variables c1 and c2 should be identical, but the calculation of c2 is
twice faster # since it runs on half of the pairs. PRINT ARG=c1,c2 STRIDE=10
\endplumedfile



*/
//+ENDPLUMEDOC

// these constant will be used within the kernels
template <typename calculateFloat> struct rationalSwitchParameters {
  calculateFloat dmaxSQ = std::numeric_limits<double>::max();
  calculateFloat invr0_2 = 1.0; // r0=1
  calculateFloat stretch = 1.0;
  calculateFloat shift = 0.0;
  int nn = 6;
  int mm = 12;
};

template <typename calculateFloat> struct ortoPBCs {
  calculateFloat invX = 1.0;
  calculateFloat invY = 1.0;
  calculateFloat invZ = 1.0;
  calculateFloat X = 1.0;
  calculateFloat Y = 1.0;
  calculateFloat Z = 1.0;
};

template <typename calculateFloat>
__device__ calculateFloat pbcClamp(calculateFloat x) {
  return 0.0;
}

template <> __device__ double pbcClamp<double>(double x) {
  // convert a double to a signed int in round-to-nearest-even mode.
  return __double2int_rn(x) - x;
  // return x - floor(x+0.5);
  // Round argument x to an integer value in single precision floating-point
  // format.
  // Uses round to nearest rounding, with ties rounding to even.
  // return nearbyint(x) - x;
}

template <> __device__ float pbcClamp<float>(float x) {
  // convert a double to a signed int in round-to-nearest-even mode.
  return __float2int_rn(x) - x;
  // return x - floorf(x+0.5f);
  // return nearbyintf(x) - x;
}

template <typename T> constexpr cudaDataType cudaPrecision() { return 0; }
template <> constexpr cudaDataType cudaPrecision<double>() {
  return CUDA_R_64F;
}
template <> constexpr cudaDataType cudaPrecision<float>() { return CUDA_R_32F; }

// does not inherit from coordination base because nl is private
template <typename calculateFloat>
class CudaCoordination : public Colvar {
  const cudaDataType USE_CUDA_SPARSE_PRECISION =
      cudaPrecision<calculateFloat>();
  // constexpr cudaDataType USE_CUDA_SPARSE_PRECISION = precision;
  // std::unique_ptr<NeighborList> nl;
  /// the pointer to the coordinates on the GPU
  CUDAHELPERS::memoryHolder<calculateFloat> cudaPositions;
  /// the pointer to the nn list on the GPU
  CUDAHELPERS::memoryHolder<calculateFloat> cudaCoordination;
  CUDAHELPERS::memoryHolder<calculateFloat> cudaDerivatives;
  CUDAHELPERS::memoryHolder<calculateFloat> cudaVirial;
  CUDAHELPERS::memoryHolder<calculateFloat> reductionMemoryVirial;
  CUDAHELPERS::memoryHolder<calculateFloat> reductionMemoryCoord;
  CUDAHELPERS::memoryHolder<unsigned> cudaTrueIndexes;

  cudaStream_t streamDerivatives;
  cudaStream_t streamVirial;
  cudaStream_t streamCoordination;

  unsigned maxNumThreads = 512;
  // SwitchingFunction switchingFunction;
  rationalSwitchParameters<calculateFloat> switchingParameters;
  ortoPBCs<calculateFloat> myPBC;

  bool pbc{true};
  // bool serial{false};
  // bool invalidateList{true};
  // bool firsttime{true};
  void setUpPermanentGPUMemory();

public:
  explicit CudaCoordination(const ActionOptions &);
  virtual ~CudaCoordination();
  // active methods:
  static void registerKeywords(Keywords &keys);
  // void prepare() override;
  void calculate() override;
};
using CudaCoordination_d = CudaCoordination<double>;
using CudaCoordination_f = CudaCoordination<float>;
PLUMED_REGISTER_ACTION(CudaCoordination_d, "CUDACOORDINATION")
PLUMED_REGISTER_ACTION(CudaCoordination_f, "CUDACOORDINATIONFLOAT")

template <typename calculateFloat>
void CudaCoordination<calculateFloat>::setUpPermanentGPUMemory() {
  auto nat = getPositions().size();
  cudaPositions.resize(3 * nat);
  cudaDerivatives.resize(3 * nat);
  cudaTrueIndexes.resize(nat);
  std::vector<unsigned> trueIndexes(nat);
  for (size_t i = 0; i < nat; ++i) {
    trueIndexes[i] = getAbsoluteIndex(i).index();
  }
  cudaTrueIndexes.copyToCuda(trueIndexes.data(), streamDerivatives);
}

// template <typename calculateFloat>
// void CudaCoordination<calculateFloat>::prepare() {
//   if (nl->getStride() > 0) {
//     if (firsttime || (getStep() % nl->getStride() == 0)) {
//       requestAtoms(nl->getFullAtomList());
//       setUpPermanentGPUMemory();
//       invalidateList = true;
//       firsttime = false;
//     } else {
//       requestAtoms(nl->getReducedAtomList());
//       setUpPermanentGPUMemory();
//       invalidateList = false;
//       if (getExchangeStep())
//         error("Neighbor lists should be updated on exchange steps - choose a
//         "
//               "NL_STRIDE which divides the exchange stride!");
//     }
//     if (getExchangeStep())
//       firsttime = true;
//   }
// }

template <typename calculateFloat>
void CudaCoordination<calculateFloat>::registerKeywords(Keywords &keys) {
  Colvar::registerKeywords(keys);
  // keys.addFlag("SERIAL", false,
  //              "Perform the calculation in serial - for debug purpose");
  // keys.addFlag("PAIR", false,
  //              "Pair only 1st element of the 1st group with 1st "
  //              "element in the second, etc");
  // keys.addFlag("NLIST", false,
  //              "Use a neighbor list to speed up the calculation");
  // keys.add("optional", "NL_CUTOFF", "The cutoff for the neighbor list");
  // keys.add("optional", "NL_STRIDE",
  //          "The frequency with which we are updating the "
  //          "atoms in the neighbor list");
  keys.add("optional", "THREADS", "The upper limit of the number of threads");
  keys.add("atoms", "GROUPA", "First list of atoms");
  // keys.add("atoms", "GROUPB",
  //          "Second list of atoms (if empty, N*(N-1)/2 pairs in "
  //          "GROUPA are counted)");
  keys.add("compulsory", "NN", "6",
           "The n parameter of the switching function ");
  keys.add("compulsory", "MM", "0",
           "The m parameter of the switching function; 0 implies 2*NN");
  keys.add("compulsory", "D_0", "0.0",
           "The d_0 parameter of the switching function");
  keys.add("compulsory", "R_0", "The r_0 parameter of the switching function");
}

// these constant will be used within the kernels
//__constant__ calculateFloat cu_epsilon;
#define cu_epsilon 1e-14
//^This deserves a more intelligent solution

template <typename calculateFloat>
__device__ calculateFloat pcuda_fastpow(calculateFloat base, int expo) {
  if (expo < 0) {
    expo = -expo;
    base = 1.0 / base;
  }
  calculateFloat result = 1.0;
  while (expo) {
    if (expo & 1)
      result *= base;
    expo >>= 1;
    base *= base;
  }
  return result;
}

template <typename calculateFloat>
__device__ calculateFloat pcuda_Rational(const calculateFloat rdist, int NN,
                                         int MM, calculateFloat &dfunc) {
  calculateFloat result;
  if (2 * NN == MM) {
    // if 2*N==M, then (1.0-rdist^N)/(1.0-rdist^M) = 1.0/(1.0+rdist^N)
    calculateFloat rNdist = pcuda_fastpow(rdist, NN - 1);
    calculateFloat iden = 1.0 / (1 + rNdist * rdist);
    dfunc = -NN * rNdist * iden * iden;
    result = iden;
  } else {
    if (rdist > (1. - 100.0 * cu_epsilon) && rdist < (1 + 100.0 * cu_epsilon)) {
      result = NN / MM;
      dfunc = 0.5 * NN * (NN - MM) / MM;
    } else {
      calculateFloat rNdist = pcuda_fastpow(rdist, NN - 1);
      calculateFloat rMdist = pcuda_fastpow(rdist, MM - 1);
      calculateFloat num = 1. - rNdist * rdist;
      calculateFloat iden = 1.0 / (1.0 - rMdist * rdist);
      calculateFloat func = num * iden;
      result = func;
      dfunc = ((-NN * rNdist * iden) + (func * (iden * MM) * rMdist));
    }
  }
  return result;
}

template <typename calculateFloat>
__global__ void getpcuda_Rational(const calculateFloat *rdists, const int NN,
                                  const int MM, calculateFloat *dfunc,
                                  calculateFloat *res) {
  const int i = threadIdx.x + blockIdx.x * blockDim.x;
  if (rdists[i] <= 0.) {
    res[i] = 1.;
    dfunc[i] = 0.0;
  } else
    res[i] = pcuda_Rational(rdists[i], NN, MM, dfunc[i]);
}

// __global__ void getConst() {
//   printf("Cuda: cu_epsilon = %f\n", cu_epsilon);
// }

template <typename calculateFloat>
CudaCoordination<calculateFloat>::CudaCoordination(
    const ActionOptions &ao)
    : PLUMED_COLVAR_INIT(ao) {
  // parseFlag("SERIAL", serial);

  std::vector<AtomNumber> GroupA;
  parseAtomList("GROUPA", GroupA);
  // std::vector<AtomNumber> GroupB;
  // parseAtomList("GROUPB", GroupB);

  bool nopbc = !pbc;
  parseFlag("NOPBC", nopbc);
  pbc = !nopbc;

  // pair stuff
  // bool dopair = false;
  // parseFlag("PAIR", dopair);

  // // neighbor list stuff
  // bool doneigh = false;
  // double nl_cut = 0.0;
  // int nl_st = 0;
  // parseFlag("NLIST", doneigh);
  // if (doneigh) {
  //   parse("NL_CUTOFF", nl_cut);
  //   if (nl_cut <= 0.0)
  //     error("NL_CUTOFF should be explicitly specified and positive");
  //   parse("NL_STRIDE", nl_st);
  //   if (nl_st <= 0)
  //     error("NL_STRIDE should be explicitly specified and positive");
  // }
  parse("THREADS", maxNumThreads);
  if (maxNumThreads <= 0)
    error("THREADS should be positive");
  addValueWithDerivatives();
  setNotPeriodic();
  // if (GroupB.size() > 0) {
  //   if (doneigh) {
  //     nl = Tools::make_unique<NeighborList>(GroupA, GroupB, serial, dopair,
  //     pbc,
  //                                           getPbc(), comm, nl_cut, nl_st);
  //   } else {
  //     nl = Tools::make_unique<NeighborList>(GroupA, GroupB, serial, dopair,
  //     pbc,
  //                                           getPbc(), comm);
  //   }
  // } else {
  //   if (doneigh) {
  //     nl = Tools::make_unique<NeighborList>(GroupA, serial, pbc, getPbc(),
  //     comm,
  //                                           nl_cut, nl_st);
  //   } else {
  //     nl =
  //         Tools::make_unique<NeighborList>(GroupA, serial, pbc, getPbc(),
  //         comm);
  //   }
  // }

  // requestAtoms(nl->getFullAtomList());
  requestAtoms(GroupA);

  // log.printf("  between two groups of %u and %u atoms\n",
  //            static_cast<unsigned>(GroupA.size()),
  //            static_cast<unsigned>(GroupB.size()));
  // log.printf("  first group:\n");
  // for (unsigned int i = 0; i < GroupA.size(); ++i) {
  //   if ((i + 1) % 25 == 0)
  //     log.printf("  \n");
  //   log.printf("  %d", GroupA[i].serial());
  // }
  // log.printf("  \n  second group:\n");
  // for (unsigned int i = 0; i < GroupB.size(); ++i) {
  //   if ((i + 1) % 25 == 0)
  //     log.printf("  \n");
  //   log.printf("  %d", GroupB[i].serial());
  // }
  log.printf("  \n");
  if (pbc)
    log.printf("  using periodic boundary conditions\n");
  else
    log.printf("  without periodic boundary conditions\n");
  // if (dopair)
  //   log.printf("  with PAIR option\n");
  // if (doneigh) {
  //   log.printf("  using neighbor lists with\n");
  //   log.printf("  update every %d steps and cutoff %f\n", nl_st, nl_cut);
  // }
  std::string sw, errors;

  { // loading data to the GPU
    int nn_ = 6;
    int mm_ = 0;
    calculateFloat d0_ = 0.0;
    calculateFloat r0_ = 0.0;
    parse("R_0", r0_);
    if (r0_ <= 0.0) {
      error("R_0 should be explicitly specified and positive");
    }
    parse("D_0", d0_);
    parse("NN", nn_);
    parse("MM", mm_);
    if (mm_ == 0) {
      mm_ = 2 * nn_;
    }

    switchingParameters.nn = nn_;
    switchingParameters.mm = mm_;
    switchingParameters.stretch = 1.0;
    switchingParameters.shift = 0.0;
    calculateFloat dmax = d0_ + r0_ * std::pow(0.00001, 1. / (nn_ - mm_));
    constexpr bool dostretch = true;
    if (dostretch) {
      std::vector<calculateFloat> inputs = {0.0, dmax};
      calculateFloat *inputsc, *dummy;
      calculateFloat *sc;
      cudaMalloc(&inputsc, 2 * sizeof(calculateFloat));
      cudaMalloc(&dummy, 2 * sizeof(calculateFloat));
      cudaMalloc(&sc, 2 * sizeof(calculateFloat));
      cudaMemcpy(inputsc, inputs.data(), 2 * sizeof(calculateFloat),
                 cudaMemcpyHostToDevice);
      getpcuda_Rational<<<1, 2>>>(inputsc, nn_, mm_, dummy, sc);
      std::vector<calculateFloat> s = {0.0, 0.0};
      cudaMemcpy(s.data(), sc, 2 * sizeof(calculateFloat),
                 cudaMemcpyDeviceToHost);
      cudaFree(inputsc);
      cudaFree(dummy);
      cudaFree(sc);
      switchingParameters.stretch = 1.0 / (s[0] - s[1]);
      switchingParameters.shift = -s[1] * switchingParameters.stretch;
    }

    // cudaMemcpyToSymbol(cu_epsilon, &epsilon, sizeof(calculateFloat));
    switchingParameters.dmaxSQ = dmax * dmax;
    calculateFloat invr0 = 1.0 / r0_;
    switchingParameters.invr0_2 = invr0 *= invr0;
  }

  checkRead();
  cudaStreamCreate(&streamDerivatives);
  cudaStreamCreate(&streamVirial);
  cudaStreamCreate(&streamCoordination);
  setUpPermanentGPUMemory();
  // log << "  contacts are counted with cutoff "
  //     << switchingFunction.description() << "\n";
  log << "  contacts are counted with cutoff (dmax)="
      << sqrt(switchingParameters.dmaxSQ)
      << ", with a rational switch with parameters: d0=0.0, r0="
      << 1.0 / sqrt(switchingParameters.invr0_2)
      // switchingParameters.stretch
      // switchingParameters.shift
      << ", N=" << switchingParameters.nn << ", M=" << switchingParameters.mm
      << ".\n";
}

template <typename calculateFloat>
CudaCoordination<calculateFloat>::~CudaCoordination() {
  cudaStreamDestroy(streamDerivatives);
  cudaStreamDestroy(streamVirial);
  cudaStreamDestroy(streamCoordination);
}

template <typename calculateFloat>
__device__ calculateFloat
calculateSqr(const calculateFloat distancesq,
             const rationalSwitchParameters<calculateFloat> switchingParameters,
             calculateFloat &dfunc) {
  calculateFloat result = 0.0;
  dfunc = 0.0;
  if (distancesq < switchingParameters.dmaxSQ) {
    const calculateFloat rdist_2 = distancesq * switchingParameters.invr0_2;
    result = pcuda_Rational(rdist_2, switchingParameters.nn / 2,
                            switchingParameters.mm / 2, dfunc);
    // chain rule:
    dfunc *= 2 * switchingParameters.invr0_2;
    // cu_stretch:
    result = result * switchingParameters.stretch + switchingParameters.shift;
    dfunc *= switchingParameters.stretch;
  }
  return result;
}

#define X(I) 3 * I
#define Y(I) 3 * I + 1
#define Z(I) 3 * I + 2

template <bool usePBC, typename calculateFloat>
__global__ void getCoord(
    const unsigned nat,
    const rationalSwitchParameters<calculateFloat> switchingParameters,
    const ortoPBCs<calculateFloat> myPBC, const calculateFloat *coordinates,
    const unsigned *trueIndexes, calculateFloat *ncoordOut,
    calculateFloat *devOut, calculateFloat *virialOut) {
  // blockDIm are the number of threads in your block
  const unsigned i = threadIdx.x + blockIdx.x * blockDim.x;
  if (i >= nat) // blocks are initializated with 'ceil (nat/threads)'
    return;
  // we try working with less global memory possible
  const unsigned idx = trueIndexes[i];
  calculateFloat mydevX = 0.0;
  calculateFloat mydevY = 0.0;
  calculateFloat mydevZ = 0.0;
  calculateFloat mycoord = 0.0;
  calculateFloat myVirial[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  calculateFloat x = coordinates[X(i)];
  calculateFloat y = coordinates[Y(i)];
  calculateFloat z = coordinates[Z(i)];
  calculateFloat d[3];
  //calculateFloat dsq;
  calculateFloat dfunc;
  calculateFloat coord;
  for (unsigned j = 0; j < nat; ++j) {
    // const unsigned j = threadIdx.y + blockIdx.y * blockDim.y;

    // Safeguard
    if (idx == trueIndexes[j])
      continue;
    // or may be better to set up an
    // const unsigned xyz = threadIdx.z
    // where the third dim is 0 1 2 ^
    if (usePBC) {
      d[0] = pbcClamp((coordinates[X(j)] - x) * myPBC.invX) * myPBC.X;
      d[1] = pbcClamp((coordinates[Y(j)] - y) * myPBC.invY) * myPBC.Y;
      d[2] = pbcClamp((coordinates[Z(j)] - z) * myPBC.invZ) * myPBC.Z;

      // d[0] = remainder(x - coordinates[X(j)], myPBC.X);
      // d[1] = remainder(y - coordinates[Y(j)], myPBC.Y);
      // d[2] = remainder(z - coordinates[Z(j)], myPBC.Z);
    } else {
      d[0] = coordinates[X(j)] - x;
      d[1] = coordinates[Y(j)] - y;
      d[2] = coordinates[Z(j)] - z;
    }

    dfunc = 0.;
    coord = calculateSqr((d[0] * d[0] + d[1] * d[1] + d[2] * d[2]),
     switchingParameters, dfunc);
    mydevX -= dfunc * d[0];
    mydevY -= dfunc * d[1];
    mydevZ -= dfunc * d[2];

    if (i < j) {
      mycoord += coord;
      myVirial[0] -= dfunc * d[0] * d[0];
      myVirial[1] -= dfunc * d[0] * d[1];
      myVirial[2] -= dfunc * d[0] * d[2];
      myVirial[3] -= dfunc * d[1] * d[0];
      myVirial[4] -= dfunc * d[1] * d[1];
      myVirial[5] -= dfunc * d[1] * d[2];
      myVirial[6] -= dfunc * d[2] * d[0];
      myVirial[7] -= dfunc * d[2] * d[1];
      myVirial[8] -= dfunc * d[2] * d[2];
    }
  }
  // working in global memory ONLY at the end
  devOut[X(i)] = mydevX;
  devOut[Y(i)] = mydevY;
  devOut[Z(i)] = mydevZ;
  ncoordOut[i] = mycoord;
  virialOut[nat * 0 + i] = myVirial[0];
  virialOut[nat * 1 + i] = myVirial[1];
  virialOut[nat * 2 + i] = myVirial[2];
  virialOut[nat * 3 + i] = myVirial[3];
  virialOut[nat * 4 + i] = myVirial[4];
  virialOut[nat * 5 + i] = myVirial[5];
  virialOut[nat * 6 + i] = myVirial[6];
  virialOut[nat * 7 + i] = myVirial[7];
  virialOut[nat * 8 + i] = myVirial[8];
}

#define getCoordOrthoPBC getCoord<true>
#define getCoordNoPBC getCoord<true>

template <typename calculateFloat>
void CudaCoordination<calculateFloat>::calculate() {
  auto positions = getPositions();
  auto nat = positions.size();

  Tensor virial;
  double coordination;
  auto deriv = std::vector<Vector>(nat);

  // if (nl->getStride() > 0 && invalidateList)
  //   nl->update(getPositions());

  constexpr unsigned nthreads = 512;

  unsigned ngroups = ceil(double(nat) / nthreads);

  /**********************allocating the memory on the GPU**********************/
  cudaPositions.copyToCuda(&positions[0][0], streamDerivatives);

  cudaCoordination.resize(nat);
  cudaVirial.resize(nat * 9);
  /**************************starting the calculations*************************/
  // this calculates the derivatives and prepare the coordination and the virial
  // for the accumulation
  if (pbc) {
    //Ortho as now
    auto box = getBox();

    myPBC.X = box(0, 0);
    myPBC.Y = box(1, 1);
    myPBC.Z = box(2, 2);
    myPBC.invX = 1.0 / myPBC.X;
    myPBC.invY = 1.0 / myPBC.Y;
    myPBC.invZ = 1.0 / myPBC.X;

    getCoordOrthoPBC<<<ngroups, nthreads, 0, 0>>>(
        nat, switchingParameters, myPBC, cudaPositions.pointer(),
        cudaTrueIndexes.pointer(), cudaCoordination.pointer(),
        cudaDerivatives.pointer(), cudaVirial.pointer());
  } else {
    getCoordNoPBC<<<ngroups, nthreads, 0, 0>>>(
        nat, switchingParameters, myPBC, cudaPositions.pointer(),
        cudaTrueIndexes.pointer(), cudaCoordination.pointer(),
        cudaDerivatives.pointer(), cudaVirial.pointer());
  }

  /**************************accumulating the results**************************/

  auto N = nat;
  bool first = true;
  while (N > 1) {
    size_t runningThreads = CUDAHELPERS::threadsPerBlock(N, maxNumThreads);
    unsigned nGroups = CUDAHELPERS::idealGroups(N, runningThreads);

    reductionMemoryCoord.resize(nGroups);
    reductionMemoryVirial.resize(9 * nGroups);

    if (first) {
      cudaDerivatives.copyFromCuda(&deriv[0][0], streamDerivatives);
      first = false;
    }

    dim3 ngroupsVirial(nGroups, 9);
    CUDAHELPERS::doReductionND(cudaVirial.pointer(),
                               reductionMemoryVirial.pointer(), N,
                               ngroupsVirial, runningThreads, streamVirial);

    CUDAHELPERS::doReduction1D(
        cudaCoordination.pointer(),     // reduceScalarIn->pointer(),
        reductionMemoryCoord.pointer(), // reduceSOut->pointer(),
        N, nGroups, runningThreads, streamCoordination);

    if (nGroups == 1) {
      reductionMemoryVirial.copyFromCuda(&virial[0][0], streamVirial);
      // reduceSOut->copyFromCuda(&coordination,streamCoordination);
      reductionMemoryCoord.copyFromCuda(&coordination, streamCoordination);
    } else {
      // std::swap(reduceScalarIn,reduceSOut);
      reductionMemoryCoord.swap(cudaCoordination);
      reductionMemoryVirial.swap(cudaVirial);
    }

    N = nGroups;
  }
  // in this way we do not resize with additional memory allocation
  if(reductionMemoryCoord.reserved() > cudaCoordination.reserved())
    reductionMemoryCoord.swap(cudaCoordination);
  if(reductionMemoryVirial.reserved() > cudaVirial.reserved())
    reductionMemoryVirial.swap(cudaVirial);
  // this ensures that the memory is fully in the host ram
  cudaDeviceSynchronize();
  for (unsigned i = 0; i < deriv.size(); ++i)
    setAtomsDerivatives(i, deriv[i]);
  

  setValue(coordination);
  setBoxDerivatives(virial);
}
#undef getCoordOrthoPBC
#undef getCoordNoPBC

} // namespace colvar
} // namespace PLMD
