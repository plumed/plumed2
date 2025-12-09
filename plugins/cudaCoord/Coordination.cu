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

#include "plumed/core/ActionRegister.h"
#include "plumed/core/Colvar.h"
#include "plumed/tools/Communicator.h"
#include "plumed/tools/SwitchingFunction.h"
#include "plumed/tools/LinkCells.h"


#include "cudaHelpers.cuh"

#include "Coordination.cuh"

#include <algorithm>
#include <cmath>
#include <cub/block/block_load.cuh>
#include <cub/block/block_reduce.cuh>
#include <memory>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

// cfloat for DLB_EPSILON and FLT_EPSILON
#include <cfloat>

#include <utility>

namespace PLMD {

//this trio of things is used to keep the backward compatibility of this:
//get_data can be used to drop all the weight of setting up the parameter to SwitchingFunction in future iterations
template<typename, typename = void>
constexpr bool has_get_data = false;

template<typename T>
constexpr bool has_get_data <T, std::void_t<decltype(std::declval<T>().get_data())> > = true;

template <typename T>
std::pair<int,int> getNNandMM(T sf, std::string_view NNs, std::string_view MMs) {
  int nn_ = 6;
  int mm_ = 0;
  if constexpr (has_get_data<T>) {
    nn_=sf.get_data().nn;
    mm_=sf.get_data().mm;
  } else {
    Tools::parse(NNs,nn_);
    Tools::parse(MMs,mm_);
    if (mm_ == 0) {
      mm_ = 2 * nn_;
    }
  }
  return {nn_,mm_};
}

namespace colvar {
//+PLUMEDOC COLVAR CUDACOORDINATION
/*
Calculate coordination numbers. Like coordination, but on nvdia gpu and with
limited switching.

CUDACOORDINATION can be invoked with CUDACOORDINATIONFLOAT, but that version
will use single floating point precision, while being faster and compatible with
desktop-based Nvidia cards.

This keyword can be used to calculate the number of contacts between two groups
of atoms and is defined as \f[ \sum_{i\in A} \sum_{i\in B} s_{ij} \f] where
\f$s_{ij}\f$ is 1 if the contact between atoms \f$i\f$ and \f$j\f$ is formed,
zero otherwise.
In actuality, \f$s_{ij}\f$ is replaced with a switching function so as to ensure
that the calculated CV has continuous derivatives. The default switching
function is: \f[ s_{ij} = \frac{ 1 - \left(\frac{{\bf r}_{ij}}{r_0}\right)^n
} { 1 - \left(\frac{{\bf r}_{ij}}{r_0}\right)^m } \f].


\par Examples

Here's an example that shows what happens when providing COORDINATION with
a single group:
\plumedfile
# define some huge group:
group: GROUP ATOMS=1-1000
# Here's coordination within a single group:
CUDACOORDINATION GROUPA=group R_0=0.3

\endplumedfile

*/
//+ENDPLUMEDOC

template <typename precision>
struct worker {
  virtual void runSelf(
    unsigned ngroups,
    unsigned threads,
    unsigned nat,
    unsigned memSize,
    cudaStream_t mystream,
    PLMD::GPU::rationalSwitchParameters<precision> switchingParameters,
    PLMD::GPU::ortoPBCs<precision> myPBC,
    unsigned maxAtomsInCells,
    unsigned const * nat_InCells,
    unsigned const * atomsInCellsIdxs,
    unsigned maxOtherAtoms,
    unsigned const * nat_OtherAtoms,
    unsigned const * otherAtomsIdxs,
    precision const * coordinates,
    unsigned const * trueIndexes,
    precision* ncoordOut,
    precision* devOut,
    precision* virialOut ) =0;

  virtual void runDual(
    unsigned ngroups,
    unsigned threads,
    unsigned nat,
    unsigned memSize,
    cudaStream_t mystream,
    PLMD::GPU::rationalSwitchParameters<precision> switchingParameters,
    PLMD::GPU::ortoPBCs<precision> myPBC,
    unsigned maxAtomsInCells,
    unsigned const * nat_InCells,
    unsigned const * atomsInCellsIdxs,
    unsigned maxOtherAtoms,
    unsigned const * nat_OtherAtoms,
    unsigned const * otherAtomsIdxs,
    precision const * coordinates,
    unsigned const * trueIndexes,
    precision* ncoordOut,
    precision* devOut,
    precision* virialOut ) =0;

  virtual void runDualDev(
    unsigned  ngroups,
    unsigned  threads,
    unsigned  nat,
    unsigned  memSize,
    cudaStream_t mystream,
    PLMD::GPU::rationalSwitchParameters<precision> switchingParameters,
    PLMD::GPU::ortoPBCs<precision> myPBC,
    unsigned  maxAtomsInCells,
    unsigned const * nat_InCells,
    unsigned const * atomsInCellsIdxs,
    unsigned const maxOtherAtoms,
    unsigned const *nat_OtherAtoms,
    unsigned const *otherAtomsIdxs,
    precision const * coordinates,
    unsigned const * trueIndexes,
    precision* devOut)=0;

  virtual void runPair (
    unsigned const ngroups,
    unsigned const threads,
    cudaStream_t mystream,
    unsigned couples,
    unsigned totThreads,
    PLMD::GPU::rationalSwitchParameters<precision> switchingParameters,
    PLMD::GPU::ortoPBCs<precision> myPBC,
    precision const * coordinatesI, // const pointer to const values
    precision const * coordinatesJ, // const pointer to const values
    unsigned const* trueIndexes,
    precision *ncoordOut, // const pointer
    precision *devOutI,
    precision *devOutJ,
    precision *virialOut)=0;

//virtual void runDevDual()=0;
  virtual void getSelfAttr(cudaFuncAttributes& attr)=0;
  virtual void getDualAttr(cudaFuncAttributes& attr)=0;
  virtual void getDualDevAttr(cudaFuncAttributes& attr)=0;
  virtual void getPairAttr(cudaFuncAttributes& attr)=0;
};

#define X(I) 3 * I
#define Y(I) 3 * I + 1
#define Z(I) 3 * I + 2

template <bool usePBC = false, typename T>
T __device__ __forceinline__ calculatePBC (T const val,
    PLMD::GPU::invData<T> const pbc) {
  if constexpr (usePBC) {
    return PLMD::GPU::pbcClamp (val * pbc.inv) * pbc.val;
  } else {
    return val;
  }
}

template <typename calcType,typename calculateFloat>
__global__ void getCoord (const unsigned nat,
                          const unsigned max_inCell,
                          const unsigned max_neigh,
                          const PLMD::GPU::rationalSwitchParameters<calculateFloat> switchingParameters,
                          const PLMD::GPU::ortoPBCs<calculateFloat> myPBC,
                          unsigned const * n_inCell,
                          unsigned const * cellIndexes,
                          unsigned const * n_neigh,
                          unsigned const * neighIndexes,
                          calculateFloat const* coordinates,
                          unsigned const* trueIndexes,
                          calculateFloat* ncoordOut,
                          calculateFloat* devOut,
                          calculateFloat* virialOut) {
  unsigned cellId=blockIdx.x;
  CUDAHELPERS::sharedArena arena;
  auto sPos = arena.get_shared_memory<calculateFloat>(3 * n_neigh[cellId]);
  auto realIndexes = arena.get_shared_memory<unsigned>(n_neigh[cellId]);
  // loading shared memory
  //mocking a View
  auto nnI=neighIndexes+cellId*max_neigh;
  for (auto k = threadIdx.x; k < n_neigh[cellId]; k += blockDim.x) {
    sPos[X (k)] = coordinates[X (nnI[k])];
    sPos[Y (k)] = coordinates[Y (nnI[k])];
    sPos[Z (k)] = coordinates[Z (nnI[k])];
    realIndexes[k] = trueIndexes[nnI[k]];
  }

  __syncthreads();

  if (threadIdx.x >= n_inCell[cellId]) { // blocks are initializated with 'ceil (nat/threads)'
    return;
  }
  const unsigned i = cellIndexes[threadIdx.x + cellId * max_inCell];
  // we try working with less global memory possible, so we set up a bunch of
  // temporary variables
  const unsigned idx = trueIndexes[i];
  // local results
  calculateFloat mydev[3]  = {0.0,0.0,0.0};
  calculateFloat mycoord = 0.0;
//static array became regsters with the compiler (cuda 12)
  calculateFloat myVirial[9] = {0.0, 0.0, 0.0,
                                0.0, 0.0, 0.0,
                                0.0, 0.0, 0.0
                               };
  // local calculation aid
  const calculateFloat xyz[3] = {coordinates[X (i)],
                                 coordinates[Y (i)],
                                 coordinates[Z (i)]
                                };
  for (unsigned j = 0; j < n_neigh[cellId]; ++j) {
    // Safeguard
    if (idx == realIndexes[j] ) {
      continue;
    }
    mycoord+=calcType:: coordLoop(j,
                                  i,
                                  nnI,
                                  xyz,
                                  sPos,
                                  myVirial,
                                  mydev,
                                  myPBC,
                                  switchingParameters);

  }
  // working in global memory ONLY at the end
  devOut[X (i)] = mydev[0];
  devOut[Y (i)] = mydev[1];
  devOut[Z (i)] = mydev[2];
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

template <bool usePBC, typename mySwitch, typename calculateFloat>
__global__ void
getDualDev (
  const unsigned max_inCell,
  const unsigned max_neigh,
  const PLMD::GPU::rationalSwitchParameters<calculateFloat> switchingParameters,
  const PLMD::GPU::ortoPBCs<calculateFloat> myPBC,
  unsigned const * n_inCell,
  unsigned const * cellIndexes,
  unsigned const * n_neigh,
  unsigned const * neighIndexes,
  calculateFloat const* coordinates,
  unsigned const* trueIndexes,
  calculateFloat* devOut
) {
  unsigned cellId=blockIdx.x;
  CUDAHELPERS::sharedArena arena;
  auto sPos = arena.get_shared_memory<calculateFloat>(3 * n_neigh[cellId]);
  auto realIndexes = arena.get_shared_memory<unsigned>(n_neigh[cellId]);
  // loading shared memory
  //mocking a View
  auto nnI=neighIndexes+cellId*max_neigh;
  for (auto k = threadIdx.x; k < n_neigh[cellId]; k += blockDim.x) {
    sPos[X (k)] = coordinates[X (nnI[k])];
    sPos[Y (k)] = coordinates[Y (nnI[k])];
    sPos[Z (k)] = coordinates[Z (nnI[k])];
    realIndexes[k] = trueIndexes[nnI[k]];
  }

  __syncthreads();

  if (threadIdx.x >= n_inCell[cellId]) { // blocks are initializated with 'ceil (nat/threads)'
    return;
  }
  const unsigned i = cellIndexes[threadIdx.x + cellId * max_inCell];
  // we try working with less global memory possible, so we set up a bunch of
  // temporary variables
  const unsigned idx = trueIndexes[i];
  // local results
  calculateFloat mydevX = 0.0;
  calculateFloat mydevY = 0.0;
  calculateFloat mydevZ = 0.0;

  // local calculation aid
  const calculateFloat x = coordinates[X (i)];
  const calculateFloat y = coordinates[Y (i)];
  const calculateFloat z = coordinates[Z (i)];
  calculateFloat d[3];
  calculateFloat t;
  calculateFloat dfunc;

  for (unsigned j = 0; j < n_neigh[cellId]; ++j) {
    // Safeguard
    if (idx == realIndexes[j] ) {
      continue;
    }

    d[0] = calculatePBC<usePBC> (sPos[X(j)] - x, myPBC.X);
    d[1] = calculatePBC<usePBC> (sPos[Y(j)] - y, myPBC.Y);
    d[2] = calculatePBC<usePBC> (sPos[Z(j)] - z, myPBC.Z);

    dfunc = 0.;
    t = GPU::calculateSqr<mySwitch> (
          d[0] * d[0] + d[1] * d[1] + d[2] * d[2],
          switchingParameters,
          dfunc);

    mydevX -= dfunc * d[0];
    mydevY -= dfunc * d[1];
    mydevZ -= dfunc * d[2];
  }
  // working in global memory ONLY at the end
  devOut[X (i)] = mydevX;
  devOut[Y (i)] = mydevY;
  devOut[Z (i)] = mydevZ;
}

constexpr unsigned couplesPerThreads = 8;
template <bool usePBC, typename mySwitch, typename calculateFloat>
__global__ void getCoordPair (
  const unsigned couples,
  const unsigned totThreads,
  const PLMD::GPU::rationalSwitchParameters<calculateFloat> switchingParameters,
  const PLMD::GPU::ortoPBCs<calculateFloat> myPBC,
  const calculateFloat *const coordinatesI, // const pointer to const values
  const calculateFloat *const coordinatesJ, // const pointer to const values
  const unsigned *const trueIndexes,
  calculateFloat *const ncoordOut, // const pointer
  calculateFloat *const devOutI,
  calculateFloat *const devOutJ,
  calculateFloat *const virialOut) {
  // blockDIm are the number of threads in your block
  const unsigned tid = threadIdx.x + blockIdx.x * blockDim.x;
  const unsigned firstCouple = couplesPerThreads * tid;

  if (tid >= totThreads) { // blocks are initializated with 'ceil (nat/threads)'
    return;
  }

  calculateFloat myVirial[9] = {0.0, 0.0, 0.0,
                                0.0, 0.0, 0.0,
                                0.0, 0.0, 0.0
                               };

  calculateFloat d[3];
  // calculateFloat t;
  calculateFloat dfunc;
  calculateFloat t;
  calculateFloat myCoord = 0.0;
  for (unsigned k = 0; k < couplesPerThreads; ++k) {
    const unsigned i = k + firstCouple;

    if (i >= couples) { // blocks are initializated with 'ceil (nat/threads)'
      break;
    }
    // I do not understand why I need to invert
    d[0] = calculatePBC<usePBC> (coordinatesJ[X (i)] - coordinatesI[X (i)],
                                 myPBC.X);
    d[1] = calculatePBC<usePBC> (coordinatesJ[Y (i)] - coordinatesI[Y (i)],
                                 myPBC.Y);
    d[2] = calculatePBC<usePBC> (coordinatesJ[Z (i)] - coordinatesI[Z (i)],
                                 myPBC.Z);
    dfunc = 0.;
    if (trueIndexes[i] != trueIndexes[i + couples]) {
      myCoord += GPU::calculateSqr<mySwitch> (
                   d[0] * d[0] + d[1] * d[1] + d[2] * d[2], switchingParameters, dfunc);
    } else {
      d[0] = d[1] = d[2] = 0.0;
    }

    t = -dfunc * d[0];
    myVirial[0] += t * d[0];
    myVirial[1] += t * d[1];
    myVirial[2] += t * d[2];
    devOutI[X (i)] = t;
    devOutJ[X (i)] = -t;

    t = -dfunc * d[1];
    myVirial[3] += t * d[0];
    myVirial[4] += t * d[1];
    myVirial[5] += t * d[2];
    devOutI[Y (i)] = t;
    devOutJ[Y (i)] = -t;

    t = -dfunc * d[2];
    myVirial[6] += t * d[0];
    myVirial[7] += t * d[1];
    myVirial[8] += t * d[2];
    devOutI[Z (i)] = t;
    devOutJ[Z (i)] = -t;
  }
  ncoordOut[tid] = myCoord;
  virialOut[totThreads * 0 + tid] = myVirial[0];
  virialOut[totThreads * 1 + tid] = myVirial[1];
  virialOut[totThreads * 2 + tid] = myVirial[2];
  virialOut[totThreads * 3 + tid] = myVirial[3];
  virialOut[totThreads * 4 + tid] = myVirial[4];
  virialOut[totThreads * 5 + tid] = myVirial[5];
  virialOut[totThreads * 6 + tid] = myVirial[6];
  virialOut[totThreads * 7 + tid] = myVirial[7];
  virialOut[totThreads * 8 + tid] = myVirial[8];
}

template <bool usePBC, typename mySwitch, typename calculateFloat>
struct selfCalc {
  static __device__ __forceinline__ calculateFloat coordLoop(
    unsigned const j,
    unsigned const i,
    unsigned const * nnI,
    calculateFloat const * xyz,
    calculateFloat const * sPos,
    calculateFloat * myVirial,
    calculateFloat * myDev,
    const PLMD::GPU::ortoPBCs<calculateFloat> myPBC,
    const PLMD::GPU::rationalSwitchParameters<calculateFloat> switchingParameters
  ) {
    calculateFloat d[3],t[3];

    d[0] = calculatePBC<usePBC> (sPos[X(j)] - xyz[0], myPBC.X);
    d[1] = calculatePBC<usePBC> (sPos[Y(j)] - xyz[1], myPBC.Y);
    d[2] = calculatePBC<usePBC> (sPos[Z(j)] - xyz[2], myPBC.Z);

    calculateFloat dfunc = 0.;
    calculateFloat coord = GPU::calculateSqr<mySwitch> (
                             d[0] * d[0] + d[1] * d[1] + d[2] * d[2],
                             switchingParameters,
                             dfunc);

    t[0] = -dfunc * d[0];
    t[1] = -dfunc * d[1];
    t[2] = -dfunc * d[2];
    myDev[0] += t[0];
    myDev[1] += t[1];
    myDev[2] += t[2];
    if (i < nnI[j]) {
      myVirial[0] += t[0] * d[0];
      myVirial[1] += t[0] * d[1];
      myVirial[2] += t[0] * d[2];
      myVirial[3] += t[1] * d[0];
      myVirial[4] += t[1] * d[1];
      myVirial[5] += t[1] * d[2];
      myVirial[6] += t[2] * d[0];
      myVirial[7] += t[2] * d[1];
      myVirial[8] += t[2] * d[2];
      return coord;
    } else {
      return 0.0;
    }
  }
};

template <bool usePBC, typename mySwitch, typename calculateFloat>
struct dualCalc {
  static __device__ __forceinline__ calculateFloat coordLoop(
    unsigned const j,
    unsigned const i,
    unsigned const * /*nnI*/,
    calculateFloat const * xyz,
    calculateFloat const * sPos,
    calculateFloat * myVirial,
    calculateFloat * myDev,
    const PLMD::GPU::ortoPBCs<calculateFloat> myPBC,
    const PLMD::GPU::rationalSwitchParameters<calculateFloat> switchingParameters
  ) {
    calculateFloat d[3],t;

    d[0] = calculatePBC<usePBC> (sPos[X(j)] - xyz[0], myPBC.X);
    d[1] = calculatePBC<usePBC> (sPos[Y(j)] - xyz[1], myPBC.Y);
    d[2] = calculatePBC<usePBC> (sPos[Z(j)] - xyz[2], myPBC.Z);

    calculateFloat dfunc = 0.;
    calculateFloat coord = GPU::calculateSqr<mySwitch> (
                             d[0] * d[0] + d[1] * d[1] + d[2] * d[2],
                             switchingParameters,
                             dfunc);

    t = -dfunc * d[0];
    myDev[0] += t;

    myVirial[0] += t * d[0];
    myVirial[1] += t * d[1];
    myVirial[2] += t * d[2];

    t = -dfunc * d[1];
    myDev[1] += t;

    myVirial[3] += t * d[0];
    myVirial[4] += t * d[1];
    myVirial[5] += t * d[2];

    t = -dfunc * d[2];
    myDev[2] += t;

    myVirial[6] += t * d[0];
    myVirial[7] += t * d[1];
    myVirial[8] += t * d[2];
    return coord;
  }
};

template <bool pbc, typename function, typename precision>
struct runFunction:
  public worker<precision> {
  using self = selfCalc<pbc,function,precision>;
  using dual = dualCalc<pbc,function,precision>;
  void runSelf(
    unsigned const ngroups,
    unsigned const threads,
    unsigned const nat,
    unsigned const memSize,
    cudaStream_t mystream,
    PLMD::GPU::rationalSwitchParameters<precision> switchingParameters,
    PLMD::GPU::ortoPBCs<precision> myPBC,
    unsigned const maxAtomsInCells,
    unsigned const * nat_InCells,
    unsigned const * atomsInCellsIdxs,
    unsigned const maxOtherAtoms,
    unsigned const * nat_OtherAtoms,
    unsigned const * otherAtomsIdxs,
    precision const * coordinates,
    unsigned const * trueIndexes,
    precision* ncoordOut,
    precision* devOut,
    precision* virialOut) override {
    getCoord<self,precision>
    <<<ngroups,threads,
    memSize,
    mystream
    >>> (nat,
         maxAtomsInCells,
         maxOtherAtoms,
         switchingParameters,
         myPBC,
         nat_InCells,
         atomsInCellsIdxs,
         nat_OtherAtoms,
         otherAtomsIdxs,
         coordinates,
         trueIndexes,
         ncoordOut,
         devOut,
         virialOut
        );
  }

  void runDual(
    unsigned const ngroups,
    unsigned const threads,
    unsigned const nat,
    unsigned const memSize,
    cudaStream_t mystream,
    PLMD::GPU::rationalSwitchParameters<precision> switchingParameters,
    PLMD::GPU::ortoPBCs<precision> myPBC,
    unsigned const maxAtomsInCells,
    unsigned const * nat_InCells,
    unsigned const * atomsInCellsIdxs,
    unsigned const maxOtherAtoms,
    unsigned const * nat_OtherAtoms,
    unsigned const * otherAtomsIdxs,
    precision const * coordinates,
    unsigned const * trueIndexes,
    precision* ncoordOut,
    precision* devOut,
    precision* virialOut) override {
    getCoord<dual,precision>
    <<<ngroups,threads,
    memSize,
    mystream
    >>> (nat,
         maxAtomsInCells,
         maxOtherAtoms,
         switchingParameters,
         myPBC,
         nat_InCells,
         atomsInCellsIdxs,
         nat_OtherAtoms,
         otherAtomsIdxs,
         coordinates,
         trueIndexes,
         ncoordOut,
         devOut,
         virialOut
        );
  }

  void runDualDev(
    unsigned const ngroups,
    unsigned const threads,
    unsigned const /*nat*/,
    unsigned const memSize,
    cudaStream_t mystream,
    PLMD::GPU::rationalSwitchParameters<precision> switchingParameters,
    PLMD::GPU::ortoPBCs<precision> myPBC,
    unsigned const maxAtomsInCells,
    unsigned const * nat_InCells,
    unsigned const * atomsInCellsIdxs,
    unsigned const maxOtherAtoms,
    unsigned const * nat_OtherAtoms,
    unsigned const * otherAtomsIdxs,
    precision const * coordinates,
    unsigned const * trueIndexes,
    precision* devOut) override {
    getDualDev<pbc,function,precision>
    <<<ngroups,threads,
    memSize,
    mystream
    >>> (
      maxAtomsInCells,
      maxOtherAtoms,
      switchingParameters,
      myPBC,
      nat_InCells,
      atomsInCellsIdxs,
      nat_OtherAtoms,
      otherAtomsIdxs,
      coordinates,
      trueIndexes,
      devOut
    );
  }

  void runPair(
    unsigned const ngroups,
    unsigned const threads,
    cudaStream_t mystream,
    unsigned const couples,
    unsigned const totThreads,
    PLMD::GPU::rationalSwitchParameters<precision> const switchingParameters,
    PLMD::GPU::ortoPBCs<precision> const myPBC,
    precision const *const coordinatesI,
    precision const *const coordinatesJ,
    unsigned const *const trueIndexes,
    precision *const ncoordOut,
    precision *const devOutI,
    precision *const devOutJ,
    precision *const virialOut) override {
    getCoordPair<pbc,function,precision>
    <<<ngroups, threads,
    0,
    mystream>>> (
      couples,
      totThreads,
      switchingParameters,
      myPBC,
      coordinatesI,
      coordinatesJ,
      trueIndexes,
      ncoordOut,
      devOutI,
      devOutJ,
      virialOut);
  }

  void getSelfAttr(cudaFuncAttributes& attr) override {
    cudaFuncGetAttributes (&attr, &getCoord<self,precision>);
  }
  void getDualAttr(cudaFuncAttributes& attr)override {
    cudaFuncGetAttributes (&attr, &getCoord<dual,precision>);
  }
  void getDualDevAttr(cudaFuncAttributes& attr) override {
    cudaFuncGetAttributes (&attr, &getDualDev<pbc,function,precision>);
  }
  void getPairAttr(cudaFuncAttributes& attr) override {
    cudaFuncGetAttributes (&attr, &getCoordPair<pbc,function,precision>);
  }
};

template <typename calculateFloat>
std::unique_ptr<worker<calculateFloat>> make_worker(bool usePBC,
                                     //for now this is here for deducing the precision
                                     PLMD::GPU::rationalSwitchParameters<calculateFloat> /*switchingParameters*/
) {
  if (usePBC) {
    return std::make_unique<runFunction<true,GPU::Rational,calculateFloat>>();
  } else {
    return std::make_unique<runFunction<false,GPU::Rational,calculateFloat>>();
  }
}

//A way of moving all the informations about the cell setup in the least possible numbers of data movements
class cellSetup {
  struct lims {
    unsigned start=0;
    unsigned size=0;
  };
  lims atomsInCells;
  lims otherAtoms;
  lims nat_InCells;
  lims nat_otherAtoms;
  std::vector<unsigned> data{};
  unsigned maxNeighbors=0;
  unsigned maxAtomsPerCell=0;
public:
  unsigned maxExpected()const {
    return maxNeighbors;
  }
  unsigned biggestCell()const {
    return maxAtomsPerCell;
  }
  View<const unsigned> dataView() const {
    return make_const_view(data);
  }
  void reset(
    const unsigned ncells,
    const unsigned biggest_cell, //the number of atoms in the biggest cell
    const unsigned cellUpper,// the error value for the cell
    const unsigned max_expected,// the maximum number of atoms in the neibourhood
    const unsigned otherUpper//the error value for the neighbourood
  ) {
    maxNeighbors=max_expected;
    maxAtomsPerCell=biggest_cell;
    //this resizes and set to zero the part of the array that stores the number
    // of atoms in cell or in the neigbor list for the cells
    data.assign((maxAtomsPerCell+maxNeighbors+1+1)*ncells,0);
    atomsInCells= {0,ncells*maxAtomsPerCell};
    otherAtoms= {ncells*maxAtomsPerCell,maxNeighbors*ncells};
    nat_InCells= {(maxAtomsPerCell+maxNeighbors)*ncells,ncells};
    nat_otherAtoms= {(maxAtomsPerCell+maxNeighbors+1)*ncells,ncells};
    //we test in production
    plumed_assert(data.size() ==
                  get_atomsInCells().size()
                  +get_otherAtoms().size()
                  +get_nat_InCells().size()
                  +get_nat_otherAtoms().size()) << "the view in cellSetup are not built correcly";

    auto tmp_nat_otherAtoms = get_nat_otherAtoms();
    plumed_assert(&data[data.size()-1] == &tmp_nat_otherAtoms[tmp_nat_otherAtoms.size()-1])
        << "the view may not be correctly set up";
    auto aic=get_atomsInCells();
    std::fill(aic.begin(),aic.end(),cellUpper);
    auto oa=get_otherAtoms();
    std::fill(oa.begin(),oa.end(),otherUpper);
  }
#define getter(tp) PLMD::View<unsigned> get_##tp() {return {data.data()+tp.start,tp.size};}
  getter(atomsInCells)
  getter(otherAtoms)
  getter(nat_InCells)
  getter(nat_otherAtoms)
#undef getter
#define deviceMap(tp) unsigned* deviceMap_##tp(thrust::device_vector<unsigned> &deviceData) const { \
    return thrust::raw_pointer_cast (deviceData.data()) + tp.start; }
  deviceMap(atomsInCells)
  deviceMap(otherAtoms)
  deviceMap(nat_InCells)
  deviceMap(nat_otherAtoms)
#undef deviceMap
};

// companion function for cellSetup
void updateCellists(
  const PLMD::LinkCells& cells,
  const PLMD::LinkCells::CellCollection& listCell,
  const PLMD::LinkCells::CellCollection& listNeigh,
  const unsigned upperCheck,// for checking
  const unsigned maxNNExpected,
  const unsigned biggestCell,
  const bool pbc,
  PLMD::View<unsigned> atomsIncells,
  PLMD::View<unsigned> otherAtoms,
  PLMD::View<unsigned> nat_InCells,
  PLMD::View<unsigned> nat_otherAtoms
) {
  std::vector<unsigned> cells_required(27);
  for(unsigned c =0; c < cells.getNumberOfCells() ; ++c) {
    auto atomsInC= listCell.getCellIndexes(c);
    if(atomsInC.size()>0) {
      plumed_assert(atomsInC.size() <= biggestCell)
          << "Unexpected number of atoms in cell " << c << " (" <<atomsInC.size() <<"!<=" << biggestCell<< ")";
      nat_InCells[c]=atomsInC.size();
      //todo: setup a check to see if the cell needs pbcs (inner)
      auto cell = cells.findMyCell(c);
      unsigned ncells_required=0;
      cells.addRequiredCells(cell,ncells_required, cells_required,pbc);
      std::copy(atomsInC.begin(),atomsInC.end(),atomsIncells.begin()+biggestCell*c);
      unsigned otherAtomsId=c*maxNNExpected;
      for (unsigned cb=0; cb <ncells_required ; ++cb) {
        for (auto B : listNeigh.getCellIndexes(cells_required[cb])) {
          otherAtoms[otherAtomsId]=B;
          ++otherAtomsId;
        }
      }
      nat_otherAtoms[c]=otherAtomsId - c*maxNNExpected;
      if(nat_otherAtoms[c]>0) {
        //no sense in risking to check the address before the start of otherAtoms or to check if no atoms have been added
        plumed_assert(otherAtoms[otherAtomsId-1]<upperCheck)
            << "otherAtoms["<<otherAtomsId-1<<"]"<<"("<<otherAtoms[otherAtomsId-1]<<")"<<"<"<<upperCheck;
      }
      plumed_assert(nat_otherAtoms[c] <= maxNNExpected);
      plumed_assert(nat_InCells[c] <= biggestCell);
    }
  }
}

template <typename calculateFloat>
class CudaCoordination : public Colvar {
  /// the pointer to the coordinates on the GPU
  thrust::device_vector<calculateFloat> cudaPositions;
  /// the pointer to the nn list on the GPU
  thrust::device_vector<calculateFloat> cudaCoordination;
  thrust::device_vector<calculateFloat> cudaDerivatives;
  thrust::device_vector<calculateFloat> cudaVirial;
  thrust::device_vector<calculateFloat> reductionMemoryVirial;
  thrust::device_vector<calculateFloat> reductionMemoryCoord;
  thrust::device_vector<unsigned> cudaTrueIndexes;

  LinkCells cells;
  struct cellsSettings {
    LinkCells::CellCollection listA;
    LinkCells::CellCollection listB;
  } cs;

  cudaStream_t streamDerivatives;
  cudaStream_t streamVirial;
  cudaStream_t streamCoordination;

  unsigned maxNumThreads = 512;
//used in dual mode for the derivative
  unsigned maxNumThreads_dv = 512;
  unsigned maxReductionNumThreads = 512;
  unsigned maxDynamicSharedMemory=0;
  unsigned maxDynamicSharedMemory_dualDev=0;
  unsigned maxNumberOfNeighborsPerCell=0;
  unsigned maxNumberOfNeighborsPerCell_dv=0;
  unsigned atomsInA = 0;
  unsigned atomsInB = 0;
  struct nlSettings {
    calculateFloat cutoff=-1.0;
    int stride=-1;
  } NL;

  std::unique_ptr<worker<calculateFloat>> runner=nullptr;

  PLMD::GPU::rationalSwitchParameters<calculateFloat> switchingParameters;
  PLMD::GPU::ortoPBCs<calculateFloat> myPBC;

  bool pbc{true};
  /// THis is used to not launch 3 times the same thing with MPI
  bool const mpiActive{true};
  void setUpPermanentGPUMemory();

  enum class calculationMode { self, dual, pair, none };
  calculationMode mode = calculationMode::none;
  //used in self and dual mode
  cellSetup cellConfiguration_coord;
  //used in dual mode for the derivatives of group b
  cellSetup cellConfiguration_dv;
  size_t doSelf();
  size_t doDual();
  size_t doPair();

public:
  explicit CudaCoordination (const ActionOptions &);
  virtual ~CudaCoordination();
  // active methods:
  static void registerKeywords (Keywords &keys);
  void calculate() override;
};

using CudaCoordination_d = CudaCoordination<double>;
using CudaCoordination_f = CudaCoordination<float>;
PLUMED_REGISTER_ACTION (CudaCoordination_d, "CUDACOORDINATION")
PLUMED_REGISTER_ACTION (CudaCoordination_f, "CUDACOORDINATIONFLOAT")

template <typename calculateFloat>
void CudaCoordination<calculateFloat>::setUpPermanentGPUMemory() {
  auto nat = getPositions().size();
  cudaPositions.resize (3 * nat);
  cudaDerivatives.resize (3 * nat);
  cudaTrueIndexes.resize (nat);
  std::vector<unsigned> trueIndexes (nat);
  for (size_t i = 0; i < nat; ++i) {
    trueIndexes[i] = getAbsoluteIndex (i).index();
  }
  cudaTrueIndexes = trueIndexes;
}

constexpr auto NNdefault="6";
constexpr auto MMdefault="0";
constexpr auto D_0default="0.0";
constexpr auto D_MAXdefault="-1.0";
template <typename calculateFloat>
void CudaCoordination<calculateFloat>::registerKeywords (Keywords &keys) {
  Colvar::registerKeywords (keys);

  keys.add ("optional", "THREADS", "The upper limit of the number of threads");
  keys.add ("atoms", "GROUPA", "First list of atoms");
  keys.add ("atoms", "GROUPB", "Second list of atoms, optional");
  keys.addFlag ("PAIR",
                false,
                "Pair only 1st element of the 1st group with 1st element in "
                "the second, etc");

  keys.add (
    "compulsory", "NN", NNdefault, "The n parameter of the switching function ");
  keys.add ("compulsory",
            "MM",
            MMdefault,
            "The m parameter of the switching function; 0 implies 2*NN");
  keys.add ("compulsory", "R_0", "The r_0 parameter of the switching function");
  keys.add (
    "compulsory", "D_MAX", D_MAXdefault, "The cut off of the switching function");
  keys.add (
    "compulsory", "D_0", D_0default, "The value of d_0 in the switching function");
  keys.add("compulsory","NL_CUTOFF","-1.0","The cutoff for the neighbor list");
  keys.add("compulsory","NL_STRIDE","-1","The frequency with which we are updating the atoms in the neighbor list");
}

template <typename calculateFloat>
CudaCoordination<calculateFloat>::~CudaCoordination() {
  if(mpiActive) {
    cudaStreamDestroy (streamDerivatives);
    cudaStreamDestroy (streamVirial);
    cudaStreamDestroy (streamCoordination);
  }
}

template <typename calculateFloat>
void CudaCoordination<calculateFloat>::calculate() {
  auto deriv = std::vector<Vector> (getPositions().size());
  if(getStep()%NL.stride==0&& mode != calculationMode::pair) {
    cells.setCutoff(NL.cutoff);
    if (pbc) {
      cells.setupCells(getPbc());
    } else {
      cells.setupCells(getPositions());
    }
    std::vector<unsigned> indexesForCells(getPositions().size());
    std::iota(indexesForCells.begin(),indexesForCells.end(),0);
    const size_t nat = cudaPositions.size() / 3;
    switch (mode) {
    case calculationMode::self: {
      //this needs to be done BY ALL the mpi processes
      cells.resetCollection(cs.listA,
                            make_const_view(getPositions()),
                            make_const_view(indexesForCells));
      if (mpiActive) {
        const auto maxExpected = cs.listA.getMaximimumCombination(27);
        const auto biggestCell = (*std::max_element(cs.listA.lcell_tots.begin(),cs.listA.lcell_tots.end()));

        plumed_assert(maxExpected <= maxNumberOfNeighborsPerCell)
            << "The number of neighbors ("<< maxExpected<<") exceed the current compute capability("
            <<maxNumberOfNeighborsPerCell<<"), try by reducing the cell dimensions";
        plumed_assert(biggestCell <= maxNumThreads)
            << "the number of atoms in a single cell exceds the maximum number of threads avaiable, try by reducing the cell dimensions or by increase the THREADS asked";

        cellConfiguration_coord.reset(cells.getNumberOfCells(),
                                      biggestCell,
                                      nat,
                                      maxExpected,
                                      nat);

        updateCellists(cells,
                       cs.listA,
                       cs.listA,
                       nat,
                       maxExpected,
                       biggestCell,
                       pbc,
                       cellConfiguration_coord.get_atomsInCells(),
                       cellConfiguration_coord.get_otherAtoms(),
                       cellConfiguration_coord.get_nat_InCells(),
                       cellConfiguration_coord.get_nat_otherAtoms());
      }
    }
    break;
    case calculationMode::dual: {
      //this needs to be done BY ALL the mpi processes
      cells.resetCollection(cs.listA,
                            View{getPositions().data(), atomsInA},
                            View<const unsigned> {indexesForCells.data(), atomsInA});
      cells.resetCollection(cs.listB,
                            View{getPositions().data() + atomsInA, atomsInB},
                            View<const unsigned> {indexesForCells.data() + atomsInA, atomsInB});

      if (mpiActive) {
        //the largest cell combination in listA
        const auto maxExpected_dv=cs.listA.getMaximimumCombination(27);
        const auto biggestCell_coord = (*std::max_element(cs.listA.lcell_tots.begin(),cs.listA.lcell_tots.end()));
        plumed_assert(maxExpected_dv <= maxNumberOfNeighborsPerCell_dv)
            << "The number of neighbors ("<< maxExpected_dv<<") exceed the current compute capability("
            << maxNumberOfNeighborsPerCell_dv <<"), try by reducing the cell dimensions";
        plumed_assert(biggestCell_coord <= maxNumThreads)
            << "the number of atoms in a single cell exceds the maximum number of threads avaiable, try by reducing the cell dimensions or by increase the THREADS asked";

        //the largest cell combination in listA
        const auto maxExpected_coord=cs.listB.getMaximimumCombination(27);
        const auto biggestCell_dv = (*std::max_element(cs.listB.lcell_tots.begin(),cs.listB.lcell_tots.end()));
        plumed_assert(maxExpected_coord <= maxNumberOfNeighborsPerCell)
            << "The number of neighbors ("<< maxExpected_coord<<") exceed the current compute capability("
            << maxNumberOfNeighborsPerCell <<"), try by reducing the cell dimensions";
        plumed_assert(biggestCell_dv <= maxNumThreads_dv)
            << "the number of atoms in a single cell exceds the maximum number of threads avaiable, try by reducing the cell dimensions or by increase the THREADS asked";

        cellConfiguration_coord.reset(cells.getNumberOfCells(),
                                      biggestCell_coord,
                                      //the higher index in listA
                                      atomsInA,
                                      maxExpected_coord,
                                      //the higher index in listB
                                      nat);
        cellConfiguration_dv.reset(cells.getNumberOfCells(),
                                   biggestCell_dv,
                                   //the higher index in listB
                                   nat,
                                   maxExpected_dv,
                                   //the higher index in listA
                                   atomsInA);
        updateCellists(cells,
                       cs.listA,
                       cs.listB,
                       nat,
                       maxExpected_coord,
                       biggestCell_coord,
                       pbc,
                       cellConfiguration_coord.get_atomsInCells(),
                       cellConfiguration_coord.get_otherAtoms(),
                       cellConfiguration_coord.get_nat_InCells(),
                       cellConfiguration_coord.get_nat_otherAtoms());

        updateCellists(cells,
                       cs.listB,
                       cs.listA,
                       atomsInA,//the higher index +1 of the first list
                       maxExpected_dv,
                       biggestCell_dv,
                       pbc,
                       cellConfiguration_dv.get_atomsInCells(),
                       cellConfiguration_dv.get_otherAtoms(),
                       cellConfiguration_dv.get_nat_InCells(),
                       cellConfiguration_dv.get_nat_otherAtoms());
      }
    }
    break;
    default:
    {}
    }
  }

  Tensor virial;
  double coordination;
  if(mpiActive) {
    constexpr unsigned dataperthread = 4;
    if (pbc) {
      //myPBC is used in the three calculation mode functions
      // Only ortho as now
      auto box = getBox();

      myPBC.X = box (0, 0);
      myPBC.Y = box (1, 1);
      myPBC.Z = box (2, 2);
      makeWhole();
    }
    auto positions = getPositions();
    /***************************copying data on the GPU**************************/
    CUDAHELPERS::plmdDataToGPU (cudaPositions, positions, streamDerivatives);
    /***************************copying data on the GPU**************************/

    // number of things to be reduced
    size_t t2br = 0;

    switch (mode) {
    case calculationMode::self:
      t2br = doSelf();
      break;
    case calculationMode::dual:
      t2br = doDual();
      break;
    case calculationMode::pair:
      t2br = doPair();
      break;
    case calculationMode::none:
      // throw"this should not have been happened"
      break;
    }

    /**************************accumulating the results**************************/

    cudaDeviceSynchronize();
    CUDAHELPERS::plmdDataFromGPU (cudaDerivatives, deriv, streamDerivatives);

    auto N = t2br;
    // if (N>1){
    //maybe if N=1 it is more efficient to not do the reduction
    //or to sneakily invert the two groups
    do {
      size_t runningThreads = CUDAHELPERS::threadsPerBlock (
                                ceil (double (N) / dataperthread), maxReductionNumThreads);

      unsigned nGroups = ceil (double (N) / (runningThreads * dataperthread));

      reductionMemoryVirial.resize (9 * nGroups);
      reductionMemoryCoord.resize (nGroups);

      dim3 ngroupsVirial (nGroups, 9);
      CUDAHELPERS::doReductionND<dataperthread> (
        thrust::raw_pointer_cast (cudaVirial.data()),
        thrust::raw_pointer_cast (reductionMemoryVirial.data()),
        N,
        ngroupsVirial,
        runningThreads,
        streamVirial);

      CUDAHELPERS::doReduction1D<dataperthread> (
        thrust::raw_pointer_cast (cudaCoordination.data()),
        thrust::raw_pointer_cast (reductionMemoryCoord.data()),
        N,
        nGroups,
        runningThreads,
        streamCoordination);

      if (nGroups == 1) {
        CUDAHELPERS::plmdDataFromGPU (
          reductionMemoryVirial, virial, streamVirial);
        // TODO:find a way to stream this
        coordination = reductionMemoryCoord[0];
      } else {
        reductionMemoryVirial.swap (cudaVirial);
        reductionMemoryCoord.swap (cudaCoordination);
      }

      N = nGroups;
    } while (N > 1);

    // in this way we do not resize with additional memory allocation
    if (reductionMemoryCoord.size() > cudaCoordination.size()) {
      reductionMemoryCoord.swap (cudaCoordination);
    }
    if (reductionMemoryVirial.size() > cudaVirial.size()) {
      reductionMemoryVirial.swap (cudaVirial);
    }
    // this ensures that the memory is fully in the host ram
    cudaDeviceSynchronize();
  }
  comm.Bcast (coordination, 0);
  comm.Bcast (virial, 0);
  comm.Bcast (deriv, 0);
  for (unsigned i = 0; i < deriv.size(); ++i) {
    setAtomsDerivatives (i, deriv[i]);
  }

  setValue (coordination);
  setBoxDerivatives (virial);
}

template <typename calculateFloat>
size_t CudaCoordination<calculateFloat>::doSelf() {
  const size_t nat = cudaPositions.size() / 3;

  /**********************allocating the memory on the GPU**********************/
  cudaCoordination.resize (nat);
  cudaVirial.resize (nat * 9);
  unsigned processedAtoms=0;
  /**************************starting the calculations*************************/
  // this calculates the derivatives and prepare the coordination and the
  // virial for the accumulation

  //Here I am assuming that cs.listA has been updated by cells, if that contract is broke this cannot work
  const auto memSize = cellConfiguration_coord.maxExpected() * (3*sizeof(calculateFloat) + sizeof(unsigned));
  thrust::device_vector<unsigned> deviceData;
  CUDAHELPERS::plmdDataToGPU(deviceData,
                             cellConfiguration_coord.dataView(),
                             streamCoordination
                            );
  const unsigned ngroups = cells.getNumberOfCells();

  //I have alreadyy checked that the biggest cell has less atoms than the limit of threads
  const unsigned threads = cellConfiguration_coord.biggestCell();

  runner->runSelf (
    ngroups,
    threads,
    nat,
    memSize,
    streamCoordination,
    switchingParameters,
    myPBC,
    cellConfiguration_coord.biggestCell(),
    cellConfiguration_coord.deviceMap_nat_InCells(deviceData),
    cellConfiguration_coord.deviceMap_atomsInCells(deviceData),
    cellConfiguration_coord.maxExpected(),
    cellConfiguration_coord.deviceMap_nat_otherAtoms(deviceData),
    cellConfiguration_coord.deviceMap_otherAtoms(deviceData),
    thrust::raw_pointer_cast (cudaPositions.data()),
    thrust::raw_pointer_cast (cudaTrueIndexes.data()),
    thrust::raw_pointer_cast (cudaCoordination.data()),
    thrust::raw_pointer_cast (cudaDerivatives.data()),
    thrust::raw_pointer_cast (cudaVirial.data())
  );
  return nat;
}

template <typename calculateFloat>
size_t CudaCoordination<calculateFloat>::doDual() {
  const size_t nat = cudaPositions.size() / 3;

  //shared memory needed in the derivative loop
  const auto memSize_dv = cellConfiguration_dv.maxExpected() * (3*sizeof(calculateFloat) + sizeof(unsigned));
  //shared memory needed in the coord loop
  const auto memSize_coord = cellConfiguration_coord.maxExpected() * (3*sizeof(calculateFloat) + sizeof(unsigned));
  /**********************allocating the memory on the GPU**********************/
  cudaCoordination.resize (atomsInA);
  cudaVirial.resize (atomsInA * 9);
  /**************************starting the calculations*************************/

  //allocating all the memory
  thrust::device_vector<unsigned> deviceData_coord;
  thrust::device_vector<unsigned> deviceData_dv;

  const unsigned ngroups = cells.getNumberOfCells();
  //I have alreadyy checked that the biggest cell has less atoms than the limit of threads
  const unsigned threads_coord = cellConfiguration_coord.biggestCell();
  const unsigned threads_dv = cellConfiguration_dv.biggestCell();

  //the order here can affect the performance
  //by running kernels concurrently on the gpu
  CUDAHELPERS::plmdDataToGPU(deviceData_coord,
                             cellConfiguration_coord.dataView(),
                             streamCoordination);

  CUDAHELPERS::plmdDataToGPU(deviceData_dv,
                             cellConfiguration_dv.dataView(),
                             streamDerivatives);
  runner->runDual(
    ngroups,
    threads_coord,
    atomsInA,
    memSize_coord,
    streamCoordination,
    switchingParameters,
    myPBC,
    cellConfiguration_coord.biggestCell(),
    cellConfiguration_coord.deviceMap_nat_InCells(deviceData_coord),
    cellConfiguration_coord.deviceMap_atomsInCells(deviceData_coord),
    cellConfiguration_coord.maxExpected(),
    cellConfiguration_coord.deviceMap_nat_otherAtoms(deviceData_coord),
    cellConfiguration_coord.deviceMap_otherAtoms(deviceData_coord),
    thrust::raw_pointer_cast (cudaPositions.data()),
    thrust::raw_pointer_cast (cudaTrueIndexes.data()),
    thrust::raw_pointer_cast (cudaCoordination.data()),
    thrust::raw_pointer_cast (cudaDerivatives.data()),
    thrust::raw_pointer_cast (cudaVirial.data())
  );

  runner->runDualDev(
    ngroups,
    threads_dv,
    atomsInB,//not used
    memSize_dv,
    streamDerivatives,
    switchingParameters,
    myPBC,
    cellConfiguration_dv.biggestCell(),
    cellConfiguration_dv.deviceMap_nat_InCells(deviceData_dv),
    cellConfiguration_dv.deviceMap_atomsInCells(deviceData_dv),
    cellConfiguration_dv.maxExpected(),
    cellConfiguration_dv.deviceMap_nat_otherAtoms(deviceData_dv),
    cellConfiguration_dv.deviceMap_otherAtoms(deviceData_dv),
    thrust::raw_pointer_cast (cudaPositions.data()),
    thrust::raw_pointer_cast (cudaTrueIndexes.data()),
    thrust::raw_pointer_cast (cudaDerivatives.data())
  );
  return atomsInA;
}

template <typename calculateFloat>
size_t CudaCoordination<calculateFloat>::doPair() {
  size_t couples = cudaPositions.size() / 6;
  size_t neededThreads = ceil (double (couples) / couplesPerThreads);
  unsigned ngroups = ceil (double (neededThreads) / maxNumThreads);
  /**********************allocating the memory on the GPU**********************/
  cudaCoordination.resize (neededThreads);
  cudaVirial.resize (neededThreads * 9);
  /**************************starting the calculations*************************/
  runner->runPair(
    ngroups,
    maxNumThreads,
    streamDerivatives,
    couples,
    neededThreads,
    switchingParameters,
    myPBC,
    thrust::raw_pointer_cast (cudaPositions.data()),
    thrust::raw_pointer_cast (cudaPositions.data()) + 3 * couples,
    thrust::raw_pointer_cast (cudaTrueIndexes.data()),
    thrust::raw_pointer_cast (cudaCoordination.data()),
    thrust::raw_pointer_cast (cudaDerivatives.data()),
    thrust::raw_pointer_cast (cudaDerivatives.data()) + 3 * couples,
    thrust::raw_pointer_cast (cudaVirial.data()));
  return neededThreads;
}

template <typename calculateFloat>
CudaCoordination<calculateFloat>::CudaCoordination (const ActionOptions &ao)
  : PLUMED_COLVAR_INIT (ao),
    //mpiActive is const
    mpiActive( comm.Get_rank()== 0),
    cells(comm) {

  std::vector<AtomNumber> GroupA;
  parseAtomList ("GROUPA", GroupA);
  std::vector<AtomNumber> GroupB;
  parseAtomList ("GROUPB", GroupB);

  if (GroupB.size() == 0) {
    mode = calculationMode::self;
    atomsInA = GroupA.size();
  } else {
    mode = calculationMode::dual;
    atomsInA = GroupA.size();
    atomsInB = GroupB.size();
    bool dopair = false;
    parseFlag ("PAIR", dopair);
    if (dopair) {
      if (GroupB.size() == GroupA.size()) {
        mode = calculationMode::pair;
      } else {
        error ("GROUPA and GROUPB must have the same number of atoms if the "
               "PAIR keyword is present");
      }
    }
    GroupA.insert (GroupA.end(), GroupB.begin(), GroupB.end());
  }
  if (mode == calculationMode::none) {
    error (
      R"(implementation error in constructor: calculation mode cannot be "none")");
  }
  bool nopbc = !pbc;
  parseFlag ("NOPBC", nopbc);
  pbc = !nopbc;

  parse ("THREADS", maxNumThreads);
  maxNumThreads_dv =  maxNumThreads;
  if (maxNumThreads <= 0) {
    error ("THREADS should be positive");
  }
  addValueWithDerivatives();
  setNotPeriodic();
  requestAtoms (GroupA);

  log.printf ("  \n");
  if (pbc) {
    log.printf ("  using periodic boundary conditions\n");
  } else {
    log.printf ("  without periodic boundary conditions\n");
  }
  if (mode == calculationMode::pair) {
    log.printf ("  with PAIR option\n");
  }
  std::string sw, errors;

  PLMD::SwitchingFunction sf;

  int nn_ = 6;
  int mm_ = 0;
  calculateFloat dmax;
  calculateFloat d0 = 0.0;
  calculateFloat r0_ = 0.0;
  std::string dmaxs = "-1.0";
  parse ("D_MAX", dmaxs);
  if (dmaxs==D_MAXdefault) {

    parse ("R_0", r0_);
    if (r0_ <= 0.0) {
      error ("R_0 should be explicitly specified and positive");
    }

    parse ("NN", nn_);
    parse ("MM", mm_);
    if (mm_ == 0) {
      mm_ = 2 * nn_;
    }
    parse ("D_0", d0);
    sf.set(nn_,mm_,r0_,d0);
  } else {
    std::string R_0s="";
    parse ("R_0", R_0s);
    if (R_0s=="") {
      error ("R_0 should be explicitly specified and positive");
    }
    std::string switchs="RATIONAL R_0="+R_0s+" D_MAX="+dmaxs;
#define defParse(pn) std::string pn ## s=pn ## default;\
      parse(#pn,pn ## s); \
      if (pn ## s!=pn ## default) {\
        switchs+=" " #pn "="+pn ## s;\
      }
    defParse(NN);
    defParse(MM);
    defParse(D_0);
    std::string errmsg;
    sf.set(switchs,errmsg);
    if (errmsg!="") {
      error(errmsg);
    }
    r0_ = sf.get_r0();
    d0=sf.get_d0();
    std::tie (nn_,mm_) = getNNandMM(sf,NNs,MMs);
  }
  dmax = sf.get_dmax();

  plumed_assert (!(d0<0.0)) << "d0 should be >=0, d0="<<d0;
  switchingParameters.nn = nn_;
  switchingParameters.mm = mm_;
  switchingParameters.stretch = 1.0;
  switchingParameters.shift = 0.0;
  switchingParameters.calcSquared= (! d0 > calculateFloat(0.0) ) && (nn_%2 == 0 && mm_%2 == 0);
  switchingParameters.d0=d0;
  switchingParameters.dmaxSQ=sf.get_dmax2();
  calculateFloat invr0 = 1.0 / r0_;
  if (switchingParameters.calcSquared) {
    switchingParameters.invr0_2 = invr0 * invr0;
  } else {
    switchingParameters.invr0_2 = invr0;
  }
  constexpr bool dostretch = true;
  if (dostretch && mpiActive) {
    std::vector<calculateFloat> inputs = {0.0, (dmax-switchingParameters.d0) * invr0};

    thrust::device_vector<calculateFloat> inputZeroMax = inputs;
    thrust::device_vector<calculateFloat> dummydfunc (2);
    thrust::device_vector<calculateFloat> resZeroMax (2);

    PLMD::GPU::getpcuda_func<PLMD::GPU::Rational><<<1, 2>>> (
      thrust::raw_pointer_cast (inputZeroMax.data()),
      switchingParameters,
      thrust::raw_pointer_cast (dummydfunc.data()),
      thrust::raw_pointer_cast (resZeroMax.data()));

    switchingParameters.stretch = 1.0 / (resZeroMax[0] - resZeroMax[1]);
    switchingParameters.shift = -resZeroMax[1] * switchingParameters.stretch;
  }
  if (switchingParameters.calcSquared) {
    switchingParameters.nn/=2;
    switchingParameters.mm/=2;
  }
  comm.Bcast (switchingParameters.dmaxSQ,0);
  comm.Bcast (switchingParameters.invr0_2,0);
  comm.Bcast (switchingParameters.d0,0);
  comm.Bcast (switchingParameters.stretch,0);
  comm.Bcast (switchingParameters.shift,0);
  comm.Bcast (switchingParameters.nn,0);
  comm.Bcast (switchingParameters.mm,0);

  parse("NL_CUTOFF",NL.cutoff);
  parse("NL_STRIDE",NL.stride);
  if(NL.cutoff < dmax) {
    log << "  (Forcing the cutoff of the neighbor list to be equal to DMAX)\n";
    NL.cutoff = dmax;
    if(NL.stride > 0) {
      error("You must setup manually the cutoff of the neigbor list to be higher than DMAX if you want to setup manually the stride");
    }
  }
  if(NL.stride < 0) {
    log << "  (Forcing the neighbor list to be calculated at each step)\n";
    NL.stride=1;
  }

  checkRead();

  runner = make_worker(pbc,switchingParameters);
  int registers=0;
  int register_dv=0;
  if (mpiActive) {

    cudaStreamCreate (&streamDerivatives);
    cudaStreamCreate (&streamVirial);
    cudaStreamCreate (&streamCoordination);

    setUpPermanentGPUMemory();

    maxReductionNumThreads = min (1024, maxNumThreads);
    {
      cudaFuncAttributes attr;
      // the kernels are heavy on registers, this adjust the maximum number of
      // threads accordingly
      switch (mode) {
      case calculationMode::self:
        runner->getSelfAttr(attr);
        break;
      case calculationMode::dual:
        runner->getDualAttr(attr);
        break;
      case calculationMode::pair:
        runner->getPairAttr(attr);
        break;
      case calculationMode::none:
        plumed_assert(false) << "No calculation mode has been selected, this should not be happened";
        break;
      }
      maxNumThreads = min (attr.maxThreadsPerBlock, maxNumThreads);
      maxDynamicSharedMemory=attr.maxDynamicSharedSizeBytes;
      registers=attr.numRegs;
    }
    if(mode == calculationMode::dual) {
      //getting the attributes for the second loop with two groups
      cudaFuncAttributes attr;
      runner->getDualDevAttr(attr);
      maxNumThreads_dv = min (attr.maxThreadsPerBlock, maxNumThreads_dv);
      maxDynamicSharedMemory_dualDev=attr.maxDynamicSharedSizeBytes;
      maxNumberOfNeighborsPerCell_dv = maxDynamicSharedMemory_dualDev / (3*sizeof(calculateFloat)+sizeof(unsigned));
      register_dv = attr.numRegs;
    }
  }
  maxNumberOfNeighborsPerCell = maxDynamicSharedMemory / (3*sizeof(calculateFloat)+sizeof(unsigned));
  comm.Bcast (maxNumThreads, 0);
  comm.Bcast (maxReductionNumThreads, 0);
  comm.Bcast (maxDynamicSharedMemory, 0);

  log << "  contacts are counted with cutoff (dmax)="
      << sqrt (switchingParameters.dmaxSQ)
      << ", with a rational switch with parameters: "
      "d0="<< switchingParameters.d0
      <<", r0="<< 1.0 / ((switchingParameters.calcSquared)?(sqrt (switchingParameters.invr0_2)):switchingParameters.invr0_2)
      << ", N=" << switchingParameters.nn * ((switchingParameters.calcSquared)?2:1)
      << ", M=" << switchingParameters.mm * ((switchingParameters.calcSquared)?2:1)
      << ".\n";

  log << "Neighbor list info:\n"
      << "\tcutoff: "<< NL.cutoff << "\n"
      << "\tstride: "<< NL.stride << "\n";

  log << "GPU info:\n"
      << "\t max threads per coordination -upper limit on inner cells-: " << maxNumThreads << "\n";
  if (mode == calculationMode::dual) {
    log<< "\t max threads per coordination(derivatives loop): " << maxNumThreads_dv << "\n";
  }

  if (mode == calculationMode::dual || mode == calculationMode::self) {
    log <<"\t max number of atoms in neighbour cells: " << maxNumberOfNeighborsPerCell << "\n";
    if (mode == calculationMode::dual) {
      log <<"\t max number of atoms in neighbour cells(derivatives loop): " << maxNumberOfNeighborsPerCell_dv << "\n";
    }
  }
  log <<"\t registers: " <<registers << "\n";
  if (mode == calculationMode::dual) {
    log <<"\t register(derivatives loop): " << register_dv << "\n";
  }
  log << "\t max threads per reduction: " << maxReductionNumThreads << "\n";
}

} // namespace colvar
} // namespace PLMD
