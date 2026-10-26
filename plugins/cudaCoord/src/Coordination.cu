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

/// the number of scalars accumulated by the kernels:
/// [0] the coordination, [1-9] the virial (row major)
constexpr unsigned nAccumulators = 10;

/// the candidates are loaded in shared memory in tiles of at most this size:
/// keeping it small lets more blocks be resident on each multiprocessor
constexpr unsigned maxTileSize = 512;
/// the number of threads per block used by the cell kernels (a warp per atom):
/// enough warps to load the tiles efficiently, small enough to spread the
/// atoms of a crowded cell over many blocks
constexpr unsigned cellKernelThreads = 256;

/// rounds up to a multiple of the warp size
constexpr unsigned roundUpToWarp(unsigned const n) {
  return ((n + 31u) / 32u) * 32u;
}

/// the number of threads per block for the cell kernels (a warp per atom):
/// cellKernelThreads, unless the cells are small or the user asked for less
/// (maxThreads is a multiple of 32)
constexpr unsigned cellThreads(unsigned const maxAtomsInCell, unsigned const maxThreads) {
  return std::min({roundUpToWarp(32u * std::max(1u, maxAtomsInCell)), cellKernelThreads, maxThreads});
}

/// the number of blocks needed to cover the biggest cell with a warp per atom
constexpr unsigned cellChunks(unsigned const maxAtomsInCell, unsigned const threads) {
  return (maxAtomsInCell + threads / 32u - 1u) / (threads / 32u);
}
/// the chunks go in gridDim.y, that is limited to 65535 blocks
constexpr unsigned maxCellChunks = 65535u;

/// the kernels that compute the coordination write a block-reduced partial
/// sum of the coordination and of the virial in
/// `partials[k*ngroups+blockIdx.x]` (k in [0,nAccumulators) )
/// these partials then need to be reduced with CUDAHELPERS::reducePartials
/// the device view of a set of cell lists, in CSR format:
/// the atoms in the cell `c` are `cellAtoms[cellStart[c]:cellStart[c+1]]` and
/// the atoms in its neighbourhood (the 27 cells around it, itself included)
/// are `neighAtoms[neighStart[c]:neighStart[c+1]]`.
/// Only the non empty cells are stored.
struct cellListsDevice {
  /// number of (non empty) cells
  unsigned nCells=0;
  /// number of atoms in the most populated cell
  unsigned maxAtomsInCell=0;
  /// number of atoms in the most populated neighbourhood
  unsigned maxNeighbours=0;
  unsigned const * cellStart=nullptr;
  unsigned const * cellAtoms=nullptr;
  unsigned const * neighStart=nullptr;
  unsigned const * neighAtoms=nullptr;
};

template <typename precision>
struct worker {
  virtual ~worker() = default;
  // threads must be a multiple of 32, chunks the number of blocks per cell
  virtual void runSelf(
    unsigned threads,
    unsigned chunks,
    cudaStream_t mystream,
    PLMD::GPU::rationalSwitchParameters<precision> switchingParameters,
    PLMD::GPU::pbcData<precision> myPBC,
    cellListsDevice const& lists,
    precision const * coordinates,
    unsigned const * trueIndexes,
    precision* devOut,
    precision* partials ) =0;

  // threads must be a multiple of 32, chunks the number of blocks per cell
  virtual void runDual(
    unsigned threads,
    unsigned chunks,
    cudaStream_t mystream,
    PLMD::GPU::rationalSwitchParameters<precision> switchingParameters,
    PLMD::GPU::pbcData<precision> myPBC,
    cellListsDevice const& lists,
    precision const * coordinates,
    unsigned const * trueIndexes,
    precision* devOut,
    precision* partials ) =0;

  // threads must be a multiple of 32, chunks the number of blocks per cell
  virtual void runDualDev(
    unsigned  threads,
    unsigned  chunks,
    cudaStream_t mystream,
    PLMD::GPU::rationalSwitchParameters<precision> switchingParameters,
    PLMD::GPU::pbcData<precision> myPBC,
    cellListsDevice const& lists,
    precision const * coordinates,
    unsigned const * trueIndexes,
    precision* devOut)=0;

  // threads must be a multiple of 32
  virtual void runPair (
    unsigned const ngroups,
    unsigned const threads,
    cudaStream_t mystream,
    unsigned couples,
    unsigned totThreads,
    PLMD::GPU::rationalSwitchParameters<precision> switchingParameters,
    PLMD::GPU::pbcData<precision> myPBC,
    precision const * coordinatesI, // const pointer to const values
    precision const * coordinatesJ, // const pointer to const values
    unsigned const* trueIndexes,
    precision *devOutI,
    precision *devOutJ,
    precision *partials)=0;

//virtual void runDevDual()=0;
  virtual void getSelfAttr(cudaFuncAttributes& attr)=0;
  virtual void getDualAttr(cudaFuncAttributes& attr)=0;
  virtual void getDualDevAttr(cudaFuncAttributes& attr)=0;
  virtual void getPairAttr(cudaFuncAttributes& attr)=0;
};

#define X(I) 3 * I
#define Y(I) 3 * I + 1
#define Z(I) 3 * I + 2

/// applies the minimum image convention to the distance vector d (in place)
template <bool usePBC = false, typename T>
__device__ __forceinline__ void minimumImage (T (&d)[3],
    PLMD::GPU::pbcData<T> const & pbc) {
  if constexpr (usePBC) {
    // pbc.ortho is a kernel argument: the branch is uniform
    if (pbc.ortho) {
      d[0] = PLMD::GPU::pbcClamp (d[0] * pbc.X.inv) * pbc.X.val;
      d[1] = PLMD::GPU::pbcClamp (d[1] * pbc.Y.inv) * pbc.Y.val;
      d[2] = PLMD::GPU::pbcClamp (d[2] * pbc.Z.inv) * pbc.Z.val;
    } else {
      // scaled coordinates: s = d * invBox, with d a row vector (PLMD convention)
      const T s0 = PLMD::GPU::pbcClamp (d[0] * pbc.invBox[0] + d[1] * pbc.invBox[3] + d[2] * pbc.invBox[6]);
      const T s1 = PLMD::GPU::pbcClamp (d[0] * pbc.invBox[1] + d[1] * pbc.invBox[4] + d[2] * pbc.invBox[7]);
      const T s2 = PLMD::GPU::pbcClamp (d[0] * pbc.invBox[2] + d[1] * pbc.invBox[5] + d[2] * pbc.invBox[8]);
      d[0] = s0 * pbc.box[0] + s1 * pbc.box[3] + s2 * pbc.box[6];
      d[1] = s0 * pbc.box[1] + s1 * pbc.box[4] + s2 * pbc.box[7];
      d[2] = s0 * pbc.box[2] + s1 * pbc.box[5] + s2 * pbc.box[8];
    }
  }
}

/// reduces the 3 components of a vector within the warp (the result is valid in lane 0)
template <typename calculateFloat>
__device__ __forceinline__ void warpReduce3 (calculateFloat (&v)[3]) {
#pragma unroll
  for (unsigned offset = 16; offset > 0; offset /= 2) {
    v[0] += __shfl_down_sync (0xffffffffu, v[0], offset);
    v[1] += __shfl_down_sync (0xffffffffu, v[1], offset);
    v[2] += __shfl_down_sync (0xffffffffu, v[2], offset);
  }
}

// Each block works on a chunk of (blockDim.x/32) atoms of a (non empty) cell:
// blockIdx.x is the cell, blockIdx.y the chunk. Each warp takes one atom of the
// chunk and its 32 lanes split the atoms in the neighbourhood of the cell, that
// are loaded in shared memory in tiles of `tile_size` atoms. In this way the
// work of a single atom is spread over the lanes and cells with many atoms are
// spread over many blocks: the whole GPU is used also with few, crowded cells.
// blockDim.x must be a multiple of 32
template <typename calcType,typename calculateFloat>
__global__ void getCoord (const unsigned tile_size,
                          const PLMD::GPU::rationalSwitchParameters<calculateFloat> switchingParameters,
                          const PLMD::GPU::pbcData<calculateFloat> myPBC,
                          // a kernel argument: a reference would point to host memory
                          // cppcheck-suppress passedByValue
                          const cellListsDevice lists,
                          calculateFloat const* coordinates,
                          unsigned const* trueIndexes,
                          calculateFloat* devOut,
                          calculateFloat* partials) {
  const unsigned cellId=blockIdx.x;
  const unsigned nWarps = blockDim.x / 32u;
  const unsigned warp = threadIdx.x / 32u;
  const unsigned lane = threadIdx.x % 32u;
  CUDAHELPERS::sharedArena arena;
  auto sPos = arena.get_shared_memory<calculateFloat>(3*tile_size);
  auto realIndexes = arena.get_shared_memory<unsigned>(tile_size);

  const unsigned atomsStart = lists.cellStart[cellId];
  const unsigned nInCell = lists.cellStart[cellId+1] - atomsStart;
  const unsigned neighStart = lists.neighStart[cellId];
  const unsigned nNeigh = lists.neighStart[cellId+1] - neighStart;
  unsigned const * neighbours = lists.neighAtoms + neighStart;
  // the atoms of the cell processed by this block
  const unsigned base = blockIdx.y * nWarps;

  // local results, accumulated by each thread over its pairs
  calculateFloat mycoord = 0.0;
//static array became registers with the compiler (cuda 12)
  calculateFloat myVirial[9] = {0.0, 0.0, 0.0,
                                0.0, 0.0, 0.0,
                                0.0, 0.0, 0.0
                               };
  // uniform within the block (the __syncthreads() inside are safe)
  if (base < nInCell) {
    // uniform within the warp
    const bool activeWarp = base + warp < nInCell;
    const unsigned i = (activeWarp) ?
                       lists.cellAtoms[atomsStart + base + warp] : 0;
    const unsigned idx = trueIndexes[i];
    calculateFloat mydev[3]  = {0.0, 0.0, 0.0};
    // local calculation aid
    const calculateFloat xyz[3] = {coordinates[X (i)],
                                   coordinates[Y (i)],
                                   coordinates[Z (i)]
                                  };
    for(unsigned tile=0; tile < nNeigh; tile+=tile_size) {
      // loading shared memory (all the threads of the block)
      const unsigned activeCandidates=min(tile_size,nNeigh-tile);
      const auto nnI=neighbours+tile;
      for (auto k = threadIdx.x; k < activeCandidates; k += blockDim.x) {
        sPos[X (k)] = coordinates[X (nnI[k])];
        sPos[Y (k)] = coordinates[Y (nnI[k])];
        sPos[Z (k)] = coordinates[Z (nnI[k])];
        realIndexes[k] = trueIndexes[nnI[k]];
      }

      __syncthreads();

      if (activeWarp) {
        // the lanes of the warp split the candidates
        for (unsigned j = lane; j < activeCandidates; j += 32u) {
          // Safeguard
          if (idx != realIndexes[j] ) {
            mycoord+=calcType:: coordLoop(j,
                                          idx,
                                          realIndexes,
                                          xyz,
                                          sPos,
                                          myVirial,
                                          mydev,
                                          myPBC,
                                          switchingParameters);
          }
        }
      }
      __syncthreads();
    }
    // the derivative of the atom is spread over the lanes
    warpReduce3 (mydev);
    if (activeWarp && lane == 0) {
      // working in global memory ONLY at the end
      devOut[X (i)] = mydev[0];
      devOut[Y (i)] = mydev[1];
      devOut[Z (i)] = mydev[2];
    }
  }
  // reducing coordination and virial within the block (inactive threads contribute zeros)
  // the shared memory used for the tiles is reused: the last iteration of
  // the tile loop ended with a __syncthreads()
  const calculateFloat acc[nAccumulators] = {
    mycoord,
    myVirial[0], myVirial[1], myVirial[2],
    myVirial[3], myVirial[4], myVirial[5],
    myVirial[6], myVirial[7], myVirial[8]
  };
  CUDAHELPERS::blockReduceToPartials (acc, sPos, partials,
                                      blockIdx.y * gridDim.x + blockIdx.x,
                                      gridDim.x * gridDim.y);
}

// same structure of getCoord, but calculates only the derivatives
// (used for the atoms of the second group in dual mode)
template <bool usePBC, typename mySwitch, typename calculateFloat>
__global__ void
getDualDev (
  const unsigned tile_size,
  const PLMD::GPU::rationalSwitchParameters<calculateFloat> switchingParameters,
  const PLMD::GPU::pbcData<calculateFloat> myPBC,
  // a kernel argument: a reference would point to host memory
  // cppcheck-suppress passedByValue
  const cellListsDevice lists,
  calculateFloat const* coordinates,
  unsigned const* trueIndexes,
  calculateFloat* devOut
) {
  const unsigned cellId=blockIdx.x;
  const unsigned nWarps = blockDim.x / 32u;
  const unsigned warp = threadIdx.x / 32u;
  const unsigned lane = threadIdx.x % 32u;
  CUDAHELPERS::sharedArena arena;
  auto sPos = arena.get_shared_memory<calculateFloat>(3 * tile_size);
  auto realIndexes = arena.get_shared_memory<unsigned>(tile_size);

  const unsigned atomsStart = lists.cellStart[cellId];
  const unsigned nInCell = lists.cellStart[cellId+1] - atomsStart;
  const unsigned neighStart = lists.neighStart[cellId];
  const unsigned nNeigh = lists.neighStart[cellId+1] - neighStart;
  unsigned const * neighbours = lists.neighAtoms + neighStart;
  const unsigned base = blockIdx.y * nWarps;
  if (base >= nInCell) {
    // uniform within the block
    return;
  }
  const bool activeWarp = base + warp < nInCell;
  const unsigned i = (activeWarp) ?
                     lists.cellAtoms[atomsStart + base + warp] : 0;
  const unsigned idx = trueIndexes[i];
  // local results
  calculateFloat mydev[3] = {0.0, 0.0, 0.0};

  // local calculation aid
  const calculateFloat x = coordinates[X (i)];
  const calculateFloat y = coordinates[Y (i)];
  const calculateFloat z = coordinates[Z (i)];

  for(unsigned tile=0; tile < nNeigh; tile+=tile_size) {
    // loading shared memory
    const unsigned activeCandidates=min(tile_size,nNeigh-tile);
    const auto nnI=neighbours+tile;
    for (auto k = threadIdx.x; k < activeCandidates; k += blockDim.x) {
      sPos[X (k)] = coordinates[X (nnI[k])];
      sPos[Y (k)] = coordinates[Y (nnI[k])];
      sPos[Z (k)] = coordinates[Z (nnI[k])];
      realIndexes[k] = trueIndexes[nnI[k]];
    }

    __syncthreads();

    if (activeWarp) {
      calculateFloat d[3];
      calculateFloat dfunc;

      for (unsigned j = lane; j < activeCandidates; j += 32u) {
        // Safeguard
        if (idx != realIndexes[j] ) {

          d[0] = sPos[X(j)] - x;
          d[1] = sPos[Y(j)] - y;
          d[2] = sPos[Z(j)] - z;
          minimumImage<usePBC> (d, myPBC);

          dfunc = 0.;
          GPU::calculateSqr<mySwitch> (
            d[0] * d[0] + d[1] * d[1] + d[2] * d[2],
            switchingParameters,
            dfunc);

          mydev[0] -= dfunc * d[0];
          mydev[1] -= dfunc * d[1];
          mydev[2] -= dfunc * d[2];
        }
      }
    }
    __syncthreads();
  }

  warpReduce3 (mydev);
  if (activeWarp && lane == 0) {
    // working in global memory ONLY at the end
    devOut[X (i)] = mydev[0];
    devOut[Y (i)] = mydev[1];
    devOut[Z (i)] = mydev[2];
  }
}

constexpr unsigned couplesPerThreads = 8;
template <bool usePBC, typename mySwitch, typename calculateFloat>
__global__ void getCoordPair (
  const unsigned couples,
  const unsigned totThreads,
  const PLMD::GPU::rationalSwitchParameters<calculateFloat> switchingParameters,
  const PLMD::GPU::pbcData<calculateFloat> myPBC,
  const calculateFloat *const coordinatesI, // const pointer to const values
  const calculateFloat *const coordinatesJ, // const pointer to const values
  const unsigned *const trueIndexes,
  calculateFloat *const devOutI,
  calculateFloat *const devOutJ,
  calculateFloat *const partials) {
  // blockDIm are the number of threads in your block, must be a multiple of 32
  const unsigned tid = threadIdx.x + blockIdx.x * blockDim.x;
  const unsigned firstCouple = couplesPerThreads * tid;
  // NB: no early return here: all the threads must take part to the final
  // block reduction. Threads with tid >= totThreads have firstCouple >= couples
  // and skip the loop below

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
    d[0] = coordinatesJ[X (i)] - coordinatesI[X (i)];
    d[1] = coordinatesJ[Y (i)] - coordinatesI[Y (i)];
    d[2] = coordinatesJ[Z (i)] - coordinatesI[Z (i)];
    minimumImage<usePBC> (d, myPBC);
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
  // the dynamic shared memory is used only by the reduction
  auto smem = CUDAHELPERS::dynamicSharedMemory<calculateFloat>();
  const calculateFloat acc[nAccumulators] = {
    myCoord,
    myVirial[0], myVirial[1], myVirial[2],
    myVirial[3], myVirial[4], myVirial[5],
    myVirial[6], myVirial[7], myVirial[8]
  };
  CUDAHELPERS::blockReduceToPartials (acc, smem, partials, blockIdx.x, gridDim.x);
}

template <bool usePBC, typename mySwitch, typename calculateFloat>
struct selfCalc {
  /// idx is the true index of the atom and trueIndexes the true indexes of
  /// the candidates in the tile (in shared memory): idx != trueIndexes[j]
  static __device__ __forceinline__ calculateFloat coordLoop(
    unsigned const j,
    unsigned const idx,
    unsigned const * trueIndexes,
    calculateFloat const * xyz,
    calculateFloat const * sPos,
    calculateFloat * myVirial,
    calculateFloat * myDev,
    const PLMD::GPU::pbcData<calculateFloat> myPBC,
    const PLMD::GPU::rationalSwitchParameters<calculateFloat> switchingParameters
  ) {
    calculateFloat d[3],t[3];

    d[0] = sPos[X(j)] - xyz[0];
    d[1] = sPos[Y(j)] - xyz[1];
    d[2] = sPos[Z(j)] - xyz[2];
    minimumImage<usePBC> (d, myPBC);

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
    // each pair is counted once, from the side with the lower true index
    // (equal true indexes are excluded by the caller)
    if (idx < trueIndexes[j]) {
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
    unsigned const /*idx*/,
    unsigned const * /*trueIndexes*/,
    calculateFloat const * xyz,
    calculateFloat const * sPos,
    calculateFloat * myVirial,
    calculateFloat * myDev,
    const PLMD::GPU::pbcData<calculateFloat> myPBC,
    const PLMD::GPU::rationalSwitchParameters<calculateFloat> switchingParameters
  ) {
    calculateFloat d[3],t;

    d[0] = sPos[X(j)] - xyz[0];
    d[1] = sPos[Y(j)] - xyz[1];
    d[2] = sPos[Z(j)] - xyz[2];
    minimumImage<usePBC> (d, myPBC);

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
class runFunction:
  public worker<precision> {
// these are called "maxXperCell" but are used for the tiling settings
  const unsigned maxAtomsPerCell;
  const unsigned maxAtomsPerCell_dual;
  const unsigned maxAtomsPerCell_dualdv;
public:
  using self = selfCalc<pbc,function,precision>;
  using dual = dualCalc<pbc,function,precision>;
  runFunction(): worker<precision>(),
    //getting the tiling informations
    maxAtomsPerCell{[this]{
      cudaFuncAttributes attr;
      // the qualified call makes explicit that the virtual function of this class is used
      runFunction::getSelfAttr(attr);
      auto maxDynamicSharedMemory=attr.maxDynamicSharedSizeBytes;
      return static_cast<unsigned>(maxDynamicSharedMemory / (3*sizeof(precision)+sizeof(unsigned)));
    }()},
  maxAtomsPerCell_dual{[this]{
      cudaFuncAttributes attr;
      // the qualified call makes explicit that the virtual function of this class is used
      runFunction::getDualAttr(attr);
      auto maxDynamicSharedMemory=attr.maxDynamicSharedSizeBytes;
      return static_cast<unsigned>(maxDynamicSharedMemory / (3*sizeof(precision)+sizeof(unsigned)));
    }()},
  maxAtomsPerCell_dualdv {[this]{
      cudaFuncAttributes attr;
      // the qualified call makes explicit that the virtual function of this class is used
      runFunction::getDualDevAttr(attr);
      auto maxDynamicSharedMemory=attr.maxDynamicSharedSizeBytes;
      return static_cast<unsigned>(maxDynamicSharedMemory / (3*sizeof(precision)+sizeof(unsigned)));
    }()}
  {}
  void runSelf(
    unsigned const threads,
    unsigned const chunks,
    cudaStream_t mystream,
    PLMD::GPU::rationalSwitchParameters<precision> switchingParameters,
    PLMD::GPU::pbcData<precision> myPBC,
    cellListsDevice const& lists,
    precision const * coordinates,
    unsigned const * trueIndexes,
    precision* devOut,
    precision* partials) override {
    //how many atoms accumulate in shared memory,
    //so that the kernel will loop on them with no risk of exceeding the maximum memory allowed
    const auto tile_size= std::min({maxAtomsPerCell, lists.maxNeighbours, maxTileSize});
    //the shared memory is reused for the final block reduction
    const auto memSize=std::max(tile_size
                                * (3*sizeof(precision) + sizeof(unsigned)),
                                CUDAHELPERS::blockReduceMemory<nAccumulators,precision>(threads));
    getCoord<self,precision>
    <<<dim3(lists.nCells,chunks),threads,
    memSize,
    mystream
    >>> (tile_size,
         switchingParameters,
         myPBC,
         lists,
         coordinates,
         trueIndexes,
         devOut,
         partials
        );
  }

  void runDual(
    unsigned const threads,
    unsigned const chunks,
    cudaStream_t mystream,
    PLMD::GPU::rationalSwitchParameters<precision> switchingParameters,
    PLMD::GPU::pbcData<precision> myPBC,
    cellListsDevice const& lists,
    precision const * coordinates,
    unsigned const * trueIndexes,
    precision* devOut,
    precision* partials) override {
    //how many atoms accumulate in shared memory,
    //so that the kernel will loop on them with no risk of exceeding the maximum memory allowed
    const auto tile_size= std::min({maxAtomsPerCell_dual, lists.maxNeighbours, maxTileSize});
    //the shared memory is reused for the final block reduction
    const auto memSize=std::max(tile_size
                                * (3*sizeof(precision) + sizeof(unsigned)),
                                CUDAHELPERS::blockReduceMemory<nAccumulators,precision>(threads));
    getCoord<dual,precision>
    <<<dim3(lists.nCells,chunks),threads,
    memSize,
    mystream
    >>> (tile_size,
         switchingParameters,
         myPBC,
         lists,
         coordinates,
         trueIndexes,
         devOut,
         partials
        );
  }

  void runDualDev(
    unsigned const threads,
    unsigned const chunks,
    cudaStream_t mystream,
    PLMD::GPU::rationalSwitchParameters<precision> switchingParameters,
    PLMD::GPU::pbcData<precision> myPBC,
    cellListsDevice const& lists,
    precision const * coordinates,
    unsigned const * trueIndexes,
    precision* devOut) override {
    //how many atoms accumulate in shared memory,
    //so that the kernel will loop on them with no risk of exceeding the maximum memory allowed
    const auto tile_size= std::min({maxAtomsPerCell_dualdv, lists.maxNeighbours, maxTileSize});
    const auto memSize=tile_size
                       * (3*sizeof(precision) + sizeof(unsigned));
    getDualDev<pbc,function,precision>
    <<<dim3(lists.nCells,chunks),threads,
    memSize,
    mystream
    >>> (
      tile_size,
      switchingParameters,
      myPBC,
      lists,
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
    PLMD::GPU::pbcData<precision> const myPBC,
    precision const *const coordinatesI,
    precision const *const coordinatesJ,
    unsigned const *const trueIndexes,
    precision *const devOutI,
    precision *const devOutJ,
    precision *const partials) override {
    getCoordPair<pbc,function,precision>
    <<<ngroups, threads,
    CUDAHELPERS::blockReduceMemory<nAccumulators,precision>(threads),
    mystream>>> (
      couples,
      totThreads,
      switchingParameters,
      myPBC,
      coordinatesI,
      coordinatesJ,
      trueIndexes,
      devOutI,
      devOutJ,
      partials);
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

//The cell lists of a group of atoms, in CSR format (see cellListsDevice).
//
//Only the cells that are non-empty AND have a non-empty neighbourhood are
//stored ("active cells"): the kernels are launched with one block per active
//cell. The data lives in a pinned host buffer and is uploaded on a persistent
//device buffer only when the neighbor list is updated.
//The lists are built with the indexes of the atoms in the action, that can
//then be remapped (see forEachAtomIndex) to the compact indexes of the
//"active atoms", the only ones that are moved to/from the GPU
class cellSetup {
  // layout of data:
  // [cellStart (nCells+1) | neighStart (nCells+1) | cellAtoms (nAtoms) | neighAtoms (nNeigh)]
  CUDAHELPERS::pinnedBuffer<unsigned> data{};
  CUDAHELPERS::deviceBuffer<unsigned> deviceData{};
  // scratch space, kept here to avoid reallocations
  std::vector<unsigned> activeCells{};
  std::vector<unsigned> nRequired{};
  std::vector<unsigned> requiredCells{};
  std::vector<unsigned> cells_required=std::vector<unsigned>(27);
  unsigned nCells=0;
  unsigned nAtoms=0;
  unsigned nNeigh=0;
  unsigned maxAtomsPerCell=0;
  unsigned maxNeighbours=0;
public:
  unsigned activeCellsCount()const {
    return nCells;
  }
  unsigned biggestCell()const {
    return maxAtomsPerCell;
  }
  unsigned maxExpected()const {
    return maxNeighbours;
  }
  ///builds the lists: for each non empty cell of listCell, the atoms in that
  ///cell and the atoms of listNeigh that are in the 27 cells around it
  void build(const PLMD::LinkCells& cells,
             const PLMD::LinkCells::CellCollection& listCell,
             const PLMD::LinkCells::CellCollection& listNeigh,
             const bool pbc) {
    activeCells.clear();
    nRequired.clear();
    requiredCells.clear();
    nAtoms=0;
    nNeigh=0;
    maxAtomsPerCell=0;
    maxNeighbours=0;
    //first pass: counting
    for(unsigned c =0; c < cells.getNumberOfCells() ; ++c) {
      const unsigned n=listCell.lcell_tots[c];
      if(n==0) {
        continue;
      }
      auto cell = cells.findMyCell(c);
      unsigned ncells_required=0;
      cells.addRequiredCells(cell,ncells_required, cells_required,pbc);
      unsigned nn=0;
      for (unsigned cb=0; cb <ncells_required ; ++cb) {
        nn+=listNeigh.lcell_tots[cells_required[cb]];
      }
      if(nn==0) {
        //nothing to compute for the atoms in this cell (their derivatives are zero)
        continue;
      }
      activeCells.push_back(c);
      nRequired.push_back(ncells_required);
      requiredCells.insert(requiredCells.end(),
                           cells_required.begin(),
                           cells_required.begin()+ncells_required);
      nAtoms+=n;
      nNeigh+=nn;
      maxAtomsPerCell=std::max(maxAtomsPerCell,n);
      maxNeighbours=std::max(maxNeighbours,nn);
    }
    nCells=activeCells.size();
    //second pass: filling
    //some slack avoids reallocating at each update of the lists
    data.resize(2*(nCells+1)+nAtoms+nNeigh,1.2);
    unsigned* cellStart=data.data();
    unsigned* neighStart=cellStart+nCells+1;
    unsigned* cellAtoms=neighStart+nCells+1;
    unsigned* neighAtoms=cellAtoms+nAtoms;
    unsigned atomsOffset=0;
    unsigned neighOffset=0;
    unsigned requiredOffset=0;
    for(unsigned k=0; k<nCells; ++k) {
      cellStart[k]=atomsOffset;
      neighStart[k]=neighOffset;
      auto atomsInC= listCell.getCellIndexes(activeCells[k]);
      std::copy(atomsInC.begin(),atomsInC.end(),cellAtoms+atomsOffset);
      atomsOffset+=atomsInC.size();
      for (unsigned cb=0; cb <nRequired[k] ; ++cb) {
        for (auto B : listNeigh.getCellIndexes(requiredCells[requiredOffset+cb])) {
          neighAtoms[neighOffset]=B;
          ++neighOffset;
        }
      }
      requiredOffset+=nRequired[k];
    }
    cellStart[nCells]=atomsOffset;
    neighStart[nCells]=neighOffset;
    //we test in production
    plumed_assert(atomsOffset==nAtoms && neighOffset==nNeigh)
        << "the cell lists are not built correctly";
  }
  ///applies f(unsigned&) to each atom index stored in the lists (the atoms
  ///of the cells first, then the neighbours), in order: it can be used to
  ///remap the indexes in place
  template <typename F>
  void forEachAtomIndex(F&& f) {
    unsigned* cellAtoms=data.data()+2*(nCells+1);
    for (unsigned k=0; k<nAtoms+nNeigh; ++k) {
      f(cellAtoms[k]);
    }
  }
  ///asynchronously copies the data on the GPU
  ///
  ///the host data must not be modified until the stream has been synchronized
  void upload(cudaStream_t stream) {
    deviceData.resize(data.size(),1.2);
    if (data.size()>0) {
      CUDAHELPERS::cudaCheck(cudaMemcpyAsync(deviceData.data(),
                                             data.data(),
                                             data.size()*sizeof(unsigned),
                                             cudaMemcpyHostToDevice,
                                             stream),
                             "uploading the cell lists");
    }
  }
  ///the pointers to the data on the GPU, valid after upload()
  cellListsDevice deviceView() const {
    cellListsDevice v;
    v.nCells=nCells;
    v.maxAtomsInCell=maxAtomsPerCell;
    v.maxNeighbours=maxNeighbours;
    v.cellStart=deviceData.data();
    v.neighStart=v.cellStart+nCells+1;
    v.cellAtoms=v.neighStart+nCells+1;
    v.neighAtoms=v.cellAtoms+nAtoms;
    return v;
  }
};

template <typename calculateFloat>
class CudaCoordination : public Colvar {
  /// the number of atoms
  unsigned nat=0;
  /// the coordinates on the GPU (3*nat)
  CUDAHELPERS::deviceBuffer<calculateFloat> devPositions;
  /// the derivatives (3*nat) followed by the nAccumulators reduced values on the GPU
  CUDAHELPERS::deviceBuffer<calculateFloat> devOutput;
  /// the per-block partial sums of coordination and virial on the GPU
  CUDAHELPERS::deviceBuffer<calculateFloat> devPartials;
  CUDAHELPERS::deviceBuffer<unsigned> devTrueIndexes;
  /// pinned host mirror of devPositions
  CUDAHELPERS::pinnedBuffer<calculateFloat> hostPositions;
  /// pinned host mirror of devOutput
  CUDAHELPERS::pinnedBuffer<calculateFloat> hostOutput;
  /// the derivatives, kept here to avoid reallocating them each step
  std::vector<Vector> deriv;
  std::vector<unsigned> indexesForCells;
  /// the absolute index of each atom of the action
  std::vector<unsigned> trueIndexes;
  /// Only the "active atoms", i.e. the ones that appear in the cell lists,
  /// are moved to/from the GPU: activeAtoms[k] is the index (in the action)
  /// of the k-th active atom, globalToActive its inverse (-1 if not active)
  unsigned nActive=0;
  CUDAHELPERS::pinnedBuffer<unsigned> activeAtoms;
  std::vector<int> globalToActive;
  /// pinned host mirror of devTrueIndexes (for the active atoms)
  CUDAHELPERS::pinnedBuffer<unsigned> hostTrueIndexes;

  LinkCells cells;
  struct cellsSettings {
    LinkCells::CellCollection listA;
    LinkCells::CellCollection listB;
  } cs;

  /// the main stream: non-blocking (does not synchronize with the legacy
  /// default stream, and so with other codes, like GROMACS, that may be using
  /// the same GPU) and with the highest priority
  cudaStream_t stream=nullptr;
  /// an auxiliary stream, used to overlap the two kernels in dual mode
  cudaStream_t streamAux=nullptr;
  cudaEvent_t inputReady=nullptr;
  cudaEvent_t auxDone=nullptr;

  unsigned maxNumThreads = 512;
//used in dual mode for the derivative
  unsigned maxNumThreads_dv = 512;
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
  /// true until the neighbor lists have been built for the first time
  bool firstNLUpdate{true};
  /// if true the lists are also updated when an atom moves more than half the skin
  bool checkDisplacement{false};
  /// (half the skin)^2, where the skin is NL_CUTOFF-D_MAX
  double maxDisplacementSQ{0.0};
  /// the positions at the last update of the neighbor lists
  std::vector<Vector> positionsAtLastUpdate;

  std::unique_ptr<worker<calculateFloat>> runner=nullptr;

  PLMD::GPU::rationalSwitchParameters<calculateFloat> switchingParameters;
  PLMD::GPU::pbcData<calculateFloat> myPBC;

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
  /// assigns compact indexes to the atoms in the lists and remaps the lists
  void remapActiveAtoms();
  /// decides (consistently on all the MPI processes) if the lists must be updated
  bool needsListUpdate();
  void updateNeighborLists();
  /// passes the current box to the kernels (and checks that the minimum image is exact)
  void setupPBC();
  // these return the number of the partial sums to be reduced
  unsigned doSelf();
  unsigned doDual();
  unsigned doPair();

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
  devPositions.resize (3 * nat);
  hostPositions.resize (3 * nat);
  devOutput.resize (3 * nat + nAccumulators);
  hostOutput.resize (3 * nat + nAccumulators);
  devTrueIndexes.resize (nat);
  hostTrueIndexes.resize (nat);
  activeAtoms.resize (nat);
  globalToActive.assign (nat, -1);
  trueIndexes.resize (nat);
  // initially all the atoms are active, with the identity mapping
  // (this stays so in pair mode, self and dual modes remap on each list update)
  nActive = nat;
  for (size_t i = 0; i < nat; ++i) {
    trueIndexes[i] = getAbsoluteIndex (i).index();
    hostTrueIndexes[i] = trueIndexes[i];
    activeAtoms[i] = i;
    globalToActive[i] = i;
  }
  //done only once: a synchronous copy is fine here
  CUDAHELPERS::cudaCheck (cudaMemcpy (devTrueIndexes.data(),
                                      hostTrueIndexes.data(),
                                      nat * sizeof (unsigned),
                                      cudaMemcpyHostToDevice),
                          "copying the atom indexes to the GPU");
}

template <typename calculateFloat>
void CudaCoordination<calculateFloat>::remapActiveAtoms() {
  // the stream must be idle here: the previous step ended with a
  // cudaStreamSynchronize, so the buffers rewritten below are not in use
  std::fill (globalToActive.begin(), globalToActive.end(), -1);
  nActive = 0;
  auto assign = [this] (unsigned & g) {
    int & c = globalToActive[g];
    if (c < 0) {
      c = nActive;
      activeAtoms[nActive] = g;
      hostTrueIndexes[nActive] = trueIndexes[g];
      ++nActive;
    }
    g = c;
  };
  // the atoms of the cells come first, so that the active atoms are sorted by
  // cell: this improves the locality of the loads in the kernels
  cellConfiguration_coord.forEachAtomIndex (assign);
  if (mode == calculationMode::dual) {
    cellConfiguration_dv.forEachAtomIndex (assign);
  }
  if (nActive > 0) {
    CUDAHELPERS::cudaCheck (cudaMemcpyAsync (devTrueIndexes.data(),
                            hostTrueIndexes.data(),
                            nActive * sizeof (unsigned),
                            cudaMemcpyHostToDevice,
                            stream),
                            "copying the atom indexes to the GPU");
  }
  // the derivatives of the atoms that are not active are zero: the active
  // ones are overwritten at each step, so zeroing here is enough
  std::fill (deriv.begin(), deriv.end(), Vector (0.0, 0.0, 0.0));
}

constexpr auto NNdefault="6";
constexpr auto MMdefault="0";
constexpr auto D_0default="0.0";
constexpr auto D_MAXdefault="-1.0";
template <typename calculateFloat>
void CudaCoordination<calculateFloat>::registerKeywords (Keywords &keys) {
  Colvar::registerKeywords (keys);

  keys.add ("optional", "THREADS", "The upper limit of the number of threads per block (default 512, the cell kernels use at most 256)");
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
  keys.add("compulsory","NL_CUTOFF","-1.0","The cutoff for the neighbor list. "
           "If larger than D_MAX the lists are updated automatically when an atom moves more than (NL_CUTOFF-D_MAX)/2");
  keys.add("compulsory","NL_STRIDE","-1","The frequency with which we are updating the atoms in the neighbor list, "
           "if NL_CUTOFF is larger than D_MAX the lists are also updated when needed");
}

template <typename calculateFloat>
CudaCoordination<calculateFloat>::~CudaCoordination() {
  if(mpiActive) {
    //wait for any pending work before releasing the resources
    if (stream) {
      cudaStreamSynchronize (stream);
    }
    if (streamAux) {
      cudaStreamSynchronize (streamAux);
      cudaStreamDestroy (streamAux);
    }
    if (stream) {
      cudaStreamDestroy (stream);
    }
    if (inputReady) {
      cudaEventDestroy (inputReady);
    }
    if (auxDone) {
      cudaEventDestroy (auxDone);
    }
  }
}

template <typename calculateFloat>
void CudaCoordination<calculateFloat>::setupPBC() {
  // no makeWhole() here: the kernels apply the minimum image convention
  // themselves and the cells wrap the coordinates, so it is not needed and
  // it would cost an O(nat) serial pass at every step
  const Pbc & myPbc = getPbc();
  const Tensor & box = myPbc.getBox();
  myPBC.ortho = myPbc.isOrthorombic();
  if (myPBC.ortho) {
    myPBC.X = box (0, 0);
    myPBC.Y = box (1, 1);
    myPBC.Z = box (2, 2);
  } else {
    const Tensor & invBox = myPbc.getInvBox();
    for (unsigned i = 0; i < 3; ++i) {
      for (unsigned j = 0; j < 3; ++j) {
        myPBC.box[3 * i + j] = box (i, j);
        myPBC.invBox[3 * i + j] = invBox (i, j);
      }
    }
    // The minimum image in scaled coordinates is exact for any pair closer
    // than half the smallest height of the box (the distance between two
    // opposite faces), which is what matters since the switching function is
    // zero beyond dmax. For a box with orthogonal vectors (e.g. an
    // orthorhombic box rotated by FIT_TO_TEMPLATE) it is exact for any
    // distance, like Pbc::distance.
    bool orthogonalVectors = true;
    for (unsigned i = 0; i < 3 && orthogonalVectors; ++i) {
      for (unsigned j = i + 1; j < 3; ++j) {
        const Vector a = box.getRow (i);
        const Vector b = box.getRow (j);
        if (std::abs (dotProduct (a, b)) > 1e-6 * a.modulo() * b.modulo()) {
          orthogonalVectors = false;
          break;
        }
      }
    }
    if (!orthogonalVectors) {
      // the heights are the inverse of the norms of the reciprocal vectors,
      // the columns of invBox
      double minHeight = std::numeric_limits<double>::max();
      for (unsigned j = 0; j < 3; ++j) {
        const Vector reciprocal (invBox (0, j), invBox (1, j), invBox (2, j));
        minHeight = std::min (minHeight, 1.0 / reciprocal.modulo());
      }
      const double dmax = std::sqrt (static_cast<double> (switchingParameters.dmaxSQ));
      if (2.0 * dmax > minHeight) {
        error ("with a skewed (non orthogonal) box D_MAX must be smaller than half the smallest "
               "height of the box (" + std::to_string (0.5 * minHeight) + ") for the minimum image "
               "convention to be exact; D_MAX=" + std::to_string (dmax));
      }
    }
  }
}

template <typename calculateFloat>
bool CudaCoordination<calculateFloat>::needsListUpdate() {
  if (mode == calculationMode::pair) {
    return false;
  }
  if (firstNLUpdate) {
    return true;
  }
  int update = (getStep()%NL.stride==0);
  if (!update && checkDisplacement) {
    //same condition on all the MPI processes, so the Bcast is safe
    if (mpiActive) {
      // the raw displacement is used (no pbcs): an atom that has been wrapped
      // in the box by the MD code only triggers a (harmless) extra update
      const auto & positions = getPositions();
      for (unsigned i = 0; i < nat; ++i) {
        if (delta(positionsAtLastUpdate[i],positions[i]).modulo2() > maxDisplacementSQ) {
          update = 1;
          break;
        }
      }
    }
    comm.Bcast (update, 0);
  }
  return update;
}

template <typename calculateFloat>
void CudaCoordination<calculateFloat>::updateNeighborLists() {
  cells.setCutoff(NL.cutoff);
  if (pbc) {
    cells.setupCells(getPbc());
  } else {
    cells.setupCells(getPositions());
  }
  switch (mode) {
  case calculationMode::self: {
    //this needs to be done BY ALL the mpi processes
    cells.resetCollection(cs.listA,
                          make_const_view(getPositions()),
                          make_const_view(indexesForCells));
    if (mpiActive) {
      cellConfiguration_coord.build(cells, cs.listA, cs.listA, pbc);
      remapActiveAtoms();
      cellConfiguration_coord.upload(stream);
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
      //atoms of A, with the atoms of B as neighbours: coordination, virial and derivatives of A
      cellConfiguration_coord.build(cells, cs.listA, cs.listB, pbc);
      //atoms of B, with the atoms of A as neighbours: derivatives of B
      cellConfiguration_dv.build(cells, cs.listB, cs.listA, pbc);
      remapActiveAtoms();
      //both on the main stream: they will be completed before the
      //inputReady event that the auxiliary stream waits for
      cellConfiguration_coord.upload(stream);
      cellConfiguration_dv.upload(stream);
    }
  }
  break;
  default:
  {}
  }
  if (mpiActive && checkDisplacement) {
    const auto & positions = getPositions();
    std::copy(positions.begin(), positions.end(), positionsAtLastUpdate.begin());
  }
}

template <typename calculateFloat>
void CudaCoordination<calculateFloat>::calculate() {
  if (needsListUpdate()) {
    updateNeighborLists();
    firstNLUpdate=false;
  }

  Tensor virial;
  double coordination=0.0;
  if(mpiActive) {
    if (pbc) {
      setupPBC();
    }
    /***************************copying data on the GPU**************************/
    // only the active atoms (the ones in the cell lists) are needed on the GPU
    const auto & positions = getPositions();
    for (unsigned k = 0; k < nActive; ++k) {
      const unsigned i = activeAtoms[k];
      hostPositions[3 * k]     = positions[i][0];
      hostPositions[3 * k + 1] = positions[i][1];
      hostPositions[3 * k + 2] = positions[i][2];
    }
    if (nActive > 0) {
      CUDAHELPERS::cudaCheck (cudaMemcpyAsync (devPositions.data(),
                              hostPositions.data(),
                              3 * nActive * sizeof (calculateFloat),
                              cudaMemcpyHostToDevice,
                              stream),
                              "copying the positions to the GPU");
    }
    /*******************************calculations*******************************/
    // number of partial results to be reduced
    unsigned nPartials = 0;

    switch (mode) {
    case calculationMode::self:
      nPartials = doSelf();
      break;
    case calculationMode::dual:
      nPartials = doDual();
      break;
    case calculationMode::pair:
      nPartials = doPair();
      break;
    case calculationMode::none:
      // throw"this should not have been happened"
      break;
    }

    /**************************accumulating the results**************************/
    // one block per accumulated quantity, the results go right after the derivatives
    CUDAHELPERS::reducePartials<calculateFloat, CUDAHELPERS::reducePartialsThreads>
    <<<nAccumulators, CUDAHELPERS::reducePartialsThreads, 0, stream>>> (devPartials.data(),
                                         nPartials,
                                         devOutput.data() + 3 * nActive);
    CUDAHELPERS::cudaCheck (cudaGetLastError(), "launching the CUDACOORDINATION kernels");
    // a single transfer for derivatives, coordination and virial
    CUDAHELPERS::cudaCheck (cudaMemcpyAsync (hostOutput.data(),
                            devOutput.data(),
                            (3 * nActive + nAccumulators) * sizeof (calculateFloat),
                            cudaMemcpyDeviceToHost,
                            stream),
                            "copying the results from the GPU");
    // the only synchronization point: waits only for the work of this action,
    // not for the other work on the device (like GROMACS' kernels)
    CUDAHELPERS::cudaCheck (cudaStreamSynchronize (stream),
                            "running the CUDACOORDINATION kernels");
    // the derivatives of the non active atoms are zero (set when the lists are updated)
    for (unsigned k = 0; k < nActive; ++k) {
      const unsigned i = activeAtoms[k];
      deriv[i][0] = hostOutput[3 * k];
      deriv[i][1] = hostOutput[3 * k + 1];
      deriv[i][2] = hostOutput[3 * k + 2];
    }
    const calculateFloat *results = hostOutput.data() + 3 * nActive;
    coordination = results[0];
    for (unsigned i = 0; i < 3; ++i) {
      for (unsigned j = 0; j < 3; ++j) {
        virial[i][j] = results[1 + 3 * i + j];
      }
    }
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
unsigned CudaCoordination<calculateFloat>::doSelf() {
  //Here I am assuming that cellConfiguration_coord has been updated and uploaded,
  //if that contract is broke this cannot work
  const auto lists = cellConfiguration_coord.deviceView();
  if (lists.nCells == 0) {
    return 0;
  }
  //a warp per atom: the atoms of a cell are spread over `chunks` blocks
  const unsigned threads = cellThreads (lists.maxAtomsInCell, maxNumThreads);
  const unsigned chunks = cellChunks (lists.maxAtomsInCell, threads);
  if (chunks > maxCellChunks) {
    error ("too many atoms in a single cell (" + std::to_string (lists.maxAtomsInCell)
           + "): increase THREADS or reduce D_MAX/NL_CUTOFF");
  }
  const unsigned nBlocks = lists.nCells * chunks;
  devPartials.resize (nAccumulators * nBlocks, 1.2);

  runner->runSelf (
    threads,
    chunks,
    stream,
    switchingParameters,
    myPBC,
    lists,
    devPositions.data(),
    devTrueIndexes.data(),
    devOutput.data(),
    devPartials.data()
  );
  return nBlocks;
}

template <typename calculateFloat>
unsigned CudaCoordination<calculateFloat>::doDual() {
  const auto lists_coord = cellConfiguration_coord.deviceView();
  const auto lists_dv = cellConfiguration_dv.deviceView();
  //a warp per atom: the atoms of a cell are spread over `chunks` blocks
  const unsigned threads_coord = cellThreads (lists_coord.maxAtomsInCell, maxNumThreads);
  const unsigned chunks_coord = cellChunks (lists_coord.maxAtomsInCell, threads_coord);
  const unsigned threads_dv = cellThreads (lists_dv.maxAtomsInCell, maxNumThreads_dv);
  const unsigned chunks_dv = cellChunks (lists_dv.maxAtomsInCell, threads_dv);
  if (chunks_coord > maxCellChunks || chunks_dv > maxCellChunks) {
    error ("too many atoms in a single cell (" + std::to_string (std::max (lists_coord.maxAtomsInCell, lists_dv.maxAtomsInCell))
           + "): increase THREADS or reduce D_MAX/NL_CUTOFF");
  }
  const unsigned nBlocks = lists_coord.nCells * chunks_coord;
  devPartials.resize (nAccumulators * nBlocks, 1.2);

  // the two kernels write on different parts of the derivatives (A and B atoms)
  // so they can run concurrently: the auxiliary stream waits for the inputs
  // (positions and lists) and then the main stream waits for it
  CUDAHELPERS::cudaCheck (cudaEventRecord (inputReady, stream), "cudaEventRecord");
  CUDAHELPERS::cudaCheck (cudaStreamWaitEvent (streamAux, inputReady, 0), "cudaStreamWaitEvent");
  if (lists_dv.nCells > 0) {
    runner->runDualDev(
      threads_dv,
      chunks_dv,
      streamAux,
      switchingParameters,
      myPBC,
      lists_dv,
      devPositions.data(),
      devTrueIndexes.data(),
      devOutput.data()
    );
  }
  CUDAHELPERS::cudaCheck (cudaEventRecord (auxDone, streamAux), "cudaEventRecord");

  if (lists_coord.nCells > 0) {
    runner->runDual(
      threads_coord,
      chunks_coord,
      stream,
      switchingParameters,
      myPBC,
      lists_coord,
      devPositions.data(),
      devTrueIndexes.data(),
      devOutput.data(),
      devPartials.data()
    );
  }
  CUDAHELPERS::cudaCheck (cudaStreamWaitEvent (stream, auxDone, 0), "cudaStreamWaitEvent");
  return nBlocks;
}

template <typename calculateFloat>
unsigned CudaCoordination<calculateFloat>::doPair() {
  const unsigned couples = nat / 2;
  const unsigned neededThreads = (couples + couplesPerThreads - 1) / couplesPerThreads;
  const unsigned ngroups = (neededThreads + maxNumThreads - 1) / maxNumThreads;
  if (ngroups == 0) {
    return 0;
  }
  devPartials.resize (nAccumulators * ngroups, 1.2);
  runner->runPair(
    ngroups,
    maxNumThreads,
    stream,
    couples,
    neededThreads,
    switchingParameters,
    myPBC,
    devPositions.data(),
    devPositions.data() + 3 * couples,
    devTrueIndexes.data(),
    devOutput.data(),
    devOutput.data() + 3 * couples,
    devPartials.data());
  return ngroups;
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
  if (GroupA.size() == 0) {
    error ("GROUPA must contain at least one atom");
  }

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
  if (maxNumThreads == 0) {
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
  switchingParameters.calcSquared= (d0 == calculateFloat(0.0)) && (nn_%2 == 0 && mm_%2 == 0);
  switchingParameters.d0=d0;
  switchingParameters.dmaxSQ=sf.get_dmax2();
  calculateFloat invr0 = 1.0 / r0_;
  if (switchingParameters.calcSquared) {
    switchingParameters.invr0_2 = invr0 * invr0;
  } else {
    switchingParameters.invr0_2 = invr0;
  }
  // the switching function is stretched so that it is 1 in 0 and 0 in dmax
  if (mpiActive) {
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
  if(NL.stride == 0) {
    error("NL_STRIDE must be positive (or negative to let the lists be updated automatically)");
  }
  if(NL.cutoff < dmax) {
    log << "  (Forcing the cutoff of the neighbor list to be equal to DMAX)\n";
    NL.cutoff = dmax;
  }
  if(!(NL.cutoff > dmax) && NL.stride > 1) {
    // without a skin a list is valid only at the step it was built
    error("You must setup manually the cutoff of the neigbor list to be higher than DMAX if you want to setup manually the stride");
  }
  const double skin = NL.cutoff - dmax;
  if (skin > 0.0) {
    //with a skin the lists are exact as long as no atom moves more than skin/2
    //since the last update: we check that, so that the lists are always valid
    checkDisplacement = true;
    maxDisplacementSQ = 0.25 * skin * skin;
    if (NL.stride < 0) {
      NL.stride = std::numeric_limits<int>::max();
      log << "  (the neighbor list will be updated only when an atom moves more than half the skin, "
          << 0.5 * skin << ")\n";
    } else {
      log << "  (the neighbor list will be updated every " << NL.stride
          << " steps or when an atom moves more than half the skin, " << 0.5 * skin << ")\n";
    }
  } else if(NL.stride < 0) {
    log << "  (Forcing the neighbor list to be calculated at each step)\n";
    NL.stride=1;
  }

  checkRead();

  runner = make_worker(pbc,switchingParameters);
  nat = getPositions().size();
  if (checkDisplacement) {
    positionsAtLastUpdate.resize(nat);
  }
  //these are needed by all the MPI processes
  deriv.resize(nat);
  indexesForCells.resize(nat);
  std::iota(indexesForCells.begin(),indexesForCells.end(),0);
  int registers=0;
  int register_dv=0;
  int streamPriority=0;
  if (mpiActive) {
    //non blocking streams with the highest priority:
    // - they do not synchronize implicitly with the legacy default stream
    // - our (small) kernels get scheduled as soon as possible even if the GPU
    //   is busy with the work of the MD code
    int leastPriority=0;
    CUDAHELPERS::cudaCheck(cudaDeviceGetStreamPriorityRange(&leastPriority,&streamPriority),
                           "cudaDeviceGetStreamPriorityRange");
    CUDAHELPERS::cudaCheck(cudaStreamCreateWithPriority(&stream,cudaStreamNonBlocking,streamPriority),
                           "cudaStreamCreateWithPriority");
    CUDAHELPERS::cudaCheck(cudaStreamCreateWithPriority(&streamAux,cudaStreamNonBlocking,streamPriority),
                           "cudaStreamCreateWithPriority");
    CUDAHELPERS::cudaCheck(cudaEventCreateWithFlags(&inputReady,cudaEventDisableTiming),
                           "cudaEventCreateWithFlags");
    CUDAHELPERS::cudaCheck(cudaEventCreateWithFlags(&auxDone,cudaEventDisableTiming),
                           "cudaEventCreateWithFlags");

    setUpPermanentGPUMemory();

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
      //the kernels that reduce the results need a multiple of the warp size
      maxNumThreads = std::max(32u, (maxNumThreads / 32u) * 32u);
      maxDynamicSharedMemory=attr.maxDynamicSharedSizeBytes;
      registers=attr.numRegs;
    }
    if(mode == calculationMode::dual) {
      //getting the attributes for the second loop with two groups
      cudaFuncAttributes attr;
      runner->getDualDevAttr(attr);
      maxNumThreads_dv = min (attr.maxThreadsPerBlock, maxNumThreads_dv);
      maxNumThreads_dv = std::max(32u, (maxNumThreads_dv / 32u) * 32u);
      maxDynamicSharedMemory_dualDev=attr.maxDynamicSharedSizeBytes;
      maxNumberOfNeighborsPerCell_dv = maxDynamicSharedMemory_dualDev / (3*sizeof(calculateFloat)+sizeof(unsigned));
      register_dv = attr.numRegs;
    }
  }
  maxNumberOfNeighborsPerCell = maxDynamicSharedMemory / (3*sizeof(calculateFloat)+sizeof(unsigned));
  comm.Bcast (maxNumThreads, 0);
  comm.Bcast (maxNumThreads_dv, 0);
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
      << "\tstride: "<< ((NL.stride == std::numeric_limits<int>::max())?
                         std::string("automatic (displacement based)"):
                         std::to_string(NL.stride)) << "\n";
  if (NL.stride == 1 && mode != calculationMode::pair) {
    log << "\t(the neighbor list is rebuilt on the CPU at every step: consider setting"
        " NL_CUTOFF larger than D_MAX and NL_STRIDE>1 for better performance)\n";
  }

  log << "GPU info:\n"
      << "\t threads per block (a warp per atom): " << std::min(cellKernelThreads, maxNumThreads) << "\n";
  if (mode == calculationMode::dual) {
    log<< "\t threads per block (derivatives loop): " << std::min(cellKernelThreads, maxNumThreads_dv) << "\n";
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
  log << "\t stream priority: " << streamPriority << " (non-blocking stream)\n";
}

} // namespace colvar
} // namespace PLMD
