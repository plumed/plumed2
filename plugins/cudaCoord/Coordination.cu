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
#include "plumed/tools/NeighborList.h"
#include "plumed/tools/SwitchingFunction.h"
#include "plumed/tools/Communicator.h"
#include "plumed/tools/LinkCells.h"


#include "cudaHelpers.cuh"
// #include "ndReduction.h"

#include "Coordination.cuh"

#include <algorithm>
#include <cmath>
#include <cub/block/block_load.cuh>
#include <cub/block/block_reduce.cuh>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

// cfloat for DLB_EPSILON and FLT_EPSILON
#include <cfloat>

#include <iostream>
#include <utility>
#include <vector>

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

// does not inherit from coordination base because nl is private
template <typename calculateFloat> class CudaCoordination : public Colvar {
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
  unsigned maxReductionNumThreads = 512;
  unsigned maxDynamiSharedMemory=0;
  unsigned atomsInA = 0;
  unsigned atomsInB = 0;
  struct nlSettings {
    calculateFloat cutoff=-1.0;
    unsigned stride=1;
  } NL;

  PLMD::GPU::rationalSwitchParameters<calculateFloat> switchingParameters;
  PLMD::GPU::ortoPBCs<calculateFloat> myPBC;

  bool pbc{true};
  /// THis is used to not launch 3 times the same thing with MPI
  bool const mpiActive{true};
  void setUpPermanentGPUMemory();

  enum class calculationMode { self, dual, pair, none };
  calculationMode mode = calculationMode::none;
  size_t doSelf();
  void runSelf(
    unsigned ngroups,
    unsigned threads,
    unsigned nat,
    unsigned memSize,
    cudaStream_t mystream,
    unsigned maxAtomsInCells,
    unsigned * const nat_InCells,
    unsigned * const atomsInCellsIdxs,
    unsigned maxOtherAtoms,
    unsigned * const nat_OtherAtoms,
    unsigned * const otherAtomsIdxs
  );
  size_t doDual();
  void runDualCoord(
    unsigned ngroups,
    unsigned threads,
    unsigned nat,
    unsigned memSize,
    cudaStream_t mystream,
    unsigned maxAtomsInCells,
    unsigned * const nat_InCells,
    unsigned * const atomsInCellsIdxs,
    unsigned maxOtherAtoms,
    unsigned * const nat_OtherAtoms,
    unsigned * const otherAtomsIdxs
  );
  void runDualDev(
    unsigned ngroups,
    unsigned threads,
    unsigned nat,
    unsigned memSize,
    cudaStream_t mystream,
    unsigned maxAtomsInCells,
    unsigned * const nat_InCells,
    unsigned * const atomsInCellsIdxs,
    unsigned maxOtherAtoms,
    unsigned * const nat_OtherAtoms,
    unsigned * const otherAtomsIdxs
  );
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
//the enourmous default value  is there becasue there are some problems with  NOPBC
  keys.add("compulsory","NL_CUTOFF","-1.0","The cutoff for the neighbor list");
  keys.add("compulsory","NL_STRIDE","1","The frequency with which we are updating the atoms in the neighbor list");
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
    //cells.setupCells(getPositions(),getPbc());
    std::vector<unsigned> indexesForCells(getPositions().size());
    std::iota(indexesForCells.begin(),indexesForCells.end(),0);
    switch (mode) {
    case calculationMode::self: {
      cells.resetCollection(cs.listA,
                            make_const_view(getPositions()),
                            make_const_view(indexesForCells));
    }
    break;
    case calculationMode::dual: {
      cells.resetCollection(cs.listA,
                            View{getPositions().data(), atomsInA},
                            View<const unsigned> {indexesForCells.data(), atomsInA});
      cells.resetCollection(cs.listB,
                            View{getPositions().data() + atomsInA, atomsInB},
                            View<const unsigned> {indexesForCells.data() + atomsInA, atomsInB});
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
    /*    auto tmp = std::vector<calculateFloat> (cudaCoordination.size());
        cudaMemcpy (tmp.data(),
                    thrust::raw_pointer_cast (cudaCoordination.data()),
                    tmp.size() * sizeof (calculateFloat),
                    cudaMemcpyDeviceToHost);
        std::copy(tmp.begin(),tmp.end(),std::ostream_iterator<calculateFloat>(std::cerr, " "));
        std::cerr <<"\n";
    */

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
    // }else {
    //   CUDAHELPERS::plmdDataFromGPU (
    //       cudaVirial, virial, streamVirial);
    //     // TODO:find a way to stream this
    //     coordination = cudaCoordination[0];
    // }

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

template <bool usePBC, typename calculateFloat>
__global__ void
getSelfCoord (const unsigned nat,
              const unsigned max_inCell,
              const unsigned max_neigh,
              const PLMD::GPU::rationalSwitchParameters<calculateFloat> switchingParameters,
              const PLMD::GPU::ortoPBCs<calculateFloat> myPBC,
              const unsigned* n_inCell,
              const unsigned* cellIndexes,
              const unsigned* n_neigh,
              const unsigned* neighIndexes,
              const calculateFloat* coordinates,
              const unsigned* trueIndexes,
              calculateFloat* ncoordOut,
              calculateFloat* devOut,
              calculateFloat* virialOut) {
  unsigned cellId=blockIdx.x;
  CUDAHELPERS::sharedArena arena;
  auto sPos = arena.get_shared_memory<calculateFloat>(3 * n_neigh[cellId]);
  auto realIndexes = arena.get_shared_memory<unsigned>(n_neigh[cellId]);
  // // loading shared memory
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
  calculateFloat mydevX  = 0.0;
  calculateFloat mydevY  = 0.0;
  calculateFloat mydevZ  = 0.0;
  calculateFloat mycoord = 0.0;
  // the previous version used static array for myVirial and d
  // using explicit variables guarantees that this data will be stored in
  // registers
  calculateFloat myVirial_0 = 0.0;
  calculateFloat myVirial_1 = 0.0;
  calculateFloat myVirial_2 = 0.0;
  calculateFloat myVirial_3 = 0.0;
  calculateFloat myVirial_4 = 0.0;
  calculateFloat myVirial_5 = 0.0;
  calculateFloat myVirial_6 = 0.0;
  calculateFloat myVirial_7 = 0.0;
  calculateFloat myVirial_8 = 0.0;
  // local calculation aid
  const calculateFloat x = coordinates[X (i)];
  const calculateFloat y = coordinates[Y (i)];
  const calculateFloat z = coordinates[Z (i)];
  calculateFloat d_0, d_1, d_2;
  calculateFloat t_0, t_1, t_2;
  calculateFloat dfunc;
  calculateFloat coord;
  for (unsigned j = 0; j < n_neigh[cellId]; ++j) {
    // Safeguard
    if (idx == realIndexes[j] ) {
      continue;
    }

    d_0 = calculatePBC<usePBC> (sPos[X(j)] - x, myPBC.X);
    d_1 = calculatePBC<usePBC> (sPos[Y(j)] - y, myPBC.Y);
    d_2 = calculatePBC<usePBC> (sPos[Z(j)] - z, myPBC.Z);

    dfunc = 0.;
    coord = calculateSqr (
              d_0 * d_0 + d_1 * d_1 + d_2 * d_2,
              switchingParameters,
              dfunc);

    t_0 = -dfunc * d_0;
    t_1 = -dfunc * d_1;
    t_2 = -dfunc * d_2;
    mydevX += t_0;
    mydevY += t_1;
    mydevZ += t_2;
    if (i < nnI[j]) {
      mycoord += coord;
      myVirial_0 += t_0 * d_0;
      myVirial_1 += t_0 * d_1;
      myVirial_2 += t_0 * d_2;
      myVirial_3 += t_1 * d_0;
      myVirial_4 += t_1 * d_1;
      myVirial_5 += t_1 * d_2;
      myVirial_6 += t_2 * d_0;
      myVirial_7 += t_2 * d_1;
      myVirial_8 += t_2 * d_2;
    }
  }
  // working in global memory ONLY at the end
  devOut[X (i)] = mydevX;
  devOut[Y (i)] = mydevY;
  devOut[Z (i)] = mydevZ;
  ncoordOut[i] = mycoord;
  virialOut[nat * 0 + i] = myVirial_0;
  virialOut[nat * 1 + i] = myVirial_1;
  virialOut[nat * 2 + i] = myVirial_2;
  virialOut[nat * 3 + i] = myVirial_3;
  virialOut[nat * 4 + i] = myVirial_4;
  virialOut[nat * 5 + i] = myVirial_5;
  virialOut[nat * 6 + i] = myVirial_6;
  virialOut[nat * 7 + i] = myVirial_7;
  virialOut[nat * 8 + i] = myVirial_8;
}

template <typename calculateFloat>
void CudaCoordination<calculateFloat>::runSelf(
  unsigned const ngroups,
  unsigned const threads,
  unsigned const nat,
  unsigned const memSize,
  cudaStream_t /*mystream*/,
  unsigned const maxAtomsInCells,
  unsigned * const nat_InCells,
  unsigned * const atomsInCellsIdxs,
  unsigned const maxOtherAtoms,
  unsigned * const nat_OtherAtoms,
  unsigned * const otherAtomsIdxs
) {
  if(pbc) {
    getSelfCoord<true>
    <<<ngroups,
    threads,
    memSize,
    0//        streamDerivatives
    >>> (
      nat,
      maxAtomsInCells,
      maxOtherAtoms,
      switchingParameters,
      myPBC,
      nat_InCells,
      atomsInCellsIdxs,
      nat_OtherAtoms,
      otherAtomsIdxs,
      thrust::raw_pointer_cast (cudaPositions.data()),
      thrust::raw_pointer_cast (cudaTrueIndexes.data()),
      thrust::raw_pointer_cast (cudaCoordination.data()),
      thrust::raw_pointer_cast (cudaDerivatives.data()),
      thrust::raw_pointer_cast (cudaVirial.data())
    );
  } else {
    getSelfCoord<false>
    <<<ngroups,
    threads,
    memSize,
    0//        streamDerivatives
    >>> (
      nat,
      maxAtomsInCells,
      maxOtherAtoms,
      switchingParameters,
      myPBC,
      nat_InCells,
      atomsInCellsIdxs,
      nat_OtherAtoms,
      otherAtomsIdxs,
      thrust::raw_pointer_cast (cudaPositions.data()),
      thrust::raw_pointer_cast (cudaTrueIndexes.data()),
      thrust::raw_pointer_cast (cudaCoordination.data()),
      thrust::raw_pointer_cast (cudaDerivatives.data()),
      thrust::raw_pointer_cast (cudaVirial.data())
    );
  }
}

void updateCellists(
  const PLMD::LinkCells& cells,
  const PLMD::LinkCells::CellCollection& listCell,
  const PLMD::LinkCells::CellCollection& listNeigh,
  const unsigned nat,// for checking
  const unsigned maxExpected,
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
      nat_InCells[c]=atomsInC.size();
      //todo: setup a check to see if the cell needs pbcs (inner)
      auto cell = cells.findMyCell(c);
      unsigned ncells_required=0;
      cells.addRequiredCells(cell,ncells_required, cells_required,pbc);
      std::copy(atomsInC.begin(),atomsInC.end(),atomsIncells.begin()+biggestCell*c);
      unsigned otherAtomsId=c*maxExpected;
      for (unsigned cb=0; cb <ncells_required ; ++cb) {
        for (auto B : listNeigh.getCellIndexes(cells_required[cb])) {
          otherAtoms[otherAtomsId]=B;
          ++otherAtomsId;
        }
      }
      nat_otherAtoms[c]=otherAtomsId - c*maxExpected;
      if(nat_otherAtoms[c]>0) {
        //no sense in risking to check the address before the start of otherAtoms or to check if no atoms have been added
        plumed_assert(otherAtoms[otherAtomsId-1]<nat) << "otherAtoms["<<otherAtomsId-1<<"]"<<"("<<otherAtoms[otherAtomsId-1]<<")"<<"<"<<nat;
      }
      plumed_assert(nat_otherAtoms[c] <= maxExpected);
      plumed_assert(nat_InCells[c] <= biggestCell);
    }
  }
}

class cellSetup {
  struct lims {
    unsigned start;
    unsigned size;
  };
  lims atomsInCells;
  lims otherAtoms;
  lims nat_InCells;
  lims nat_otherAtoms;
public:
  std::vector<unsigned> data;
  cellSetup(
    const unsigned ncells,
    const unsigned biggestCell, //the number of atoms in the biggest cell
    const unsigned cellUpper,// the error value for the cell
    const unsigned maxExpected,// the maximum number of atoms in the neibourhood
    const unsigned otherUpper//the error value for the neighbourood
  ):data((biggestCell+maxExpected+1+1)*ncells,0),
    atomsInCells{0,ncells*biggestCell},
    otherAtoms{ncells*biggestCell,maxExpected*ncells},
    nat_InCells{(biggestCell+maxExpected)*ncells,ncells},
    nat_otherAtoms{(biggestCell+maxExpected+1)*ncells,ncells} {
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
    std::fill(oa.begin(),oa.end(),maxExpected);
  }
#define getter(tp) PLMD::View<unsigned> get_##tp(){return {data.data()+tp.start,tp.size};}
  getter(atomsInCells)
  getter(otherAtoms)
  getter(nat_InCells)
  getter(nat_otherAtoms)
#undef getter
#define deviceMap(tp) unsigned* deviceMap_##tp(thrust::device_vector<unsigned> &deviceData){ \
    return thrust::raw_pointer_cast (deviceData.data()) + tp.start; }
  deviceMap(atomsInCells)
  deviceMap(otherAtoms)
  deviceMap(nat_InCells)
  deviceMap(nat_otherAtoms)
#undef deviceMap
};



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
  std::vector<unsigned> cells_required(27);
  const auto maxExpected = cs.listA.getMaximimumCombination(27);
  const auto memSize = maxExpected * (3*sizeof(calculateFloat) + sizeof(unsigned));
  const auto biggestCell = (*std::max_element(cs.listA.lcell_tots.begin(),cs.listA.lcell_tots.end()));

  cellSetup cellConfiguration(cells.getNumberOfCells(),biggestCell,nat,maxExpected,nat);

  updateCellists(cells,
                 cs.listA,
                 cs.listA,
                 nat,
                 maxExpected,
                 biggestCell,
                 pbc,
                 cellConfiguration.get_atomsInCells(),
                 cellConfiguration.get_otherAtoms(),
                 cellConfiguration.get_nat_InCells(),
                 cellConfiguration.get_nat_otherAtoms());

  thrust::device_vector<unsigned> deviceData;
  CUDAHELPERS::plmdDataToGPU(deviceData,
                             make_view(cellConfiguration.data),
                             streamCoordination);

  const unsigned ngroups = cells.getNumberOfCells();
  //this should make the loading of the shared memory a little more efficient
//TODO: measure if using the number of atoms in the inner cell is faster
  const unsigned threads = std::min(maxNumThreads,maxExpected);
  //const unsigned threads = std::min(maxNumThreads,biggestCell);
  plumed_assert(biggestCell <= maxNumThreads) << "the number of atoms in a single cell exceds the maximum number of threads avaiable";
  plumed_assert(memSize <= maxDynamiSharedMemory) << "the shared memory asked exceed the limit for this GPU:"
      << memSize << " > " << maxDynamiSharedMemory;
  runSelf(
    ngroups,
    threads,
    nat,
    memSize,
    streamCoordination,
    biggestCell,
    cellConfiguration.deviceMap_nat_InCells(deviceData),
    cellConfiguration.deviceMap_atomsInCells(deviceData),
    maxExpected,
    cellConfiguration.deviceMap_nat_otherAtoms(deviceData),
    cellConfiguration.deviceMap_otherAtoms(deviceData)
  );
  return nat;
}

template <bool usePBC, typename calculateFloat>
__global__ void
getDualCoord (const unsigned nat,
              const unsigned max_inCell,
              const unsigned max_neigh,
              const PLMD::GPU::rationalSwitchParameters<calculateFloat> switchingParameters,
              const PLMD::GPU::ortoPBCs<calculateFloat> myPBC,
              const unsigned* n_inCell,
              const unsigned* cellIndexes,
              const unsigned* n_neigh,
              const unsigned* neighIndexes,
              const calculateFloat* coordinates,
              const unsigned* trueIndexes,
              calculateFloat* ncoordOut,
              calculateFloat* devOut,
              calculateFloat* virialOut) {
  unsigned cellId=blockIdx.x;
  CUDAHELPERS::sharedArena arena;
  auto sPos = arena.get_shared_memory<calculateFloat>(3 * n_neigh[cellId]);
  auto realIndexes = arena.get_shared_memory<unsigned>(n_neigh[cellId]);
  // // loading shared memory
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
  calculateFloat mycoord = 0.0;
  // the previous version used static array for myVirial and d
  // using explicit variables guarantees that this data will be stored in
  // registers
  calculateFloat myVirial_0 = 0.0;
  calculateFloat myVirial_1 = 0.0;
  calculateFloat myVirial_2 = 0.0;
  calculateFloat myVirial_3 = 0.0;
  calculateFloat myVirial_4 = 0.0;
  calculateFloat myVirial_5 = 0.0;
  calculateFloat myVirial_6 = 0.0;
  calculateFloat myVirial_7 = 0.0;
  calculateFloat myVirial_8 = 0.0;
  // local calculation aid
  const calculateFloat x = coordinates[X (i)];
  const calculateFloat y = coordinates[Y (i)];
  const calculateFloat z = coordinates[Z (i)];
  calculateFloat d_0, d_1, d_2;
  calculateFloat t;
  calculateFloat dfunc;
  for (unsigned j = 0; j < n_neigh[cellId]; ++j) {
    const unsigned jB=neighIndexes[j];
    // Safeguard
    if (idx == realIndexes[j] ) {
      continue;
    }

    d_0 = calculatePBC<usePBC> (sPos[X(j)] - x, myPBC.X);
    d_1 = calculatePBC<usePBC> (sPos[Y(j)] - y, myPBC.Y);
    d_2 = calculatePBC<usePBC> (sPos[Z(j)] - z, myPBC.Z);

    dfunc = 0.;
    mycoord += calculateSqr (
                 d_0 * d_0 + d_1 * d_1 + d_2 * d_2,
                 switchingParameters,
                 dfunc);

    t = -dfunc * d_0;
    mydevX += t;

    myVirial_0 += t * d_0;
    myVirial_1 += t * d_1;
    myVirial_2 += t * d_2;

    t = -dfunc * d_1;
    mydevY += t;

    myVirial_3 += t * d_0;
    myVirial_4 += t * d_1;
    myVirial_5 += t * d_2;

    t = -dfunc * d_2;
    mydevZ += t;

    myVirial_6 += t * d_0;
    myVirial_7 += t * d_1;
    myVirial_8 += t * d_2;
  }
  // working in global memory ONLY at the end
  devOut[X (i)] = mydevX;
  devOut[Y (i)] = mydevY;
  devOut[Z (i)] = mydevZ;
  ncoordOut[i] = mycoord;
  virialOut[nat * 0 + i] = myVirial_0;
  virialOut[nat * 1 + i] = myVirial_1;
  virialOut[nat * 2 + i] = myVirial_2;
  virialOut[nat * 3 + i] = myVirial_3;
  virialOut[nat * 4 + i] = myVirial_4;
  virialOut[nat * 5 + i] = myVirial_5;
  virialOut[nat * 6 + i] = myVirial_6;
  virialOut[nat * 7 + i] = myVirial_7;
  virialOut[nat * 8 + i] = myVirial_8;
}

template <bool usePBC, typename calculateFloat>
__global__ void
getDualDev (
  const unsigned max_inCell,
  const unsigned max_neigh,
  const PLMD::GPU::rationalSwitchParameters<calculateFloat> switchingParameters,
  const PLMD::GPU::ortoPBCs<calculateFloat> myPBC,
  const unsigned* n_inCell,
  const unsigned* cellIndexes,
  const unsigned* n_neigh,
  const unsigned* neighIndexes,
  const calculateFloat* coordinates,
  const unsigned* trueIndexes,
  calculateFloat* devOut
) {
  unsigned cellId=blockIdx.x;
  CUDAHELPERS::sharedArena arena;
  auto sPos = arena.get_shared_memory<calculateFloat>(3 * n_neigh[cellId]);
  auto realIndexes = arena.get_shared_memory<unsigned>(n_neigh[cellId]);
  // // loading shared memory
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
  calculateFloat d_0, d_1, d_2;
  calculateFloat t;
  calculateFloat dfunc;

  for (unsigned j = 0; j < n_neigh[cellId]; ++j) {
    const unsigned jB=neighIndexes[j];
    // Safeguard
    if (idx == realIndexes[j] ) {
      continue;
    }

    d_0 = calculatePBC<usePBC> (sPos[X(j)] - x, myPBC.X);
    d_1 = calculatePBC<usePBC> (sPos[Y(j)] - y, myPBC.Y);
    d_2 = calculatePBC<usePBC> (sPos[Z(j)] - z, myPBC.Z);

    dfunc = 0.;
    t = calculateSqr (
          d_0 * d_0 + d_1 * d_1 + d_2 * d_2,
          switchingParameters,
          dfunc);

    mydevX -= dfunc * d_0;
    mydevY -= dfunc * d_1;
    mydevZ -= dfunc * d_2;
  }
  // working in global memory ONLY at the end
  devOut[X (i)] = mydevX;
  devOut[Y (i)] = mydevY;
  devOut[Z (i)] = mydevZ;
}

template <typename calculateFloat>
void CudaCoordination<calculateFloat>::runDualCoord(
  unsigned const ngroups,
  unsigned const threads,
  unsigned const nat,
  unsigned const memSize,
  cudaStream_t mystream,
  unsigned const maxAtomsInCells,
  unsigned * const nat_InCells,
  unsigned * const atomsInCellsIdxs,
  unsigned const maxOtherAtoms,
  unsigned * const nat_OtherAtoms,
  unsigned * const otherAtomsIdxs
) {
  if(pbc) {
    getDualCoord<true>
    <<<ngroups,
    threads,
    memSize,
    mystream
    >>> (
      nat,
      maxAtomsInCells,
      maxOtherAtoms,
      switchingParameters,
      myPBC,
      nat_InCells,
      atomsInCellsIdxs,
      nat_OtherAtoms,
      otherAtomsIdxs,
      thrust::raw_pointer_cast (cudaPositions.data()),
      thrust::raw_pointer_cast (cudaTrueIndexes.data()),
      thrust::raw_pointer_cast (cudaCoordination.data()),
      thrust::raw_pointer_cast (cudaDerivatives.data()),
      thrust::raw_pointer_cast (cudaVirial.data())
    );
  } else {
    getDualCoord<false>
    <<<ngroups,
    threads,
    memSize,
    mystream
    >>> (
      nat,
      maxAtomsInCells,
      maxOtherAtoms,
      switchingParameters,
      myPBC,
      nat_InCells,
      atomsInCellsIdxs,
      nat_OtherAtoms,
      otherAtomsIdxs,
      thrust::raw_pointer_cast (cudaPositions.data()),
      thrust::raw_pointer_cast (cudaTrueIndexes.data()),
      thrust::raw_pointer_cast (cudaCoordination.data()),
      thrust::raw_pointer_cast (cudaDerivatives.data()),
      thrust::raw_pointer_cast (cudaVirial.data())
    );
  }
}

template <typename calculateFloat>
void CudaCoordination<calculateFloat>::runDualDev(
  unsigned const ngroups,
  unsigned const threads,
  unsigned const /*nat*/,
  unsigned const memSize,
  cudaStream_t mystream,
  unsigned const maxAtomsInCells,
  unsigned * const nat_InCells,
  unsigned * const atomsInCellsIdxs,
  unsigned const maxOtherAtoms,
  unsigned * const nat_OtherAtoms,
  unsigned * const otherAtomsIdxs
) {
  if(pbc) {
    getDualDev<true>
    <<<ngroups,
    threads,
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
      thrust::raw_pointer_cast (cudaPositions.data()),
      thrust::raw_pointer_cast (cudaTrueIndexes.data()),
      thrust::raw_pointer_cast (cudaDerivatives.data())
    );
  } else {
    getDualDev<false>
    <<<ngroups,
    threads,
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
      thrust::raw_pointer_cast (cudaPositions.data()),
      thrust::raw_pointer_cast (cudaTrueIndexes.data()),
      thrust::raw_pointer_cast (cudaDerivatives.data())
    );
  }
}



template <typename calculateFloat>
size_t CudaCoordination<calculateFloat>::doDual() {
  const size_t nat = cudaPositions.size() / 3;

  std::vector<unsigned> cells_required(27);
  const auto maxExpected_dev=cs.listA.getMaximimumCombination(27);
  const auto biggestCell_coord = (*std::max_element(cs.listA.lcell_tots.begin(),cs.listA.lcell_tots.end()));
  //shared memory needed in the derivative loop
  const auto memSize_dev = maxExpected_dev * (3*sizeof(calculateFloat) + sizeof(unsigned));
  const auto maxExpected_coord=cs.listB.getMaximimumCombination(27);
  const auto biggestCell_dev = (*std::max_element(cs.listB.lcell_tots.begin(),cs.listB.lcell_tots.end()));
  //shared memory needed in the coord loop
  const auto memSize_coord = maxExpected_coord * (3*sizeof(calculateFloat) + sizeof(unsigned));
  /**********************allocating the memory on the GPU**********************/
  cudaCoordination.resize (atomsInA);
  cudaVirial.resize (atomsInA * 9);
  /**************************starting the calculations*************************/
//allocating all the memory
  cellSetup cellConfiguration_coord(cells.getNumberOfCells(),biggestCell_coord,atomsInA,maxExpected_coord,nat);

  cellSetup cellConfiguration_dev(cells.getNumberOfCells(),biggestCell_dev,nat,maxExpected_dev,atomsInA);
  thrust::device_vector<unsigned> deviceData_coord;
  thrust::device_vector<unsigned> deviceData_dev;
  plumed_assert(biggestCell_dev <= maxNumThreads) << "the number of atoms in a single cell exceds the maximum number of threads avaiable";
  plumed_assert(biggestCell_coord <= maxNumThreads) << "the number of atoms in a single cell exceds the maximum number of threads avaiable";
  plumed_assert(memSize_coord <= maxDynamiSharedMemory) << "the shared memory asked exceed the limit for this GPU:"
      << memSize_coord << " > " << maxDynamiSharedMemory;
  plumed_assert(memSize_dev <= maxDynamiSharedMemory) << "the shared memory asked exceed the limit for this GPU:"
      << memSize_dev << " > " << maxDynamiSharedMemory;
  const unsigned ngroups = cells.getNumberOfCells();
  const unsigned threads_coord = std::min(maxNumThreads,maxExpected_coord);
  const unsigned threads_dev = std::min(maxNumThreads,maxExpected_dev);
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
                 maxExpected_dev,
                 biggestCell_dev,
                 pbc,
                 cellConfiguration_dev.get_atomsInCells(),
                 cellConfiguration_dev.get_otherAtoms(),
                 cellConfiguration_dev.get_nat_InCells(),
                 cellConfiguration_dev.get_nat_otherAtoms());
  CUDAHELPERS::plmdDataToGPU(deviceData_coord,
                             make_view(cellConfiguration_coord.data),
                             streamCoordination);

  CUDAHELPERS::plmdDataToGPU(deviceData_dev,
                             make_view(cellConfiguration_dev.data),
                             streamDerivatives);

  runDualCoord(
    ngroups,
    threads_coord,
    atomsInA,
    memSize_coord,
    streamCoordination,
    biggestCell_coord,
    cellConfiguration_coord.deviceMap_nat_InCells(deviceData_coord),
    cellConfiguration_coord.deviceMap_atomsInCells(deviceData_coord),
    maxExpected_coord,
    cellConfiguration_coord.deviceMap_nat_otherAtoms(deviceData_coord),
    cellConfiguration_coord.deviceMap_otherAtoms(deviceData_coord)
  );

  //this should make the loading of the shared memory a little more efficient
  runDualDev(
    ngroups,
    threads_dev,
    atomsInB,//not used actually
    memSize_dev,
    streamDerivatives,
    biggestCell_dev,
    cellConfiguration_dev.deviceMap_nat_InCells(deviceData_dev),
    cellConfiguration_dev.deviceMap_atomsInCells(deviceData_dev),
    maxExpected_dev,
    cellConfiguration_dev.deviceMap_nat_otherAtoms(deviceData_dev),
    cellConfiguration_dev.deviceMap_otherAtoms(deviceData_dev)
  );
  return atomsInA;
}

constexpr unsigned couplesPerThreads = 8;
template <bool usePBC, typename calculateFloat>
__global__ void getCoordPair (
  const unsigned couples,
  const unsigned totThreads,
  const PLMD::GPU::rationalSwitchParameters<calculateFloat>
  switchingParameters,
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
  // we try working with less global memory possible, so we set up a bunch of
  // temporary variables
  // const unsigned idx = ;
  // const unsigned jdx = trueIndexes[j];

  // the previous version used static array for myVirial and d
  // using explicit variables guarantees that this data will be stored in
  // registers
  calculateFloat myVirial_0 = 0.0;
  calculateFloat myVirial_1 = 0.0;
  calculateFloat myVirial_2 = 0.0;
  calculateFloat myVirial_3 = 0.0;
  calculateFloat myVirial_4 = 0.0;
  calculateFloat myVirial_5 = 0.0;
  calculateFloat myVirial_6 = 0.0;
  calculateFloat myVirial_7 = 0.0;
  calculateFloat myVirial_8 = 0.0;
  // local calculation aid

  calculateFloat d_0, d_1, d_2;
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
    d_0 = calculatePBC<usePBC> (coordinatesJ[X (i)] - coordinatesI[X (i)],
                                myPBC.X);
    d_1 = calculatePBC<usePBC> (coordinatesJ[Y (i)] - coordinatesI[Y (i)],
                                myPBC.Y);
    d_2 = calculatePBC<usePBC> (coordinatesJ[Z (i)] - coordinatesI[Z (i)],
                                myPBC.Z);
    dfunc = 0.;
    if (trueIndexes[i] != trueIndexes[i + couples]) {
      myCoord += calculateSqr (
                   d_0 * d_0 + d_1 * d_1 + d_2 * d_2, switchingParameters, dfunc);
    } else {
      d_0 = d_1 = d_2 = 0.0;
    }

    t = -dfunc * d_0;
    myVirial_0 += t * d_0;
    myVirial_1 += t * d_1;
    myVirial_2 += t * d_2;
    devOutI[X (i)] = t;
    devOutJ[X (i)] = -t;

    t = -dfunc * d_1;
    myVirial_3 += t * d_0;
    myVirial_4 += t * d_1;
    myVirial_5 += t * d_2;
    devOutI[Y (i)] = t;
    devOutJ[Y (i)] = -t;

    t = -dfunc * d_2;
    myVirial_6 += t * d_0;
    myVirial_7 += t * d_1;
    myVirial_8 += t * d_2;
    devOutI[Z (i)] = t;
    devOutJ[Z (i)] = -t;
  }
  ncoordOut[tid] = myCoord;
  virialOut[totThreads * 0 + tid] = myVirial_0;
  virialOut[totThreads * 1 + tid] = myVirial_1;
  virialOut[totThreads * 2 + tid] = myVirial_2;
  virialOut[totThreads * 3 + tid] = myVirial_3;
  virialOut[totThreads * 4 + tid] = myVirial_4;
  virialOut[totThreads * 5 + tid] = myVirial_5;
  virialOut[totThreads * 6 + tid] = myVirial_6;
  virialOut[totThreads * 7 + tid] = myVirial_7;
  virialOut[totThreads * 8 + tid] = myVirial_8;
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
  if (pbc) {
    getCoordPair<true><<<ngroups, maxNumThreads, 0, streamDerivatives>>> (
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
  } else {
    getCoordPair<false><<<ngroups, maxNumThreads, 0, streamDerivatives>>> (
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
  }
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
    NL.cutoff = dmax*1.01;
    NL.stride=1;
  }

  checkRead();
  if (mpiActive) {

    cudaStreamCreate (&streamDerivatives);
    cudaStreamCreate (&streamVirial);
    cudaStreamCreate (&streamCoordination);

    setUpPermanentGPUMemory();

    maxReductionNumThreads = min (1024, maxNumThreads);

    cudaFuncAttributes attr;
    // the kernels are heavy on registers, this adjust the maximum number of
    // threads accordingly
    switch (mode) {
    case calculationMode::self:
      if (pbc) {
        cudaFuncGetAttributes (&attr, &getSelfCoord<true, calculateFloat>);
      } else {
        cudaFuncGetAttributes (&attr, &getSelfCoord<false, calculateFloat>);
      }
      break;
    case calculationMode::dual:
      if (pbc) {
        cudaFuncGetAttributes (&attr, &getDualDev<true, calculateFloat>);
        maxNumThreads = min (attr.maxThreadsPerBlock, maxNumThreads);
        cudaFuncGetAttributes (&attr, &getDualDev<true, calculateFloat>);
      } else {
        cudaFuncGetAttributes (&attr, &getDualDev<false, calculateFloat>);
        maxNumThreads = min (attr.maxThreadsPerBlock, maxNumThreads);
        cudaFuncGetAttributes (&attr, &getDualCoord<false, calculateFloat>);
      }
      break;
    case calculationMode::pair:
      if (pbc) {
        cudaFuncGetAttributes (&attr, &getCoordPair<true, calculateFloat>);
      } else {
        cudaFuncGetAttributes (&attr, &getCoordPair<false, calculateFloat>);
      }
      break;
    case calculationMode::none:
      // throw"this should not have been happened"
      break;
    }
    maxNumThreads = min (attr.maxThreadsPerBlock, maxNumThreads);
    maxDynamiSharedMemory=attr.maxDynamicSharedSizeBytes;
  }
  comm.Bcast (maxNumThreads, 0);
  comm.Bcast (maxReductionNumThreads, 0);
  comm.Bcast (maxDynamiSharedMemory, 0);

  log << "  contacts are counted with cutoff (dmax)="
      << sqrt (switchingParameters.dmaxSQ)
      << ", with a rational switch with parameters: "
      "d0="<< switchingParameters.d0
      <<", r0="<< 1.0 / ((switchingParameters.calcSquared)?(sqrt (switchingParameters.invr0_2)):switchingParameters.invr0_2)
      << ", N=" << switchingParameters.nn * ((switchingParameters.calcSquared)?2:1)
      << ", M=" << switchingParameters.mm * ((switchingParameters.calcSquared)?2:1)
      << ".\n";
  log << "GPU info:\n"
      << "\t max threads per coordination" << maxNumThreads << "\n"
      <<"\t max shared memory" << maxDynamiSharedMemory << "\n"
      << "\t max threads per reduction" << maxReductionNumThreads << "\n";

  // cudaFuncGetAttributes (&attr, &getSelfCoord<true, calculateFloat>);
  // std::cout<< "Attributes for pbc: \n";

  // std::cout << "\tmaxDynamicSharedSizeBytes: " <<
  // attr.maxDynamicSharedSizeBytes <<"\n" ; std::cout << "\tmaxThreadsPerBlock:
  // " << attr.maxThreadsPerBlock <<"\n" ; std::cout << "\tnumRegs: " <<
  // attr.numRegs <<"\n" ; std::cout << "\tsharedSizeBytes: " <<
  // attr.sharedSizeBytes <<"\n" ;

  // cudaFuncGetAttributes ( &attr, &getSelfCoord<false,calculateFloat> );
  // std::cout<< "Attributes for no pbc: \n";

  // std::cout << "\tmaxDynamicSharedSizeBytes: " <<
  // attr.maxDynamicSharedSizeBytes <<"\n" ; std::cout << "\tmaxThreadsPerBlock:
  // " << attr.maxThreadsPerBlock <<"\n" ; std::cout << "\tnumRegs: " <<
  // attr.numRegs <<"\n" ; std::cout << "\tsharedSizeBytes: " <<
  // attr.sharedSizeBytes <<"\n" ;

  // cudaFuncGetAttributes ( &attr,
  // &CUDAHELPERS::reduction1DKernel<calculateFloat, 128, 4> ); std::cout<<
  // "Attributes for reduction kernel (128): \n";

  // std::cout << "\tmaxDynamicSharedSizeBytes: " <<
  // attr.maxDynamicSharedSizeBytes <<"\n" ; std::cout << "\tmaxThreadsPerBlock:
  // " << attr.maxThreadsPerBlock <<"\n" ; std::cout << "\tnumRegs: " <<
  // attr.numRegs <<"\n" ; std::cout << "\tsharedSizeBytes: " <<
  // attr.sharedSizeBytes <<"\n" ;
}

} // namespace colvar
} // namespace PLMD
