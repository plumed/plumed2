/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2026 The plumed team
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
/*
The sums over the neighbours of COORDINATIONNUMBER and Q1,...,Q6 (src/symfunc/NeighborSums.h)
on the GPU: COORDINATIONNUMBER_SUMSACC, Q1_SUMSACC,...,Q6_SUMSACC (whose coordination numbers,
the denominators of Q1,...,Q6, are computed in the same pass and given by Qn_DENOM_SUMS).

They are the actions of the core (keywords, outputs, cells) with the kernels of the core run on
the GPU: only the positions and the cells go to the GPU, and only the sums and the forces on the
atoms come back. The search around each row is split between several GPU threads (parts of the
rows of cells around it, or blocks of columns with few rows), so that the GPU has enough threads
also with few rows; the results of the parts of a row are summed in order by a second kernel. The link cells are built on the CPU (LinkCells) and
searched on the GPU (LinkCellsView, forEachNeighbor). With GROUP the force on each atom is
gathered from its own neighbours (no atomics); with GROUPA/GROUPB the forces on the columns are
added with atomics.

The COORDINATIONNUMBER and Q1,...,Q6 shortcuts use them with USEGPU when the plugin is loaded.
With MPI the first process computes everything and sends the results to the others; with
MULTIGPU (and without SERIAL) the rows are split in blocks between the processes, each using its
GPU, and the results are summed.
*/
#include "plumed/symfunc/NeighborSums.h"
#include "plumed/core/ActionRegister.h"
#include "plumed/tools/OpenACC.h"

#include "ACCDevice.h"
#include "ACCParallelTaskManager.h"

#include <cstdlib>
#include <memory>

namespace PLMD {
namespace openacc {

/// The device copies of the input and of the results, kept between steps
struct NeighborSumsWorkspace {
  OpenACC::memoryManager<double> positions;
  OpenACC::memoryManager<unsigned> cellStarts;
  OpenACC::memoryManager<unsigned> cellTots;
  OpenACC::memoryManager<unsigned> cellLists;
  OpenACC::memoryManager<double> cellPositions;
  AccCachedBuffer<Pbc> pbc;
  OpenACC::memoryManager<double> sums;
  OpenACC::memoryManager<double> forcesOnSums;
  OpenACC::memoryManager<double> forces;
  /// the forces on the atoms of the rows from each part of the search
  OpenACC::memoryManager<double> rowForces;
};

/// The forces on the columns (without GROUP), added with atomics
struct NeighborSumsAtomicColumns {
  double* forces;
  std::size_t colstart;
#pragma acc routine seq
  void add( const unsigned c, const double g0, const double g1, const double g2 ) {
    const std::size_t j=3*(colstart+c);
#pragma acc atomic update
    forces[j] += g0;
#pragma acc atomic update
    forces[j+1] += g1;
#pragma acc atomic update
    forces[j+2] += g2;
  }
};

template <class PF>
class NeighborSumsACC : public symfunc::NeighborSums<PF> {
  using Base = symfunc::NeighborSums<PF>;
  bool multigpu=false;
  ACC::Device device;
  std::unique_ptr<NeighborSumsWorkspace> ws;
  bool splitRows() const override {
    return multigpu && !this->serial && this->comm.Get_size()>1;
  }
  bool firstProcessOnly() const override {
    return true;
  }
  /// the geometry of this step on the device
  symfunc::NeighborSumsGeometry toDevice();
  /// enough GPU threads: each row is searched by several threads (see neighborSumsVisit)
  unsigned numberOfSlices() const override;
  void computeSums() override;
  void computeForces( std::vector<double>& outforces ) override;
public:
  static void registerKeywords( Keywords& keys );
  explicit NeighborSumsACC( const ActionOptions& );
  ~NeighborSumsACC();
};

typedef NeighborSumsACC<symfunc::ContactPairs> ContactSumsACC;
PLUMED_REGISTER_ACTION(ContactSumsACC,"COORDINATIONNUMBER_SUMSACC")
typedef NeighborSumsACC<symfunc::SphericalHarmonicPairs> SphericalHarmonicSumsACC;
PLUMED_REGISTER_ACTION(SphericalHarmonicSumsACC,"Q1_SUMSACC")
PLUMED_REGISTER_ACTION(SphericalHarmonicSumsACC,"Q3_SUMSACC")
PLUMED_REGISTER_ACTION(SphericalHarmonicSumsACC,"Q4_SUMSACC")
PLUMED_REGISTER_ACTION(SphericalHarmonicSumsACC,"Q6_SUMSACC")

template <class PF>
void NeighborSumsACC<PF>::registerKeywords( Keywords& keys ) {
  Base::registerKeywords( keys );
  keys.addFlag("MULTIGPU",false,"with MPI, share the atoms between the processes, each using its own GPU, and sum the results (by default only the first process computes, on its GPU)");
  keys.addFlag("USEGPU",false,"ignored: this action always runs on the GPU");
}

template <class PF>
NeighborSumsACC<PF>::NeighborSumsACC( const ActionOptions& ao ):
  Action(ao),
  Base(ao),
  ws(std::make_unique<NeighborSumsWorkspace>()) {
  this->parseFlag("MULTIGPU",multigpu);
  bool usegpu;
  this->parseFlag("USEGPU",usegpu);
  if( multigpu && !this->serial ) {
    this->log.printf("  with MPI the atoms are shared between the processes, each using its GPU\n");
  }
  this->checkRead();
}

template <class PF>
NeighborSumsACC<PF>::~NeighborSumsACC() {
  // the device memory is freed on its GPU
  ACC::DeviceScope scope(device.id());
  ws.reset();
}

template <class PF>
symfunc::NeighborSumsGeometry NeighborSumsACC<PF>::toDevice() {
  NeighborSumsWorkspace& w=*ws;
  const auto & pos = this->getPositions();
  w.positions.resize( 3*pos.size() );
  w.positions.copyToDeviceAsync( &pos[0][0], accQueue );
  symfunc::NeighborSumsGeometry g = this->geometry;
  g.positions = w.positions.devicePtr();
  g.pbc = w.pbc.update( &this->getPbc(), 1 );
  if( g.usecells ) {
    const LinkCellsView& v = this->geometry.cells;
    w.cellStarts.resize( v.ncellstot );
    w.cellStarts.copyToDeviceAsync( v.starts, accQueue );
    w.cellTots.resize( v.ncellstot );
    w.cellTots.copyToDeviceAsync( v.tots, accQueue );
    w.cellLists.resize( v.natoms );
    w.cellLists.copyToDeviceAsync( v.lists, accQueue );
    w.cellPositions.resize( 3*v.natoms );
    w.cellPositions.copyToDeviceAsync( &v.positions[0][0], accQueue );
    g.cells.starts = w.cellStarts.devicePtr();
    g.cells.tots = w.cellTots.devicePtr();
    g.cells.lists = w.cellLists.devicePtr();
    g.cells.positions = reinterpret_cast<const Vector*>( w.cellPositions.devicePtr() );
  }
  return g;
}

template <class PF>
unsigned NeighborSumsACC<PF>::numberOfSlices() const {
  // TEMPORARY, for the measurements on the GPU: PLUMED_NEIGHBORSUMS_SLICES (0: automatic)
  static const long env = std::getenv("PLUMED_NEIGHBORSUMS_SLICES") ? std::atol(std::getenv("PLUMED_NEIGHBORSUMS_SLICES")) : 0;
  if( env>0 ) {
    return env;
  }
  // about this many GPU threads (tasks of a row and a part of its search)
  constexpr std::size_t tasks=std::size_t(1)<<17;
  const std::size_t rows=std::max<std::size_t>( 1, this->mycount );
  std::size_t S=std::max<std::size_t>( 1, tasks/rows );
  if( this->geometry.usecells ) {
    // the rows of cells around a cell
    constexpr std::size_t side=2*LinkCellsView::maxSubdivision+1;
    S=std::min( S, side*side );
  } else {
    // at least 32 columns per part
    S=std::min( S, std::max<std::size_t>( 1, this->ncols/32 ) );
  }
  return S;
}

template <class PF>
void NeighborSumsACC<PF>::computeSums() {
  ACC::DeviceScope scope( device.resolve(*this) );
  const unsigned ncomp = this->pairs.ncomponents;
  NeighborSumsWorkspace& w=*ws;
  // the positions and the cells are also used by the backward pass of this step
  const symfunc::NeighborSumsGeometry g = toDevice();
  const std::size_t first=this->myfirst, count=this->mycount;
  // one GPU thread for each row and part of its search (row-major: the threads of a row are
  // next to each other), then the sums of the parts of each row in order
  const unsigned S=this->nslices, s0=this->slicefirst, ns=this->sliceend-this->slicefirst;
  const std::size_t ntasks=count*ns;
  w.sums.resize( ntasks*ncomp + count*ncomp );
  double* parts = w.sums.devicePtr();
  double* out = parts + ntasks*ncomp;
  // a local copy on the device for the kernel (as the task manager does with the action data)
  PF data = this->pairs;
  auto data_acc = OpenACC::fromToDataHelper( data );
#pragma acc parallel loop gang vector async(accQueue) present(data) \
                    firstprivate(g, ncomp, first, ntasks, ns, s0, S) deviceptr(parts) default(none)
  for(std::size_t t=0; t<ntasks; ++t) {
    symfunc::neighborSumsRow( data, g, first+t/ns, parts+t*ncomp, s0+t%ns, S );
  }
#pragma acc parallel loop gang vector async(accQueue) \
                    firstprivate(ncomp, count, ns) deviceptr(parts, out) default(none)
  for(std::size_t e=0; e<count*ncomp; ++e) {
    const std::size_t rr=e/ncomp, k=e%ncomp;
    double v=0.0;
    for(unsigned s=0; s<ns; ++s) {
      v += parts[(rr*ns+s)*ncomp+k];
    }
    out[e]=v;
  }
  acc_memcpy_from_device_async( this->sums.data()+first*ncomp, out, count*ncomp*sizeof(double), accQueue );
  acc_wait(accQueue);
}

template <class PF>
void NeighborSumsACC<PF>::computeForces( std::vector<double>& outforces ) {
  ACC::DeviceScope scope( device.resolve(*this) );
  const unsigned ncomp = this->pairs.ncomponents;
  const std::size_t natoms = this->getNumberOfAtoms();
  NeighborSumsWorkspace& w=*ws;
  // the positions and the cells are those of the forward pass of this step
  symfunc::NeighborSumsGeometry g = this->geometry;
  g.positions = w.positions.devicePtr();
  g.pbc = w.pbc.update( &this->getPbc(), 1 );
  if( g.usecells ) {
    g.cells.starts = w.cellStarts.devicePtr();
    g.cells.tots = w.cellTots.devicePtr();
    g.cells.lists = w.cellLists.devicePtr();
    g.cells.positions = reinterpret_cast<const Vector*>( w.cellPositions.devicePtr() );
  }
  // the forces on the sums of all the rows (with GROUP those of the neighbours are needed too)
  w.forcesOnSums.resize( this->nrows*ncomp );
  w.forcesOnSums.copyToDeviceAsync( this->forcesOnSums.data(), accQueue );
  w.forces.resize( 3*natoms );
  accZeroAsync( w.forces.devicePtr(), 3*natoms );
  const double* G = w.forcesOnSums.devicePtr();
  double* forces = w.forces.devicePtr();
  const std::size_t first=this->myfirst, count=this->mycount, rowstart=this->rowstart;
  const unsigned S=this->nslices, s0=this->slicefirst, ns=this->sliceend-this->slicefirst;
  const std::size_t ntasks=count*ns;
  w.rowForces.resize( 3*ntasks );
  double* rowForces = w.rowForces.devicePtr();
  NeighborSumsAtomicColumns cols{forces,this->colstart};
  double v0=0, v1=0, v2=0, v3=0, v4=0, v5=0, v6=0, v7=0, v8=0;
  PF data = this->pairs;
  auto data_acc = OpenACC::fromToDataHelper( data );
  // one GPU thread for each row and part of its search; the forces on the columns with atomics
#pragma acc parallel loop gang vector async(accQueue) present(data) \
                    firstprivate(g, cols, first, ntasks, ns, s0, S) deviceptr(G, rowForces) \
                    reduction(+:v0,v1,v2,v3,v4,v5,v6,v7,v8) default(none)
  for(std::size_t t=0; t<ntasks; ++t) {
    double vir[9]= {0,0,0,0,0,0,0,0,0};
    const Vector fr = symfunc::neighborForcesRow( data, g, G, first+t/ns, cols, vir, s0+t%ns, S );
    rowForces[3*t] = fr[0];
    rowForces[3*t+1] = fr[1];
    rowForces[3*t+2] = fr[2];
    v0 += vir[0];
    v1 += vir[1];
    v2 += vir[2];
    v3 += vir[3];
    v4 += vir[4];
    v5 += vir[5];
    v6 += vir[6];
    v7 += vir[7];
    v8 += vir[8];
  }
  // the forces on the atoms of the rows, summed over the parts in order (each row is a
  // different atom, and the columns are other atoms)
#pragma acc parallel loop gang vector async(accQueue) \
                    firstprivate(first, count, rowstart, ns) deviceptr(rowForces, forces) default(none)
  for(std::size_t e=0; e<3*count; ++e) {
    const std::size_t rr=e/3, k=e%3;
    double f=0.0;
    for(unsigned s=0; s<ns; ++s) {
      f += rowForces[3*(rr*ns+s)+k];
    }
    forces[3*(rowstart+first+rr)+k] += f;
  }
  w.forces.copyFromDeviceAsync( outforces.data(), accQueue );
  acc_wait(accQueue);
  const double vir[9]= {v0,v1,v2,v3,v4,v5,v6,v7,v8};
  std::copy( vir, vir+9, outforces.begin()+3*natoms );
}

} // namespace openacc
} // namespace PLMD
