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
on the GPU: COORDINATIONNUMBER_SUMSACC, Q1_SUMSACC,...,Q6_SUMSACC, Q1_DENOM_SUMSACC,...

They are the actions of the core (keywords, outputs, cells) with the kernels of the core run on
the GPU, one GPU thread per row: only the positions and the cells go to the GPU, and only the
sums and the forces on the atoms come back. The link cells are built on the CPU (LinkCells) and
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
  void computeSums() override;
  void computeForces( std::vector<double>& outforces ) override;
public:
  static void registerKeywords( Keywords& keys );
  explicit NeighborSumsACC( const ActionOptions& );
  ~NeighborSumsACC();
};

typedef NeighborSumsACC<symfunc::ContactPairs> ContactSumsACC;
PLUMED_REGISTER_ACTION(ContactSumsACC,"COORDINATIONNUMBER_SUMSACC")
PLUMED_REGISTER_ACTION(ContactSumsACC,"Q1_DENOM_SUMSACC")
PLUMED_REGISTER_ACTION(ContactSumsACC,"Q3_DENOM_SUMSACC")
PLUMED_REGISTER_ACTION(ContactSumsACC,"Q4_DENOM_SUMSACC")
PLUMED_REGISTER_ACTION(ContactSumsACC,"Q6_DENOM_SUMSACC")
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
  // cells per cutoff searched on the GPU (temporary, to compare 1 and 2 on the GPUs)
  if( const char* sub=std::getenv("PLUMED_NEIGHBORSUMS_SUBDIVISION") ) {
    this->cells.setSubdivision( std::atoi(sub) );
    this->log.printf("  %s cells per cutoff (PLUMED_NEIGHBORSUMS_SUBDIVISION)\n",sub);
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
void NeighborSumsACC<PF>::computeSums() {
  ACC::DeviceScope scope( device.resolve(*this) );
  const unsigned ncomp = this->pairs.ncomponents;
  NeighborSumsWorkspace& w=*ws;
  // the positions and the cells are also used by the backward pass of this step
  const symfunc::NeighborSumsGeometry g = toDevice();
  w.sums.resize( this->mycount*ncomp );
  double* out = w.sums.devicePtr();
  const std::size_t first=this->myfirst, count=this->mycount;
  // a local copy on the device for the kernel (as the task manager does with the action data)
  PF data = this->pairs;
  auto data_acc = OpenACC::fromToDataHelper( data );
#pragma acc parallel loop gang vector async(accQueue) present(data) \
                    firstprivate(g, ncomp, first, count) deviceptr(out) default(none)
  for(std::size_t rr=0; rr<count; ++rr) {
    symfunc::neighborSumsRow( data, g, first+rr, out+rr*ncomp );
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
  NeighborSumsAtomicColumns cols{forces,this->colstart};
  double v0=0, v1=0, v2=0, v3=0, v4=0, v5=0, v6=0, v7=0, v8=0;
  PF data = this->pairs;
  auto data_acc = OpenACC::fromToDataHelper( data );
#pragma acc parallel loop gang vector async(accQueue) present(data) \
                    firstprivate(g, cols, first, count, rowstart) deviceptr(G, forces) \
                    reduction(+:v0,v1,v2,v3,v4,v5,v6,v7,v8) default(none)
  for(std::size_t rr=0; rr<count; ++rr) {
    const std::size_t r = first+rr;
    double vir[9]= {0,0,0,0,0,0,0,0,0};
    const Vector fr = symfunc::neighborForcesRow( data, g, G, r, cols, vir );
    // each row is a different atom, and the columns are other atoms
    const std::size_t i = 3*(rowstart+r);
    forces[i] += fr[0];
    forces[i+1] += fr[1];
    forces[i+2] += fr[2];
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
  w.forces.copyFromDeviceAsync( outforces.data(), accQueue );
  acc_wait(accQueue);
  const double vir[9]= {v0,v1,v2,v3,v4,v5,v6,v7,v8};
  std::copy( vir, vir+9, outforces.begin()+3*natoms );
}

} // namespace openacc
} // namespace PLMD
