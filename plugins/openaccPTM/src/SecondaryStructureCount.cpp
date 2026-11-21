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
#include "plumed/secondarystructure/SecondaryStructureCount.h"
#include "plumed/core/ActionRegister.h"
#include "plumed/tools/OpenACC.h"
#include "ACCDevice.h"
#include "ACCParallelTaskManager.h"

#include <openacc.h>

#include <array>
#include <memory>
#include <vector>

namespace PLMD {
namespace secondarystructure {

namespace {
/// For each segment, the value of the switching function and its contribution to the virial
constexpr std::size_t segmentColumns=10;
}

/// The GPU version of ALPHARMSD_SCALAR, ANTIBETARMSD_SCALAR and PARABETARMSD_SCALAR (only with TYPE=DRMSD), created by
/// ALPHARMSD, ANTIBETARMSD and PARABETARMSD with USEGPU. A step is done in these kernels:
/// - for the sheets, the list of the pairs of strands within STRANDS_CUTOFF and whose switching function is not zero for
///   certain (one thread per pair, see sscount::segmentActive; the radii of the strands are computed first);
/// - one thread per segment of that list (or per segment, for the helices): the segment is built
///   and its distance from the references and the switching function are computed with the same
///   functions as on the CPU (sscount::segmentPositions and sscount::segmentSwitch); the derivatives
///   are added to the atoms with atomics, the value and the virial are written for the next kernel;
/// - the sum of the values and of the virial.
/// The derivatives are copied back to the host and used as in SecondaryStructureCount: there is no backward pass.
/// The pairs, the reference distances and the switching function are copied to the GPU once.
/// The order of the sums is not fixed, so results are reproducible only up to rounding.
class SecondaryStructureCountACC : public SecondaryStructureCount {
  ACC::Device device_;
  SwitchingFunctionAccelerable switchAcc;
  /// the device memory, that must be freed on its GPU
  struct Buffers {
    /// a copy of the reference distances, present on the device
    std::unique_ptr<SecondaryStructureDRMSDInput<double> > drmsd;
    /// the first atom of the two strands of each pair
    OpenACC::memoryManager<unsigned> pairs;
    OpenACC::memoryManager<double> positions;
    AccCachedBuffer<Pbc> pbc;
    AccCachedBuffer<SwitchingFunctionAccelerable> switchingFunction;
    /// the first atom of each strand and the radii of the strands, if the pairs are skipped with skipDistance
    OpenACC::memoryManager<unsigned> strands;
    OpenACC::memoryManager<double> radius;
    /// the pairs that are computed (within STRANDS_CUTOFF and skipDistance), and their number
    OpenACC::memoryManager<unsigned> active;
    OpenACC::memoryManager<unsigned> count;
    /// segmentColumns values per segment: column c of segment e is at c*nactive+e
    OpenACC::memoryManager<double> segments;
    OpenACC::memoryManager<double> derivatives;
    OpenACC::memoryManager<double> sums;
  };
  std::unique_ptr<Buffers> buf_;
  std::vector<double> hder_;
  std::vector<double> hsums_;
  /// The state of the step, between the functions below: the number of pairs of this process,
  /// the number of the ones that are computed (read from the GPU in hcount_ with a list), the device pointers
  std::size_t nmine_=0;
  std::size_t nactive_=0;
  bool uselist_=false;
  unsigned hcount_=0;
  const Pbc* dpbc_=nullptr;
  const SwitchingFunctionAccelerable* dsw_=nullptr;
  /// Queues the kernel of the segments and the sums (the results are downloaded in calculate():
  /// a download into pageable host memory returns only when it is done, so in launch() it
  /// would make launch() wait for the GPU)
  void enqueueSegments();
public:
  static void registerKeywords( Keywords& keys );
  explicit SecondaryStructureCountACC(const ActionOptions&);
  ~SecondaryStructureCountACC() override;
  SecondaryStructureCountACC(const SecondaryStructureCountACC&)=delete;
  SecondaryStructureCountACC& operator=(const SecondaryStructureCountACC&)=delete;
  /// Queues the copy of the positions and the list of the pairs or, without a list, the whole calculation
  void launch() override;
  void cancelLaunch() override;
  void calculate() override;
};

PLUMED_REGISTER_ACTION(SecondaryStructureCountACC,"ALPHARMSD_SCALAR_ACC")
PLUMED_REGISTER_ACTION(SecondaryStructureCountACC,"ANTIBETARMSD_SCALAR_ACC")
PLUMED_REGISTER_ACTION(SecondaryStructureCountACC,"PARABETARMSD_SCALAR_ACC")

void SecondaryStructureCountACC::registerKeywords( Keywords& keys ) {
  SecondaryStructureCount::registerKeywords( keys );
}

SecondaryStructureCountACC::SecondaryStructureCountACC(const ActionOptions&ao):
  Action(ao),
  SecondaryStructureCount(ao),
  buf_(std::make_unique<Buffers>()) {
  if( !usedrmsd ) {
    error("on the GPU only TYPE=DRMSD is available");
  }
  std::string errors;
  switchAcc.set( switchInput, errors );
  if( errors.length()>0 ) {
    error("problem reading SWITCH keyword : " + errors );
  }
  log.printf("  computed on the GPU\n");
  // the work on the GPU is queued as soon as the atoms are available, see Action::launch()
  setLaunchable();
}

SecondaryStructureCountACC::~SecondaryStructureCountACC() {
  // the device memory is freed here, on its GPU
  ACC::DeviceScope scope(device_.id());
  if( buf_ && buf_->drmsd ) {
    buf_->drmsd->removeFromACCDevice();
  }
  buf_.reset();
}

void SecondaryStructureCountACC::launch() {
  const unsigned npairs=pairs.size();
  const unsigned natoms=getNumberOfAtoms();
  const unsigned stride=comm.Get_size();
  const unsigned rank=comm.Get_rank();
  // the pairs rank, rank+stride, ... are computed here, as in SecondaryStructureCount
  const std::size_t nmine = nmine_ = rank<npairs ? (npairs-1-rank)/stride+1 : 0;
  hder_.assign( 3*natoms+9, 0.0 );
  nactive_=0;
  if( nmine==0 ) {
    return;
  }
  ACC::DeviceScope scope(device_.resolve(*this));
  auto& b=*buf_;
  if( !b.drmsd ) {
    b.drmsd=std::make_unique<SecondaryStructureDRMSDInput<double> >( drmsd );
    b.drmsd->toACCDevice();
    std::vector<unsigned> flat( 2*npairs );
    for(unsigned i=0; i<npairs; ++i) {
      flat[2*i]=pairs[i].first;
      flat[2*i+1]=pairs[i].second;
    }
    b.pairs.resize( flat.size() );
    b.pairs.copyToDevice( flat.data() );
    b.active.resize( nmine );
    b.count.resize( 1 );
    b.sums.resize( segmentColumns );
    if( skipDistance>=0 ) {
      b.strands.resize( strands.size() );
      b.strands.copyToDevice( strands.data() );
      b.radius.resize( natoms );
    }
  }
  const unsigned* dpairs=b.pairs.devicePtr();
  b.positions.resize( 3*natoms );
  b.positions.copyToDeviceAsync( &getPositions()[0][0], accQueue );
  const double* dx=b.positions.devicePtr();
  const Pbc* dpbc = dpbc_ = nopbc ? nullptr : b.pbc.update( &getPbc(), 1 );
  dsw_=b.switchingFunction.update( &switchAcc, 1 );
  unsigned* dactive=b.active.devicePtr();

  // the radii of the strands, to skip the pairs whose switching function is zero
  const double* drad=nullptr;
  if( skipDistance>=0 ) {
    const unsigned* dstrands=b.strands.devicePtr();
    double* dradius=b.radius.devicePtr();
    const std::size_t nstrands=strands.size();
#pragma acc parallel loop gang vector async(accQueue) \
                   firstprivate(nstrands) \
                      deviceptr(dx, dpbc, dstrands, dradius) \
                        default(none)
    for(std::size_t i=0; i<nstrands; ++i) {
      dradius[dstrands[i]]=sscount::strandRadius( dx, dstrands[i], dpbc );
    }
    drad=dradius;
  }

  // the list of the pairs that are computed
  uselist_ = strands_cutoff>0 || drad;
  if( uselist_ ) {
    unsigned* dcount=b.count.devicePtr();
    accZeroAsync( dcount, 1 );
    const double cutoff=strands_cutoff;
    const double skip=skipDistance;
#pragma acc parallel loop gang vector async(accQueue) \
                   firstprivate(nmine, rank, stride, cutoff, skip) \
                      deviceptr(dx, dpairs, dpbc, drad, dactive, dcount) \
                        default(none)
    for(std::size_t k=0; k<nmine; ++k) {
      const std::size_t p=rank+k*stride;
      if( sscount::segmentActive( dx, dpairs[2*p], dpairs[2*p+1], dpbc, cutoff, drad, skip ) ) {
        unsigned slot;
#pragma acc atomic capture
        slot=(*dcount)++;
        dactive[slot]=p;
      }
    }
    // the number of pairs in the list is needed to go on: calculate() reads it
  } else {
    nactive_=nmine;
    enqueueSegments();
  }
}

void SecondaryStructureCountACC::enqueueSegments() {
  const std::size_t nactive=nactive_;
  if( nactive==0 ) {
    return;
  }
  auto& b=*buf_;
  const unsigned natoms=getNumberOfAtoms();
  const unsigned stride=comm.Get_size();
  const unsigned rank=comm.Get_rank();
  const bool uselist=uselist_;
  const SecondaryStructureDRMSDInput<double>* dd=b.drmsd.get();
  const unsigned* dpairs=b.pairs.devicePtr();
  const double* dx=b.positions.devicePtr();
  const Pbc* dpbc=dpbc_;
  const SwitchingFunctionAccelerable* dsw=dsw_;
  unsigned* dactive=b.active.devicePtr();
  const std::size_t ncols=segmentColumns;
  b.segments.resize( ncols*nactive );
  double* dseg=b.segments.devicePtr();
  b.derivatives.resize( 3*natoms );
  double* dder=b.derivatives.devicePtr();
  accZeroAsync( dder, 3*natoms );
  const unsigned nstr=nstructures;
  const bool alignStrands=!alpha;
  const unsigned nsegatoms=maxAtomsInSegment;
  const unsigned nstrand=sscount::strandAtoms;
  const unsigned nder=sscount::segmentDerivatives;
  // one thread per segment, with all the inner loops sequential
#pragma acc parallel loop gang vector async(accQueue) present(dd[0:1]) \
                 firstprivate(nactive, uselist, rank, stride, alignStrands, \
                              nstr, nsegatoms, nstrand, nder) \
                    deviceptr(dx, dpairs, dpbc, dsw, dactive, dseg, dder) \
                      default(none)
  for(std::size_t e=0; e<nactive; ++e) {
    // the arrays are declared first in the body of the loop, before its inner loops:
    // otherwise, with default(none), nvc++ asks for a data clause for them
    std::array<Vector,maxAtomsInSegment> pos;
    std::array<double,2*sscount::segmentDerivatives> der;
    unsigned best=0;
    double g=0.0;
    const std::size_t p = uselist ? dactive[e] : rank+e*stride;
    const unsigned s1=dpairs[2*p], s2=dpairs[2*p+1];
    sscount::segmentPositions( dx, s1, s2, dpbc, alignStrands, pos );
    dseg[e]=sscount::segmentSwitch( *dd, nstr, *dsw, pos, der, best, g );
    // with g zero the derivatives were not computed
    const double* d=der.data()+best*nder;
#pragma acc loop seq
    for(unsigned j=0; j<9; ++j) {
      dseg[(j+1)*nactive+e] = (g!=0.0) ? g*d[3*nsegatoms+j] : 0.0;
    }
    if( g!=0.0 ) {
#pragma acc loop seq
      for(unsigned k=0; k<nsegatoms; ++k) {
        const unsigned a = (k<nstrand) ? s1+k : s2+k-nstrand;
#pragma acc loop seq
        for(unsigned j=0; j<3; ++j) {
          const double c=g*d[3*k+j];
#pragma acc atomic update
          dder[3*a+j] += c;
        }
      }
    }
  }

  // the sums of the values and of the virial: one gang per column
  double* dsums=b.sums.devicePtr();
#pragma acc parallel loop gang async(accQueue) \
                 firstprivate(ncols, nactive) \
                    deviceptr(dseg, dsums) \
                      default(none)
  for(std::size_t c=0; c<ncols; ++c) {
    double t=0.0;
#pragma acc loop vector reduction(+:t)
    for(std::size_t e=0; e<nactive; ++e) {
      t+=dseg[c*nactive+e];
    }
    dsums[c]=t;
  }
}

void SecondaryStructureCountACC::cancelLaunch() {
  // the queued work reads the positions of this action, that are about to change
  if( nmine_>0 ) {
    ACC::DeviceScope scope(device_.id());
    acc_wait( accQueue );
  }
}

void SecondaryStructureCountACC::calculate() {
  // unless launch() did it already
  if( !takeLaunch() ) {
    launch();
  }
  const unsigned natoms=getNumberOfAtoms();
  const unsigned stride=comm.Get_size();
  double value=0.0;
  if( nmine_>0 ) {
    ACC::DeviceScope scope(device_.id());
    auto& b=*buf_;
    if( uselist_ ) {
      b.count.copyFromDeviceAsync( &hcount_, accQueue );
      acc_wait( accQueue );
      nactive_=hcount_;
      enqueueSegments();
    }
    if( nactive_>0 ) {
      hsums_.resize( segmentColumns );
      b.sums.copyFromDeviceAsync( hsums_.data(), accQueue );
      b.derivatives.copyFromDeviceAsync( hder_.data(), accQueue );
      acc_wait( accQueue );
      value=hsums_[0];
      for(unsigned j=0; j<9; ++j) {
        hder_[3*natoms+j]=hsums_[1+j];
      }
    }
  }
  if( stride>1 ) {
    comm.Sum( value );
    comm.Sum( hder_ );
  }
  setValueAndDerivatives( value, hder_ );
}

}
}
