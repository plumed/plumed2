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
#include "plumed/secondarystructure/SheetCount.h"
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
/// Consecutive GPU threads take segments that are this many positions apart in the list
/// of the segments, so that the threads of a warp, which run together, rarely add to the
/// derivatives of the same atoms (neighbouring segments share a strand)
constexpr std::size_t interleave=32;
/// For each segment, the value of the switching function and its contribution to the virial
constexpr std::size_t segmentColumns=10;
}

/// The GPU version of ANTIBETARMSD_SCALAR and PARABETARMSD_SCALAR (only with TYPE=DRMSD), created by
/// ANTIBETARMSD and PARABETARMSD with USEGPU. A step is done in three kernels:
/// - with STRANDS_CUTOFF, the list of the pairs of strands within the cutoff (one thread per pair);
/// - one thread per segment of that list (or per pair, without STRANDS_CUTOFF): the segment is built
///   and its distance from the references and the switching function are computed with the same
///   functions as on the CPU (sheet::segmentPositions and sheet::segmentSwitch); the derivatives
///   are added to the atoms with atomics, the value and the virial are written for the next kernel;
/// - the sum of the values and of the virial.
/// The derivatives are copied back to the host and used as in SheetCount: there is no backward pass.
/// The pairs, the reference distances and the switching function are copied to the GPU once.
/// The order of the sums is not fixed, so results are reproducible only up to rounding.
class SheetCountACC : public SheetCount {
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
    /// the pairs within STRANDS_CUTOFF, and their number
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
public:
  static void registerKeywords( Keywords& keys );
  explicit SheetCountACC(const ActionOptions&);
  ~SheetCountACC() override;
  SheetCountACC(const SheetCountACC&)=delete;
  SheetCountACC& operator=(const SheetCountACC&)=delete;
  void calculate() override;
};

PLUMED_REGISTER_ACTION(SheetCountACC,"ANTIBETARMSD_SCALAR_ACC")
PLUMED_REGISTER_ACTION(SheetCountACC,"PARABETARMSD_SCALAR_ACC")

void SheetCountACC::registerKeywords( Keywords& keys ) {
  SheetCount::registerKeywords( keys );
}

SheetCountACC::SheetCountACC(const ActionOptions&ao):
  Action(ao),
  SheetCount(ao),
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
}

SheetCountACC::~SheetCountACC() {
  // the device memory is freed here, on its GPU
  ACC::DeviceScope scope(device_.id());
  if( buf_ && buf_->drmsd ) {
    buf_->drmsd->removeFromACCDevice();
  }
  buf_.reset();
}

void SheetCountACC::calculate() {
  const unsigned npairs=pairs.size();
  const unsigned natoms=getNumberOfAtoms();
  const unsigned stride=comm.Get_size();
  const unsigned rank=comm.Get_rank();
  // the pairs rank, rank+stride, ... are computed here, as in SheetCount
  const std::size_t nmine = rank<npairs ? (npairs-1-rank)/stride+1 : 0;
  hder_.assign( 3*natoms+9, 0.0 );
  double value=0.0;
  if( nmine>0 ) {
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
    }
    const SecondaryStructureDRMSDInput<double>* dd=b.drmsd.get();
    const unsigned* dpairs=b.pairs.devicePtr();
    b.positions.resize( 3*natoms );
    b.positions.copyToDeviceAsync( &getPositions()[0][0], accQueue );
    const double* dx=b.positions.devicePtr();
    const Pbc* dpbc = nopbc ? nullptr : b.pbc.update( &getPbc(), 1 );
    const SwitchingFunctionAccelerable* dsw=b.switchingFunction.update( &switchAcc, 1 );
    unsigned* dactive=b.active.devicePtr();

    // the pairs within the cutoff
    const bool usecut = strands_cutoff>0;
    std::size_t nactive=nmine;
    if( usecut ) {
      unsigned* dcount=b.count.devicePtr();
      accZeroAsync( dcount, 1 );
      const double cutoff=strands_cutoff;
#pragma acc parallel loop gang vector async(accQueue) \
                     firstprivate(nmine, rank, stride, cutoff) \
                        deviceptr(dx, dpairs, dpbc, dactive, dcount) \
                          default(none)
      for(std::size_t k=0; k<nmine; ++k) {
        const std::size_t p=rank+k*stride;
        if( sheet::withinCutoff( dx, dpairs[2*p], dpairs[2*p+1], dpbc, cutoff ) ) {
          unsigned slot;
#pragma acc atomic capture
          slot=(*dcount)++;
          dactive[slot]=p;
        }
      }
      unsigned hcount=0;
      b.count.copyFromDeviceAsync( &hcount, accQueue );
      acc_wait( accQueue );
      nactive=hcount;
    }

    if( nactive>0 ) {
      const std::size_t ncols=segmentColumns;
      b.segments.resize( ncols*nactive );
      double* dseg=b.segments.devicePtr();
      b.derivatives.resize( 3*natoms );
      double* dder=b.derivatives.devicePtr();
      accZeroAsync( dder, 3*natoms );
      const std::size_t nint=interleave;
      const std::size_t chunk=(nactive+nint-1)/nint;
      const std::size_t nthreads=nint*chunk;
      const unsigned nstr=nstructures;
      const unsigned nsegatoms=maxAtomsInSegment;
      const unsigned nstrand=sheet::strandAtoms;
      const unsigned nder=sheet::segmentDerivatives;
      // one thread per segment, with all the inner loops sequential
#pragma acc parallel loop gang vector async(accQueue) present(dd[0:1]) \
                     firstprivate(nthreads, chunk, nint, nactive, usecut, rank, stride, \
                                  nstr, nsegatoms, nstrand, nder) \
                        deviceptr(dx, dpairs, dpbc, dsw, dactive, dseg, dder) \
                          default(none)
      for(std::size_t u=0; u<nthreads; ++u) {
        const std::size_t e=(u%nint)*chunk+u/nint;
        if( e>=nactive ) {
          continue;
        }
        const std::size_t p = usecut ? dactive[e] : rank+e*stride;
        const unsigned s1=dpairs[2*p], s2=dpairs[2*p+1];
        std::array<Vector,maxAtomsInSegment> pos;
        std::array<double,2*sheet::segmentDerivatives> der;
        sheet::segmentPositions( dx, s1, s2, dpbc, pos );
        unsigned best=0;
        double g=0.0;
        dseg[e]=sheet::segmentSwitch( *dd, nstr, *dsw, pos, der, best, g );
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
