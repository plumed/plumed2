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
#include "plumed/isdb/SAXSDebyeBackend.h"
#include "plumed/tools/OpenACC.h"
#include "plumed/tools/Exception.h"

#include <openacc.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <memory>
#include <vector>

#if defined(PLUMED_OPENACC_GPU) && defined(__linux__)
#include <dlfcn.h>
#endif

namespace PLMD {
namespace isdb {

namespace {

/// All the transfers and the kernels of a call are queued here, the host waits once
constexpr int saxsQueue=1;

/// Uploads data to a device buffer only if it differs from the previous upload
/// (the form factors change only when the solvation is updated, the q values never).
class CachedUpload {
  OpenACC::memoryManager<double> device;
  std::vector<double> host;
  bool valid=false;
public:
  double* update(const std::vector<double>& data) {
    const bool same = valid && host.size()==data.size()
                      && (data.empty() || std::memcmp(host.data(), data.data(), data.size()*sizeof(double))==0);
    if(!same) {
      host=data;
      device.resize(host.size());
      device.copyToDeviceAsync(host.data(), saxsQueue);
      valid=true;
    }
    return device.devicePtr();
  }
};

/// The few functions of the CUDA driver needed to save and restore the CUDA context
/// of the caller (e.g. the MD code), see DeviceScope. They are taken from the
/// driver library that is already loaded (by the caller or by the OpenACC runtime),
/// so that the plugin does not depend on it at link time.
struct CudaDriver {
  typedef int (*getCurrent_t)(void**);
  typedef int (*setCurrent_t)(void*);
  typedef int (*getDevice_t)(int*);
  getCurrent_t getCurrent=nullptr;
  setCurrent_t setCurrent=nullptr;
  getDevice_t getDevice=nullptr;
  bool ok() const {
    return getCurrent && setCurrent && getDevice;
  }
};

/// The CUDA driver, or an empty one if it is not (yet) loaded or not available
const CudaDriver& cudaDriver() {
  static CudaDriver driver;
#if defined(PLUMED_OPENACC_GPU) && defined(__linux__)
  if(!driver.ok()) {
    // never released: the driver stays loaded until the end of the process anyway
    if(void* lib=dlopen("libcuda.so.1",RTLD_LAZY|RTLD_NOLOAD)) {
      CudaDriver d;
      d.getCurrent=reinterpret_cast<CudaDriver::getCurrent_t>(dlsym(lib,"cuCtxGetCurrent"));
      d.setCurrent=reinterpret_cast<CudaDriver::setCurrent_t>(dlsym(lib,"cuCtxSetCurrent"));
      d.getDevice=reinterpret_cast<CudaDriver::getDevice_t>(dlsym(lib,"cuCtxGetDevice"));
      if(d.ok()) {
        driver=d;
      }
    }
  }
#endif
  return driver;
}

/// The device of the CUDA context current in this thread, -1 if there is none
int currentCudaDevice() {
  const auto & cu=cudaDriver();
  void* ctx=nullptr;
  int dev=-1;
  if(cu.ok() && cu.getCurrent(&ctx)==0 && ctx && cu.getDevice(&dev)==0) {
    return dev;
  }
  return -1;
}

/// Makes a device the current one while in scope, and then gives back to the caller
/// its current CUDA context: the MD code may be using another GPU, and it does not
/// expect PLUMED to change it. The OpenACC runtime does not make its context current
/// again before each operation, so this is done here too.
/// The OpenACC device number is not restored: the runtime is linked statically in the
/// plugin, so it is not shared with the caller (and selecting the previous device
/// could create a context on a GPU that is not used at all). With device<0 nothing is done.
class DeviceScope {
  void* prevCtx_=nullptr;
  bool restore_=false;
public:
  explicit DeviceScope(const int device) {
    if(device<0) {
      return;
    }
    const auto & cu=cudaDriver();
    restore_ = cu.ok() && cu.getCurrent(&prevCtx_)==0;
    if(acc_get_device_num(acc_device_not_host)!=device) {
      acc_set_device_num(device,acc_device_not_host);
    }
    if(restore_) {
      // null if the runtime has not initialized the device yet: it makes its
      // context current when it does
      void* ctx=acc_get_current_cuda_context();
      void* cur=nullptr;
      if(ctx && cu.getCurrent(&cur)==0 && cur!=ctx) {
        cu.setCurrent(ctx);
      }
    }
  }
  ~DeviceScope() {
    if(restore_) {
      cudaDriver().setCurrent(prevCtx_);
    }
  }
  DeviceScope(const DeviceScope&)=delete;
  DeviceScope& operator=(const DeviceScope&)=delete;
};

/// Number of q values of a thread (its accumulators stay in registers)
#ifndef PLUMED_SAXS_ACC_QCHUNK
#define PLUMED_SAXS_ACC_QCHUNK 8
#endif
constexpr std::size_t qchunk=PLUMED_SAXS_ACC_QCHUNK;
/// When there are fewer threads than this, the loop over j is split in tiles
/// (whose partial results are then summed in a fixed order)
#ifndef PLUMED_SAXS_ACC_TARGET_THREADS
#define PLUMED_SAXS_ACC_TARGET_THREADS 131072
#endif
constexpr std::size_t targetThreads=PLUMED_SAXS_ACC_TARGET_THREADS;
constexpr std::size_t maxTiles=64;

/// Rows of the Debye sum, see OpenACCDebye. Thread t deals with row i=f+m*s,
/// q values k0..k0+qchunk-1 and particles j of a tile; consecutive threads have
/// consecutive rows and the same chunk and tile, so they read the same j at the same time.
/// With Exact=false sin(q_k r) and cos(q_k r) come from angle addition on the grid
/// q0+k*dq plus a Taylor correction for eps_k=q_k-q0-k*dq (|eps_k r|<=0.1 for all pairs);
/// with Exact=true they are computed directly. The loops over the chunk contain no
/// function calls, so that the compiler can unroll them and keep the per-q
/// accumulators in registers.
template<bool Exact>
void debyeRows(const double* dpos, const double* dff, const double* dq_list, const double* deps,
               double* outpartial, double* outgrad,
               const std::size_t n, const std::size_t numq, const std::size_t nmine,
               const std::size_t nchunk, const std::size_t nrows, const std::size_t tilesize,
               const std::size_t f, const std::size_t s, const std::size_t ntot,
               const double q0, const double dq) {
#pragma acc parallel loop gang vector async(saxsQueue) \
  deviceptr(dpos,dff,dq_list,deps,outpartial,outgrad) \
  firstprivate(n,numq,nmine,nchunk,nrows,tilesize,f,s,ntot,q0,dq) default(none)
  for(std::size_t t=0; t<ntot; ++t) {
    const std::size_t m=t%nmine;
    const std::size_t ct=t/nmine;
    const std::size_t c=ct%nchunk;
    const std::size_t tile=ct/nchunk;
    const std::size_t i=f+m*s;
    const std::size_t k0=c*qchunk;
    const std::size_t nk = numq-k0<qchunk ? numq-k0 : qchunk;
    const std::size_t jbegin=tile*tilesize;
    const std::size_t jend = jbegin+tilesize<n ? jbegin+tilesize : n;
    const double xi=dpos[3*i];
    const double yi=dpos[3*i+1];
    const double zi=dpos[3*i+2];
    const double qbase=q0+k0*dq;
    // beyond the last q value (l>=nk) fi is zero and the other values are those of k0
    double fi[qchunk];
    double qk[qchunk];
    double invq[qchunk];
    double ek[qchunk];
    double intensity[qchunk];
    double gx[qchunk];
    double gy[qchunk];
    double gz[qchunk];
    for(std::size_t l=0; l<qchunk; ++l) {
      const std::size_t k = l<nk ? k0+l : k0;
      fi[l] = l<nk ? dff[i*numq+k] : 0.0;
      qk[l]=dq_list[k];
      invq[l]=1./dq_list[k];
      ek[l]=deps[k];
      intensity[l]=0.;
      gx[l]=0.;
      gy[l]=0.;
      gz[l]=0.;
    }
#pragma acc loop seq
    for(std::size_t j=jbegin; j<jend; ++j) {
      if(j==i) {
        continue;
      }
      const double dx=dpos[3*j]-xi;
      const double dy=dpos[3*j+1]-yi;
      const double dz=dpos[3*j+2]-zi;
      const double r2=dx*dx+dy*dy+dz*dz;
      const double r=std::sqrt(r2);
      const double invr=1./r;
      const double twoinvr2=2.*invr*invr;
      const double* ffj=dff+j*numq+k0;
      double sk[qchunk];
      double ck[qchunk];
      if(Exact) {
        for(std::size_t l=0; l<qchunk; ++l) {
          sk[l]=std::sin(qk[l]*r);
          ck[l]=std::cos(qk[l]*r);
        }
      } else {
        // sin and cos on the grid, from the first q of the chunk by angle addition
        double sg=std::sin(qbase*r);
        double cg=std::cos(qbase*r);
        const double sd=std::sin(dq*r);
        const double cd=std::cos(dq*r);
        for(std::size_t l=0; l<qchunk; ++l) {
          const double x=ek[l]*r;
          const double x2=x*x;
          const double cx=1.-x2*(1./2-x2*(1./24-x2*(1./720-x2*(1./40320))));
          const double sx=x*(1.-x2*(1./6-x2*(1./120-x2*(1./5040-x2*(1./362880)))));
          sk[l]=sg*cx+cg*sx;
          ck[l]=cg*cx-sg*sx;
          const double s2=sg*cd+cg*sd;
          cg=cg*cd-sg*sd;
          sg=s2;
        }
      }
      for(std::size_t l=0; l<qchunk; ++l) {
        const double sinc=sk[l]*invq[l]*invr;
        const double weight=fi[l]*ffj[l<nk ? l : 0];
        intensity[l]+=weight*sinc;
        // d/dr_i of 2 f_i f_j sinc(q r_ij): the pair appears twice in the sum
        const double scale=weight*(sinc-ck[l])*twoinvr2;
        gx[l]+=dx*scale;
        gy[l]+=dy*scale;
        gz[l]+=dz*scale;
      }
    }
    const std::size_t base=tile*nrows+m*numq+k0;
    for(std::size_t l=0; l<qchunk; ++l) {
      if(l<nk) {
        outpartial[base+l]=intensity[l];
        outgrad[3*(base+l)]=gx[l];
        outgrad[3*(base+l)+1]=gy[l];
        outgrad[3*(base+l)+2]=gz[l];
      }
    }
  }
}

/// Debye sum with OpenACC: one thread per (particle i of this rank, chunk of
/// qchunk q values, tile of particles j). Each pair is computed twice (once for
/// row i and once for row j), but there are no atomics, the result does not depend
/// on the scheduling, and the derivative of row i is complete in its thread
/// (or in the sum of its tiles).
/// For every pair the distance is computed once per chunk and, as in SAXS::sincos_q,
/// sin(q_k r) and cos(q_k r) are obtained by angle addition from the first q of the
/// chunk on the uniform grid q0+k*dq, corrected for the deviations eps_k of the
/// actual q values with a Taylor expansion (exact sin and cos, in another version of
/// the kernel, when |eps_k r|>0.1 is possible for some pair or there are less than 4 q values).
///
/// Every call works on the device of the instance (DEVICEID, or by default the one
/// the MD code is using) and leaves the device of the caller unchanged, see DeviceScope.
class OpenACCDebye :
  public SAXSDebyeBackend {
  /// device of this instance, -1 until it is known (or without GPUs)
  int device_=-1;
  /// whether there are GPUs (otherwise the kernels run on the host)
  bool gpu_=false;
  /// the device memory, that must be freed on its device
  struct Buffers {
    CachedUpload ff;
    CachedUpload q;
    CachedUpload eps;
    OpenACC::memoryManager<double> pos;
    /// per row (m*numq+k): partial intensity and derivative
    OpenACC::memoryManager<double> partial;
    OpenACC::memoryManager<double> grad;
    /// the same per tile, when the loop over j is split
    OpenACC::memoryManager<double> tpartial;
    OpenACC::memoryManager<double> tgrad;
    OpenACC::memoryManager<double> sum;
    OpenACC::memoryManager<double> w;
    OpenACC::memoryManager<double> forces;
  };
  std::unique_ptr<Buffers> buf_;
  std::vector<double> heps_;
  std::vector<double> hsum_;
  std::vector<double> hgrad_;
  std::vector<double> hforces_;
  /// sizes and rows of the last compute, used by contract
  std::size_t n_=0;
  std::size_t numq_=0;
  std::size_t nmine_=0;
  std::size_t first_=0;
  std::size_t stride_=1;
  bool computed_=false;
  /// sets device_ to the device of the caller, if it has one
  void followCaller();
public:
  explicit OpenACCDebye(int deviceid);
  ~OpenACCDebye() override;
  OpenACCDebye(const OpenACCDebye&)=delete;
  OpenACCDebye& operator=(const OpenACCDebye&)=delete;
  void compute(const std::vector<Vector>& pos,
               const std::vector<double>& ff,
               const std::vector<double>& q,
               unsigned first, unsigned stride,
               std::vector<double>& sum,
               std::vector<Vector>* deriv) override;
  void contract(const std::vector<double>& w,
                std::vector<Vector>& forces) override;
};

OpenACCDebye::OpenACCDebye(const int deviceid):
  buf_(std::make_unique<Buffers>()) {
  const int ndevices=acc_get_num_devices(acc_device_not_host);
  // with ACCTARGET=host there are no devices and the kernels run on the host
  gpu_ = ndevices>0;
  if(gpu_) {
    if(deviceid>=0) {
      plumed_massert(deviceid<ndevices, "DEVICEID does not identify an available GPU");
      device_=deviceid;
    } else {
      followCaller();
    }
  }
}

void OpenACCDebye::followCaller() {
  const int dev=currentCudaDevice();
  if(dev>=0 && dev<acc_get_num_devices(acc_device_not_host)) {
    device_=dev;
  }
}

OpenACCDebye::~OpenACCDebye() {
  // the buffers are freed here, on their device (members are destroyed after the scope ends)
  DeviceScope scope(device_);
  buf_.reset();
}

void OpenACCDebye::compute(const std::vector<Vector>& pos,
                           const std::vector<double>& ff,
                           const std::vector<double>& q,
                           const unsigned first, const unsigned stride,
                           std::vector<double>& sum,
                           std::vector<Vector>* deriv) {
  const std::size_t n=pos.size();
  const std::size_t numq=q.size();
  plumed_assert(ff.size()==n*numq);
  plumed_assert(sum.size()==numq);
  plumed_assert(stride>0);
  plumed_assert(!deriv || deriv->size()==n*numq);
  n_=n;
  numq_=numq;
  first_=first;
  stride_=stride;
  nmine_ = first<n ? (n-1-first)/stride+1 : 0;
  computed_=true;
  const std::size_t nmine=nmine_;
  const std::size_t nrows=nmine*numq;
  if(nrows==0) {
    return;
  }
  if(gpu_ && device_<0) {
    // the caller may have selected its device after creating this instance;
    // if it has none, the current OpenACC device is used from now on
    followCaller();
    if(device_<0) {
      device_=acc_get_device_num(acc_device_not_host);
    }
  }
  DeviceScope scope(device_);
  auto & b=*buf_;

  // q values as a uniform grid q0+k*dq plus deviations eps_k (as in SAXS)
  const double q0=q[0];
  const double dq = numq>1 ? (q[numq-1]-q[0])/static_cast<double>(numq-1) : 0.0;
  heps_.resize(numq);
  double maxeps=0.0;
  for(std::size_t k=0; k<numq; ++k) {
    heps_[k]=q[k]-(q0+k*dq);
    maxeps=std::max(maxeps,std::fabs(heps_[k]));
  }
  // the angle addition with the Taylor correction is used if |eps_k r|<=0.1 for all the
  // pairs (r is at most the diagonal of the bounding box) and there are at least 4 q values
  Vector lo=pos[0];
  Vector hi=pos[0];
  for(const auto & p : pos) {
    for(unsigned d=0; d<3; ++d) {
      lo[d]=std::min(lo[d],p[d]);
      hi[d]=std::max(hi[d],p[d]);
    }
  }
  const bool exact = numq<4 || maxeps*(hi-lo).modulo()>0.1;

  // work split: chunks of q values and, for small systems, tiles of j
  const std::size_t nchunk=(numq+qchunk-1)/qchunk;
  const std::size_t nthreads=nmine*nchunk;
  std::size_t ntile=(targetThreads+nthreads-1)/nthreads;
  ntile=std::max<std::size_t>(1,std::min(ntile,std::min(maxTiles,n)));
  const std::size_t tilesize=(n+ntile-1)/ntile;

  b.pos.resize(3*n);
  b.pos.copyToDeviceAsync(&pos[0][0], saxsQueue);
  const double* dpos=b.pos.devicePtr();
  const double* dff=b.ff.update(ff);
  const double* dq_list=b.q.update(q);
  const double* deps=b.eps.update(heps_);
  b.partial.resize(nrows);
  b.grad.resize(3*nrows);
  b.sum.resize(numq);
  double* dpartial=b.partial.devicePtr();
  double* dgrad=b.grad.devicePtr();
  double* dsum=b.sum.devicePtr();
  // with a single tile the threads write the rows directly
  double* outpartial=dpartial;
  double* outgrad=dgrad;
  if(ntile>1) {
    b.tpartial.resize(ntile*nrows);
    b.tgrad.resize(3*ntile*nrows);
    outpartial=b.tpartial.devicePtr();
    outgrad=b.tgrad.devicePtr();
  }
  const std::size_t f=first;
  const std::size_t s=stride;
  const std::size_t ntot=ntile*nthreads;

  if(exact) {
    debyeRows<true>(dpos,dff,dq_list,deps,outpartial,outgrad,n,numq,nmine,nchunk,nrows,tilesize,f,s,ntot,q0,dq);
  } else {
    debyeRows<false>(dpos,dff,dq_list,deps,outpartial,outgrad,n,numq,nmine,nchunk,nrows,tilesize,f,s,ntot,q0,dq);
  }

  // sum of the tiles, in a fixed order
  if(ntile>1) {
    const double* dtpartial=b.tpartial.devicePtr();
    const double* dtgrad=b.tgrad.devicePtr();
#pragma acc parallel loop gang vector async(saxsQueue) deviceptr(dtpartial,dtgrad,dpartial,dgrad) \
  firstprivate(nrows,ntile) default(none)
    for(std::size_t r=0; r<nrows; ++r) {
      double p=0.;
      double x=0.;
      double y=0.;
      double z=0.;
#pragma acc loop seq
      for(std::size_t tile=0; tile<ntile; ++tile) {
        const std::size_t tr=tile*nrows+r;
        p+=dtpartial[tr];
        x+=dtgrad[3*tr];
        y+=dtgrad[3*tr+1];
        z+=dtgrad[3*tr+2];
      }
      dpartial[r]=p;
      dgrad[3*r]=x;
      dgrad[3*r+1]=y;
      dgrad[3*r+2]=z;
    }
  }

  // sum of the rows for each q value
#pragma acc parallel loop gang async(saxsQueue) deviceptr(dpartial,dsum) firstprivate(numq,nmine) default(none)
  for(std::size_t k=0; k<numq; ++k) {
    double t=0.;
#pragma acc loop vector reduction(+:t)
    for(std::size_t m=0; m<nmine; ++m) {
      t+=dpartial[m*numq+k];
    }
    dsum[k]=t;
  }

  hsum_.resize(numq);
  b.sum.copyFromDeviceAsync(hsum_.data(), saxsQueue);
  if(deriv) {
    hgrad_.resize(3*nrows);
    b.grad.copyFromDeviceAsync(hgrad_.data(), saxsQueue);
  }
  acc_wait(saxsQueue);

  for(std::size_t k=0; k<numq; ++k) {
    sum[k]+=hsum_[k];
  }
  if(deriv) {
    auto & d=*deriv;
    for(std::size_t m=0; m<nmine; ++m) {
      const std::size_t i=first+m*stride;
      for(std::size_t k=0; k<numq; ++k) {
        const std::size_t r=m*numq+k;
        d[k*n+i]+=Vector(hgrad_[3*r],hgrad_[3*r+1],hgrad_[3*r+2]);
      }
    }
  }
}

void OpenACCDebye::contract(const std::vector<double>& w,
                            std::vector<Vector>& forces) {
  plumed_massert(computed_, "contract called before compute");
  plumed_assert(w.size()==numq_);
  plumed_assert(forces.size()==n_);
  const std::size_t numq=numq_;
  const std::size_t nmine=nmine_;
  if(nmine==0) {
    return;
  }
  DeviceScope scope(device_);
  auto & b=*buf_;
  b.w.resize(numq);
  b.w.copyToDeviceAsync(w.data(), saxsQueue);
  b.forces.resize(3*nmine);
  const double* dw=b.w.devicePtr();
  const double* dgrad=b.grad.devicePtr();
  double* dforces=b.forces.devicePtr();

#pragma acc parallel loop gang vector async(saxsQueue) deviceptr(dw,dgrad,dforces) firstprivate(numq,nmine) default(none)
  for(std::size_t m=0; m<nmine; ++m) {
    double fx=0.;
    double fy=0.;
    double fz=0.;
#pragma acc loop seq
    for(std::size_t k=0; k<numq; ++k) {
      const std::size_t r=m*numq+k;
      fx+=dw[k]*dgrad[3*r];
      fy+=dw[k]*dgrad[3*r+1];
      fz+=dw[k]*dgrad[3*r+2];
    }
    dforces[3*m]=fx;
    dforces[3*m+1]=fy;
    dforces[3*m+2]=fz;
  }

  hforces_.resize(3*nmine);
  b.forces.copyFromDeviceAsync(hforces_.data(), saxsQueue);
  acc_wait(saxsQueue);
  for(std::size_t m=0; m<nmine; ++m) {
    forces[first_+m*stride_]+=Vector(hforces_[3*m],hforces_[3*m+1],hforces_[3*m+2]);
  }
}

SAXSDebyeBackendRegistration<OpenACCDebye> registration("OPENACC");

}

}
}
