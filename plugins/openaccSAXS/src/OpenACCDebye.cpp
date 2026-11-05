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

#include <cmath>
#include <cstddef>
#include <cstring>
#include <vector>

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

/// Debye sum with OpenACC: one thread per row (particle i of this rank, q value k),
/// looping over all the other particles j. Each pair is computed twice (once for
/// row i and once for row j), but there are no atomics, the result does not depend
/// on the scheduling, and the derivative of row i is complete in its thread.
class OpenACCDebye :
  public SAXSDebyeBackend {
  CachedUpload ff_;
  CachedUpload q_;
  OpenACC::memoryManager<double> pos_;
  /// per row (m*numq+k): partial intensity and derivative
  OpenACC::memoryManager<double> partial_;
  OpenACC::memoryManager<double> grad_;
  OpenACC::memoryManager<double> sum_;
  OpenACC::memoryManager<double> w_;
  OpenACC::memoryManager<double> forces_;
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
public:
  explicit OpenACCDebye(int deviceid);
  void compute(const std::vector<Vector>& pos,
               const std::vector<double>& ff,
               const std::vector<double>& q,
               unsigned first, unsigned stride,
               std::vector<double>& sum,
               std::vector<Vector>* deriv) override;
  void contract(const std::vector<double>& w,
                std::vector<Vector>& forces) override;
};

OpenACCDebye::OpenACCDebye(const int deviceid) {
  const int ndevices=acc_get_num_devices(acc_device_not_host);
  // with ACCTARGET=host there are no devices and the kernels run on the host
  if(ndevices>0) {
    if(deviceid>=0) {
      plumed_massert(deviceid<ndevices, "DEVICEID does not identify an available GPU");
      acc_set_device_num(deviceid, acc_device_not_host);
    }
  }
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

  pos_.resize(3*n);
  pos_.copyToDeviceAsync(&pos[0][0], saxsQueue);
  const double* dpos=pos_.devicePtr();
  const double* dff=ff_.update(ff);
  const double* dq=q_.update(q);
  partial_.resize(nrows);
  grad_.resize(3*nrows);
  sum_.resize(numq);
  double* dpartial=partial_.devicePtr();
  double* dgrad=grad_.devicePtr();
  double* dsum=sum_.devicePtr();
  const std::size_t f=first;
  const std::size_t s=stride;

#pragma acc parallel loop gang vector async(saxsQueue) deviceptr(dpos,dff,dq,dpartial,dgrad) \
  firstprivate(n,numq,nrows,f,s) default(none)
  for(std::size_t r=0; r<nrows; ++r) {
    const std::size_t m=r/numq;
    const std::size_t k=r-m*numq;
    const std::size_t i=f+m*s;
    const double xi=dpos[3*i];
    const double yi=dpos[3*i+1];
    const double zi=dpos[3*i+2];
    const double fi=dff[i*numq+k];
    const double qk=dq[k];
    double intensity=0.;
    double gx=0.;
    double gy=0.;
    double gz=0.;
#pragma acc loop seq
    for(std::size_t j=0; j<n; ++j) {
      if(j==i) {
        continue;
      }
      const double dx=dpos[3*j]-xi;
      const double dy=dpos[3*j+1]-yi;
      const double dz=dpos[3*j+2]-zi;
      const double r2=dx*dx+dy*dy+dz*dz;
      const double qr=qk*std::sqrt(r2);
      const double sinc=std::sin(qr)/qr;
      const double weight=fi*dff[j*numq+k];
      intensity+=weight*sinc;
      // d/dr_i of 2 f_i f_j sinc(q r_ij): the pair appears twice in the sum
      const double scale=2.*weight*(sinc-std::cos(qr))/r2;
      gx+=dx*scale;
      gy+=dy*scale;
      gz+=dz*scale;
    }
    dpartial[r]=intensity;
    dgrad[3*r]=gx;
    dgrad[3*r+1]=gy;
    dgrad[3*r+2]=gz;
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
  sum_.copyFromDeviceAsync(hsum_.data(), saxsQueue);
  if(deriv) {
    hgrad_.resize(3*nrows);
    grad_.copyFromDeviceAsync(hgrad_.data(), saxsQueue);
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
  w_.resize(numq);
  w_.copyToDeviceAsync(w.data(), saxsQueue);
  forces_.resize(3*nmine);
  const double* dw=w_.devicePtr();
  const double* dgrad=grad_.devicePtr();
  double* dforces=forces_.devicePtr();

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
  forces_.copyFromDeviceAsync(hforces_.data(), saxsQueue);
  acc_wait(saxsQueue);
  for(std::size_t m=0; m<nmine; ++m) {
    forces[first_+m*stride_]+=Vector(hforces_[3*m],hforces_[3*m+1],hforces_[3*m+2]);
  }
}

SAXSDebyeBackendRegistration<OpenACCDebye> registration("OPENACC");

}

}
}
