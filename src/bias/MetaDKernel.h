/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2026 The plumed team
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
#ifndef __PLUMED_bias_MetaDKernel_h
#define __PLUMED_bias_MetaDKernel_h

#include "core/Value.h"
#include "tools/Grid.h"
#include "tools/Matrix.h"
#include "tools/Tools.h"

#include <cmath>
#include <vector>

namespace PLMD {
namespace bias {

/// A Gaussian hill deposited by METAD.
/// For multivariate hills (ADAPTIVE) sigma contains the upper triangle (row-major) of the
/// inverse covariance matrix, otherwise it contains the widths along each CV.
struct MetaDGaussian {
  bool   multivariate; // this is required to discriminate the one dimensional case
  double height;
  std::vector<double> center;
  std::vector<double> sigma;
  std::vector<double> invsigma;
  /// multivariate hills: the full inverse covariance matrix (row-major), built from sigma
  std::vector<double> metric;
  MetaDGaussian(const bool m,
                const double h,
                const std::vector<double>& c,
                const std::vector<double>& s):
    multivariate(m),
    height(h),
    center(c),
    sigma(s),
    invsigma(s) {
    // to avoid troubles from zero element in flexible hills
    for(unsigned i=0; i<invsigma.size(); ++i) {
      if(std::abs(invsigma[i])>1.e-20) {
        invsigma[i]=1.0/invsigma[i] ;
      } else {
        invsigma[i]=0.0;
      }
    }
    if(multivariate) {
      const unsigned ncv=center.size();
      metric.resize(ncv*ncv);
      unsigned k=0;
      for(unsigned i=0; i<ncv; i++) {
        for(unsigned j=i; j<ncv; j++) {
          metric[i*ncv+j]=metric[j*ncv+i]=sigma[k];
          k++;
        }
      }
    }
  }
};

/// Evaluation of the (possibly stretched) Gaussian hills of METAD as a function of its arguments.
/// Gaussians are truncated at dp2cutoff. With INTERVAL (one dimensional only) the CV is
/// clamped to [lower,upper] and the force outside the interval is zero.
class MetaDKernel {
  /// the arguments of the action, used for (periodic) differences
  std::vector<Value*> args_;
  /// periodicity of the arguments, cached to compute the differences inline
  /// (same arithmetic as Value::difference)
  std::vector<char> periodic_;
  std::vector<double> period_;
  std::vector<double> invperiod_;
  double stretchA_=1.0;
  double stretchB_=0.0;
  bool interval_=false;
  double lower_=-1;
  double upper_=-1;
public:
  /// d2-d1, with the periodicity of argument i (as Value::difference)
  double difference(unsigned i,double d1,double d2) const {
    if(!periodic_[i]) {
      return d2-d1;
    }
    double s=(d2-d1)*invperiod_[i];
    s=Tools::pbc(s);
    return s*period_[i];
  }
  void setArguments(const std::vector<Value*>& args) {
    args_=args;
    periodic_.resize(args.size());
    period_.assign(args.size(),0.0);
    invperiod_.assign(args.size(),0.0);
    for(unsigned i=0; i<args.size(); i++) {
      periodic_[i]=args[i]->isPeriodic();
      if(periodic_[i]) {
        period_[i]=args[i]->getMaxMinusMin();
        invperiod_[i]=1.0/period_[i];
      }
    }
  }
  /// stretched Gaussians: A*exp(-dp2)+B, so that they go to zero at the cutoff
  void setStretch(double A,double B) {
    stretchA_=A;
    stretchB_=B;
  }
  void setInterval(double lower,double upper) {
    interval_=true;
    lower_=lower;
    upper_=upper;
  }
  /// value of hill in cv
  double evaluate(const std::vector<double>& cv, const MetaDGaussian& hill) const;
  /// value of hill in cv; the derivatives are subtracted to der; dp_ is a work array of size ncv
  double evaluateWithDerivatives(const std::vector<double>& cv, const MetaDGaussian& hill, std::vector<double>& der, std::vector<double>& dp_) const;
  /// value of hill at a point, given the differences dp of the point from the center of the hill
  /// (multiplied by invsigma for diagonal hills, plain differences for multivariate hills);
  /// the derivatives are subtracted to der. Outside the cutoff it returns 0, leaves der unchanged and sets
  /// inside=false. Same arithmetic as evaluateWithDerivatives, INTERVAL is not taken into account.
  double evaluateFromDifferences(const double* dp, const MetaDGaussian& hill, double* der, bool& inside) const;
  /// true if INTERVAL is used
  bool hasInterval() const {
    return interval_;
  }
  /// as evaluate, for a diagonal hill given by pointers to its center and inverse widths (no INTERVAL)
  double evaluateDiagonal(const std::vector<double>& cv, const double* center, const double* invsigma, double height) const;
  /// as evaluateWithDerivatives, for a diagonal hill given by pointers to its center and inverse widths (no INTERVAL)
  double evaluateDiagonalWithDerivatives(const std::vector<double>& cv, const double* center, const double* invsigma, double height,
                                         std::vector<double>& der, std::vector<double>& dp_) const;
  /// number of grid points along each direction that are within the cutoff of hill
  std::vector<unsigned> getSupport(const MetaDGaussian& hill, const GridBase& grid) const;
  /// squared distance of cv from the center of hill, in units of the hill widths (d^T*S^-1*d, S covariance);
  /// the hill is truncated where this is 2*dp2cutoff. INTERVAL is not taken into account.
  double getDistance2(const std::vector<double>& cv, const MetaDGaussian& hill) const;
  /// variance of hill along CV i (diagonal of the covariance matrix for multivariate hills)
  double getVariance(const MetaDGaussian& hill, unsigned i) const;
};

inline
double MetaDKernel::evaluate(const std::vector<double>& cv, const MetaDGaussian& hill) const {
  unsigned ncv=cv.size();

  // I use a pointer here because cv is const (and should be const)
  // but when using doInt it is easier to locally replace cv[0] with
  // the upper/lower limit in case it is out of range
  double tmpcv[1];
  const double *pcv=NULL; // pointer to cv
  if(ncv>0) {
    pcv=&cv[0];
  }
  if(interval_) {
    plumed_assert(ncv==1);
    tmpcv[0]=cv[0];
    if(cv[0]<lower_) {
      tmpcv[0]=lower_;
    }
    if(cv[0]>upper_) {
      tmpcv[0]=upper_;
    }
    pcv=&(tmpcv[0]);
  }

  double dp2=0.0;
  if(hill.multivariate) {
    // full inverse covariance matrix
    const double* mymatrix=hill.metric.data();
    for(unsigned i=0; i<ncv; i++) {
      double dp_i=difference(i,hill.center[i],pcv[i]);
      for(unsigned j=i; j<ncv; j++) {
        if(i==j) {
          dp2+=dp_i*dp_i*mymatrix[i*ncv+j]*0.5;
        } else {
          double dp_j=difference(j,hill.center[j],pcv[j]);
          dp2+=dp_i*dp_j*mymatrix[i*ncv+j];
        }
      }
    }
  } else {
    for(unsigned i=0; i<ncv; i++) {
      double dp=difference(i,hill.center[i],pcv[i])*hill.invsigma[i];
      dp2+=dp*dp;
    }
    dp2*=0.5;
  }

  double bias=0.0;
  if(dp2<dp2cutoff) {
    bias=hill.height*(stretchA_*std::exp(-dp2)+stretchB_);
  }

  return bias;
}

inline
double MetaDKernel::evaluateWithDerivatives(const std::vector<double>& cv, const MetaDGaussian& hill, std::vector<double>& der, std::vector<double>& dp_) const {
  unsigned ncv=cv.size();

  // I use a pointer here because cv is const (and should be const)
  // but when using doInt it is easier to locally replace cv[0] with
  // the upper/lower limit in case it is out of range
  const double *pcv=NULL; // pointer to cv
  double tmpcv[1]; // tmp array with cv (to be used with INTERVAL)
  if(ncv>0) {
    pcv=&cv[0];
  }
  if(interval_) {
    plumed_assert(ncv==1);
    tmpcv[0]=cv[0];
    if(cv[0]<lower_) {
      tmpcv[0]=lower_;
    }
    if(cv[0]>upper_) {
      tmpcv[0]=upper_;
    }
    pcv=&(tmpcv[0]);
  }

  bool int_der=false;
  if(interval_) {
    if(cv[0]<lower_ || cv[0]>upper_) {
      int_der=true;
    }
  }

  double dp2=0.0;
  double bias=0.0;
  if(hill.multivariate) {
    // full inverse covariance matrix
    const double* mymatrix=hill.metric.data();
    for(unsigned i=0; i<ncv; i++) {
      dp_[i]=difference(i,hill.center[i],pcv[i]);
      for(unsigned j=i; j<ncv; j++) {
        if(i==j) {
          dp2+=dp_[i]*dp_[i]*mymatrix[i*ncv+j]*0.5;
        } else {
          double dp_j=difference(j,hill.center[j],pcv[j]);
          dp2+=dp_[i]*dp_j*mymatrix[i*ncv+j];
        }
      }
    }
    if(dp2<dp2cutoff) {
      bias=hill.height*std::exp(-dp2);
      if(!int_der) {
        for(unsigned i=0; i<ncv; i++) {
          double tmp=0.0;
          for(unsigned j=0; j<ncv; j++) {
            tmp += dp_[j]*mymatrix[i*ncv+j]*bias;
          }
          der[i]-=tmp*stretchA_;
        }
      } else {
        for(unsigned i=0; i<ncv; i++) {
          der[i]=0.;
        }
      }
      bias=stretchA_*bias+hill.height*stretchB_;
    }
  } else {
    for(unsigned i=0; i<ncv; i++) {
      dp_[i]=difference(i,hill.center[i],pcv[i])*hill.invsigma[i];
      dp2+=dp_[i]*dp_[i];
    }
    dp2*=0.5;
    if(dp2<dp2cutoff) {
      bias=hill.height*std::exp(-dp2);
      if(!int_der) {
        for(unsigned i=0; i<ncv; i++) {
          der[i]-=bias*dp_[i]*hill.invsigma[i]*stretchA_;
        }
      } else {
        for(unsigned i=0; i<ncv; i++) {
          der[i]=0.;
        }
      }
      bias=stretchA_*bias+hill.height*stretchB_;
    }
  }

  return bias;
}

inline
std::vector<unsigned> MetaDKernel::getSupport(const MetaDGaussian& hill, const GridBase& grid) const {
  std::vector<unsigned> nneigh;
  std::vector<double> cutoff;
  unsigned ncv=args_.size();

  // traditional or flexible hill?
  if(hill.multivariate) {
    unsigned k=0;
    Matrix<double> mymatrix(ncv,ncv);
    for(unsigned i=0; i<ncv; i++) {
      for(unsigned j=i; j<ncv; j++) {
        // recompose the full inverse matrix
        mymatrix(i,j)=mymatrix(j,i)=hill.sigma[k];
        k++;
      }
    }
    // the covariance matrix is the inverse of the matrix stored in the hill
    Matrix<double> mycov(ncv,ncv);
    Invert(mymatrix,mycov);
    // the hill is truncated on the ellipsoid 0.5*d^T*mymatrix*d = dp2cutoff,
    // whose bounding box has half-width sqrt(2*dp2cutoff*cov_ii) along CV i
    for(unsigned i=0; i<ncv; i++) {
      cutoff.push_back(std::sqrt(2.0*dp2cutoff)*std::sqrt(mycov(i,i)));
    }
  } else {
    for(unsigned i=0; i<ncv; ++i) {
      cutoff.push_back(std::sqrt(2.0*dp2cutoff)*hill.sigma[i]);
    }
  }

  // number of grid points within the cutoff along CV i; it is never larger than the grid
  // (a support larger than the grid, or not finite for a degenerate hill, covers the whole grid)
  const std::vector<unsigned> nbin=grid.getNbin();
  const std::vector<double> dx=grid.getDx();
  auto npoints=[&](unsigned i) {
    const double n=std::ceil(cutoff[i]/dx[i]);
    return (n>=0.0 && n<nbin[i]) ? static_cast<unsigned>(n) : nbin[i];
  };
  if(interval_) {
    if(hill.center[0]+cutoff[0] > upper_ || hill.center[0]-cutoff[0] < lower_) {
      // in this case, we updated the entire grid to avoid problems
      return nbin;
    } else {
      nneigh.push_back(npoints(0));
      return nneigh;
    }
  } else {
    for(unsigned i=0; i<ncv; i++) {
      nneigh.push_back(npoints(i));
    }
  }

  return nneigh;
}

inline
double MetaDKernel::getDistance2(const std::vector<double>& cv, const MetaDGaussian& hill) const {
  const unsigned ncv=args_.size();
  double dist2=0.0;
  if(hill.multivariate) {
    // hill.sigma holds the upper triangle of the inverse covariance matrix
    std::vector<double> d(ncv);
    for(unsigned i=0; i<ncv; i++) {
      d[i]=difference(i,cv[i],hill.center[i]);
    }
    unsigned k=0;
    for(unsigned i=0; i<ncv; i++) {
      for(unsigned j=i; j<ncv; j++) {
        dist2+=(i==j ? 1.0 : 2.0)*d[i]*d[j]*hill.sigma[k];
        k++;
      }
    }
  } else {
    for(unsigned i=0; i<ncv; i++) {
      const double d=difference(i,cv[i],hill.center[i])/hill.sigma[i];
      dist2+=d*d;
    }
  }
  return dist2;
}

inline
double MetaDKernel::getVariance(const MetaDGaussian& hill, unsigned i) const {
  if(!hill.multivariate) {
    return hill.sigma[i]*hill.sigma[i];
  }
  const unsigned ncv=args_.size();
  Matrix<double> mymatrix(ncv,ncv);
  unsigned k=0;
  for(unsigned j=0; j<ncv; j++) {
    for(unsigned l=j; l<ncv; l++) {
      mymatrix(j,l)=mymatrix(l,j)=hill.sigma[k];
      k++;
    }
  }
  Matrix<double> mycov(ncv,ncv);
  Invert(mymatrix,mycov);
  return mycov(i,i);
}

inline
double MetaDKernel::evaluateFromDifferences(const double* dp, const MetaDGaussian& hill, double* der, bool& inside) const {
  const unsigned ncv=hill.center.size();
  double dp2=0.0;
  double bias=0.0;
  if(hill.multivariate) {
    const double* mymatrix=hill.metric.data();
    for(unsigned i=0; i<ncv; i++) {
      for(unsigned j=i; j<ncv; j++) {
        if(i==j) {
          dp2+=dp[i]*dp[i]*mymatrix[i*ncv+j]*0.5;
        } else {
          dp2+=dp[i]*dp[j]*mymatrix[i*ncv+j];
        }
      }
    }
    inside=(dp2<dp2cutoff);
    if(inside) {
      bias=hill.height*std::exp(-dp2);
      for(unsigned i=0; i<ncv; i++) {
        double tmp=0.0;
        for(unsigned j=0; j<ncv; j++) {
          tmp += dp[j]*mymatrix[i*ncv+j]*bias;
        }
        der[i]-=tmp*stretchA_;
      }
      bias=stretchA_*bias+hill.height*stretchB_;
    }
  } else {
    for(unsigned i=0; i<ncv; i++) {
      dp2+=dp[i]*dp[i];
    }
    dp2*=0.5;
    inside=(dp2<dp2cutoff);
    if(inside) {
      bias=hill.height*std::exp(-dp2);
      for(unsigned i=0; i<ncv; i++) {
        der[i]-=bias*dp[i]*hill.invsigma[i]*stretchA_;
      }
      bias=stretchA_*bias+hill.height*stretchB_;
    }
  }
  return bias;
}

inline
double MetaDKernel::evaluateDiagonal(const std::vector<double>& cv, const double* center, const double* invsigma, double height) const {
  const unsigned ncv=cv.size();
  double dp2=0.0;
  for(unsigned i=0; i<ncv; i++) {
    double dp=difference(i,center[i],cv[i])*invsigma[i];
    dp2+=dp*dp;
  }
  dp2*=0.5;
  double bias=0.0;
  if(dp2<dp2cutoff) {
    bias=height*(stretchA_*std::exp(-dp2)+stretchB_);
  }
  return bias;
}

inline
double MetaDKernel::evaluateDiagonalWithDerivatives(const std::vector<double>& cv, const double* center, const double* invsigma, double height,
    std::vector<double>& der, std::vector<double>& dp_) const {
  const unsigned ncv=cv.size();
  double dp2=0.0;
  for(unsigned i=0; i<ncv; i++) {
    dp_[i]=difference(i,center[i],cv[i])*invsigma[i];
    dp2+=dp_[i]*dp_[i];
  }
  dp2*=0.5;
  double bias=0.0;
  if(dp2<dp2cutoff) {
    bias=height*std::exp(-dp2);
    for(unsigned i=0; i<ncv; i++) {
      der[i]-=bias*dp_[i]*invsigma[i]*stretchA_;
    }
    bias=stretchA_*bias+height*stretchB_;
  }
  return bias;
}

}
}

#endif
