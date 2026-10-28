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
  }
};

/// Evaluation of the (possibly stretched) Gaussian hills of METAD as a function of its arguments.
/// Gaussians are truncated at dp2cutoff. With INTERVAL (one dimensional only) the CV is
/// clamped to [lower,upper] and the force outside the interval is zero.
class MetaDKernel {
  /// the arguments of the action, used for (periodic) differences
  std::vector<Value*> args_;
  double stretchA_=1.0;
  double stretchB_=0.0;
  bool interval_=false;
  double lower_=-1;
  double upper_=-1;
  double difference(unsigned i,double d1,double d2) const {
    return args_[i]->difference(d1,d2);
  }
public:
  void setArguments(const std::vector<Value*>& args) {
    args_=args;
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
  /// number of grid points along each direction that are within the cutoff of hill
  std::vector<unsigned> getSupport(const MetaDGaussian& hill, const GridBase& grid) const;
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
    unsigned k=0;
    // recompose the full sigma from the upper diag cholesky
    Matrix<double> mymatrix(ncv,ncv);
    for(unsigned i=0; i<ncv; i++) {
      for(unsigned j=i; j<ncv; j++) {
        mymatrix(i,j)=mymatrix(j,i)=hill.sigma[k]; // recompose the full inverse matrix
        k++;
      }
    }
    for(unsigned i=0; i<ncv; i++) {
      double dp_i=difference(i,hill.center[i],pcv[i]);
      for(unsigned j=i; j<ncv; j++) {
        if(i==j) {
          dp2+=dp_i*dp_i*mymatrix(i,j)*0.5;
        } else {
          double dp_j=difference(j,hill.center[j],pcv[j]);
          dp2+=dp_i*dp_j*mymatrix(i,j);
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
    unsigned k=0;
    // recompose the full sigma from the upper diag cholesky
    Matrix<double> mymatrix(ncv,ncv);
    for(unsigned i=0; i<ncv; i++) {
      for(unsigned j=i; j<ncv; j++) {
        mymatrix(i,j)=mymatrix(j,i)=hill.sigma[k]; // recompose the full inverse matrix
        k++;
      }
    }
    for(unsigned i=0; i<ncv; i++) {
      dp_[i]=difference(i,hill.center[i],pcv[i]);
      for(unsigned j=i; j<ncv; j++) {
        if(i==j) {
          dp2+=dp_[i]*dp_[i]*mymatrix(i,j)*0.5;
        } else {
          double dp_j=difference(j,hill.center[j],pcv[j]);
          dp2+=dp_[i]*dp_j*mymatrix(i,j);
        }
      }
    }
    if(dp2<dp2cutoff) {
      bias=hill.height*std::exp(-dp2);
      if(!int_der) {
        for(unsigned i=0; i<ncv; i++) {
          double tmp=0.0;
          for(unsigned j=0; j<ncv; j++) {
            tmp += dp_[j]*mymatrix(i,j)*bias;
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

  if(interval_) {
    if(hill.center[0]+cutoff[0] > upper_ || hill.center[0]-cutoff[0] < lower_) {
      // in this case, we updated the entire grid to avoid problems
      return grid.getNbin();
    } else {
      nneigh.push_back( static_cast<unsigned>(ceil(cutoff[0]/grid.getDx()[0])) );
      return nneigh;
    }
  } else {
    for(unsigned i=0; i<ncv; i++) {
      nneigh.push_back( static_cast<unsigned>(ceil(cutoff[i]/grid.getDx()[i])) );
    }
  }

  return nneigh;
}

}
}

#endif
