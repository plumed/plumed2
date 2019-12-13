/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2018 The VES code team
   (see the PEOPLE-VES file at the root of this folder for a list of names)

   See http://www.ves-code.org for more information.

   This file is part of VES code module.

   The VES code module is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   The VES code module is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with the VES code module.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

#include "TargetDistribution.h"

#include "core/ActionRegister.h"


namespace PLMD {
namespace ves {

//+PLUMEDOC VES_TARGETDIST TD_GAUSSIAN
/*
Target distribution given by a sum of Gaussian kernels (static).

Employ a target distribution that is given by a sum of multivariate Gaussian (or normal)
distributions, defined as
\f[
p(\mathbf{s}) = \sum_{i} \, w_{i} \, N(\mathbf{s};\mathbf{\mu}_{i},\mathbf{\Sigma}_{i})
\f]
where \f$\mathbf{\mu}_{i}=(\mu_{1,i},\mu_{2,i},\ldots,\mu_{d,i})\f$
and \f$\mathbf{\Sigma}_{i}\f$ are
the center and the covariance matrix for the \f$i\f$-th Gaussian.
The weights \f$w_{i}\f$ are normalized to 1, \f$\sum_{i}w_{i}=1\f$.

By default the Gaussian distributions are considered as separable into
independent one-dimensional Gaussian distributions. In other words,
the covariance matrix is taken as diagonal
\f$\mathbf{\Sigma}_{i}=(\sigma^2_{1,i},\sigma^2_{2,i},\ldots,\sigma^{2}_{d,i})\f$.
The Gaussian distribution is then written as
\f[
N(\mathbf{s};\mathbf{\mu}_{i},\mathbf{\sigma}_{i}) =
\prod^{d}_{k} \, \frac{1}{\sqrt{2\pi\sigma^2_{d,i}}} \,
\exp\left(
-\frac{(s_{d}-\mu_{d,i})^2}{2\sigma^2_{d,i}}
\right)
\f]
where
\f$\mathbf{\sigma}_{i}=(\sigma_{1,i},\sigma_{2,i},\ldots,\sigma_{d,i})\f$
is the standard deviation.
In this case you need to specify the centers \f$\mathbf{\mu}_{i}\f$ using the
numbered CENTER keywords and the standard deviations \f$\mathbf{\sigma}_{i}\f$
using the numbered SIGMA keywords.

For two arguments it is possible to employ
[bivariate Gaussian kernels](https://en.wikipedia.org/wiki/Multivariate_normal_distribution)
with correlation between arguments, defined as
\f[
N(\mathbf{s};\mathbf{\mu}_{i},\mathbf{\sigma}_{i},\rho_i) =
\frac{1}{2 \pi \sigma_{1,i} \sigma_{2,i} \sqrt{1-\rho_i^2}}
\,
\exp\left(
-\frac{1}{2(1-\rho_i^2)}
\left[
\frac{(s_{1}-\mu_{1,i})^2}{\sigma_{1,i}^2}+
\frac{(s_{2}-\mu_{2,i})^2}{\sigma_{2,i}^2}-
\frac{2 \rho_i (s_{1}-\mu_{1,i})(s_{2}-\mu_{2,i})}{\sigma_{1,i}\sigma_{2,i}}
\right]
\right)
\f]
where \f$\rho_i\f$ is the correlation between \f$s_{1}\f$ and \f$s_{2}\f$
that goes from -1 to 1. In this case the covariance matrix is given as
\f[
\mathbf{\Sigma}=
\left[
\begin{array}{cc}
\sigma^2_{1,i} & \rho_i \sigma_{1,i} \sigma_{2,i} \\
\rho_i \sigma_{1,i} \sigma_{2,i} & \sigma^2_{2,i}
\end{array}
\right]
\f]
The correlation \f$\rho\f$ is given using
the numbered CORRELATION keywords. A value of \f$\rho=0\f$ means
that the arguments are considered as
un-correlated, which is the default behavior.

The Gaussian distributions are always defined with the conventional
normalization factor such that they are normalized to 1 over an unbounded
region. However, in calculation within VES we normally consider bounded
region on which the target distribution is defined. Thus, if the center of
a Gaussian is close to the boundary of the region it can happen that the
tails go outside the region. In that case it might be needed to use the
NORMALIZE keyword to make sure that the target distribution is properly
normalized to 1 over the bounded region. The code will issue a warning
if that is needed.

For periodic CVs it is generally better to use \ref TD_VONMISES "Von Mises"
distributions instead of Gaussian kernels as these distributions properly
account for the periodicity of the CVs.


\par Examples

One single Gaussian kernel in one-dimension.
\plumedfile
td: TD_GAUSSIAN CENTER1=-1.5 SIGMA1=0.8
\endplumedfile

Sum of three Gaussian kernels in two-dimensions with equal weights as
no weights are given.
\plumedfile
TD_GAUSSIAN ...
 CENTER1=-1.5,+1.5 SIGMA1=0.8,0.3
 CENTER2=+1.5,-1.5 SIGMA2=0.3,0.8
 CENTER3=+1.5,+1.5 SIGMA3=0.4,0.4
 LABEL=td
... TD_GAUSSIAN
\endplumedfile

Sum of three Gaussian kernels in two-dimensions which
are weighted unequally. Note that weights are automatically
normalized to 1 so that WEIGHTS=1.0,2.0,1.0 is equal to
specifying WEIGHTS=0.25,0.50,0.25.
\plumedfile
TD_GAUSSIAN ...
 CENTER1=-1.5,+1.5 SIGMA1=0.8,0.3
 CENTER2=+1.5,-1.5 SIGMA2=0.3,0.8
 CENTER3=+1.5,+1.5 SIGMA3=0.4,0.4
 WEIGHTS=1.0,2.0,1.0
 LABEL=td
... TD_GAUSSIAN
\endplumedfile

Sum of two bivariate Gaussian kernels where there is correlation of
\f$\rho_{2}=0.75\f$ between the two arguments for the second Gaussian.
\plumedfile
TD_GAUSSIAN ...
 CENTER1=-1.5,+1.5 SIGMA1=0.8,0.3
 CENTER2=+1.5,-1.5 SIGMA2=0.3,0.8 CORRELATION2=0.75
 LABEL=td
... TD_GAUSSIAN
\endplumedfile





*/
//+ENDPLUMEDOC

class TD_Gaussian: public TargetDistribution {
  std::vector< std::vector<double> > centers_;
  std::vector< std::vector<double> > sigmas_;
  std::vector< std::vector<double> > correlation_;
  std::vector<double> weights_;
  bool diagonal_;
  unsigned int ncenters_;
  double GaussianDiagonal(const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, const bool normalize=true) const;
  double Gaussian2D(const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, const bool normalize=true) const;
public:
  static void registerKeywords(Keywords&);
  explicit TD_Gaussian(const ActionOptions& ao);
  double getValue(const std::vector<double>&) const override;
};


PLUMED_REGISTER_ACTION(TD_Gaussian,"TD_GAUSSIAN")


void TD_Gaussian::registerKeywords(Keywords& keys) {
  TargetDistribution::registerKeywords(keys);
  keys.add("numbered","CENTER","The centers of the Gaussian distributions.");
  keys.add("numbered","SIGMA","The standard deviations of the Gaussian distributions.");
  keys.add("numbered","CORRELATION","The correlation for two-dimensional bivariate Gaussian distributions. Only works for two arguments. The value should be between -1 and 1. If no value is given the Gaussian kernels is considered as un-correlated (i.e. value of 0.0).");
  keys.add("optional","WEIGHTS","The weights of the Gaussian distributions. Have to be as many as the number of centers given with the numbered CENTER keywords. If no weights are given the distributions are weighted equally. The weights are automatically normalized to 1.");
  keys.use("WELLTEMPERED_FACTOR");
  keys.use("SHIFT_TO_ZERO");
  keys.use("NORMALIZE");
}


TD_Gaussian::TD_Gaussian(const ActionOptions& ao):
  PLUMED_VES_TARGETDISTRIBUTION_INIT(ao),
  centers_(0),
  sigmas_(0),
  correlation_(0),
  weights_(0),
  diagonal_(true),
  ncenters_(0)
{
  for(unsigned int i=1;; i++) {
    std::vector<double> tmp_center;
    if(!parseNumberedVector("CENTER",i,tmp_center) ) {break;}
    centers_.push_back(tmp_center);
  }
  for(unsigned int i=1;; i++) {
    std::vector<double> tmp_sigma;
    if(!parseNumberedVector("SIGMA",i,tmp_sigma) ) {break;}
    sigmas_.push_back(tmp_sigma);
  }

  if(centers_.size()==0) {
    plumed_merror(getName()+": CENTER keywords seem to be missing. Note that numbered keywords start at CENTER1.");
  }
  //
  if(centers_.size()!=sigmas_.size()) {
    plumed_merror(getName()+": there has to be an equal amount of CENTER and SIGMA keywords");
  }
  //
  setDimension(centers_[0].size());
  ncenters_ = centers_.size();
  // check centers and sigmas
  for(unsigned int i=0; i<ncenters_; i++) {
    if(centers_[i].size()!=getDimension()) {
      plumed_merror(getName()+": one of the CENTER keyword does not match the given dimension");
    }
    if(sigmas_[i].size()!=getDimension()) {
      plumed_merror(getName()+": one of the SIGMA keyword does not match the given dimension");
    }
  }
  //
  correlation_.resize(ncenters_);

  for(unsigned int i=0; i<ncenters_; i++) {
    std::vector<double> corr;
    parseNumberedVector("CORRELATION",(i+1),corr);
    if(corr.size()>0) {
      diagonal_ = false;
    }
    else {
      corr.assign(1,0.0);
    }
    correlation_[i] = corr;
  }

  if(!diagonal_ && getDimension()!=2) {
    plumed_merror(getName()+": CORRELATION is only defined for two-dimensional Gaussians for now.");
  }
  for(unsigned int i=0; i<correlation_.size(); i++) {
    if(correlation_[i].size()!=1) {
      plumed_merror(getName()+": only one value should be given in CORRELATION");
    }
    for(unsigned int k=0; k<correlation_[i].size(); k++) {
      if(correlation_[i][k] <= -1.0 ||  correlation_[i][k] >= 1.0) {
        plumed_merror(getName()+": values given in CORRELATION should be between -1.0 and 1.0" );
      }
    }
  }
  //
  parseVector("WEIGHTS",weights_);
  if(weights_.size()==0) {weights_.assign(centers_.size(),1.0);}
  if(centers_.size()!=weights_.size()) {
    plumed_merror(getName()+": there has to be as many weights given in WEIGHTS as numbered CENTER keywords");
  }
  //
  double sum_weights=0.0;
  for(unsigned int i=0; i<weights_.size(); i++) {sum_weights+=weights_[i];}
  for(unsigned int i=0; i<weights_.size(); i++) {weights_[i]/=sum_weights;}
  //
  checkRead();
}


double TD_Gaussian::getValue(const std::vector<double>& argument) const {
  double value=0.0;
  if(diagonal_) {
    for(unsigned int i=0; i<ncenters_; i++) {
      value+=weights_[i]*GaussianDiagonal(argument, centers_[i], sigmas_[i]);
    }
  }
  else if(!diagonal_ && getDimension()==2) {
    for(unsigned int i=0; i<ncenters_; i++) {
      value+=weights_[i]*Gaussian2D(argument, centers_[i], sigmas_[i],correlation_[i]);
    }
  }
  return value;
}


double TD_Gaussian::GaussianDiagonal(const std::vector<double>& argument, const std::vector<double>& center, const std::vector<double>& sigma, bool normalize) const {
  double value = 1.0;
  for(unsigned int k=0; k<argument.size(); k++) {
    double arg=(argument[k]-center[k])/sigma[k];
    double tmp_exp = exp(-0.5*arg*arg);
    if(normalize) {tmp_exp/=(sigma[k]*sqrt(2.0*pi));}
    value*=tmp_exp;
  }
  return value;
}


double TD_Gaussian::Gaussian2D(const std::vector<double>& argument, const std::vector<double>& center, const std::vector<double>& sigma, const std::vector<double>& correlation, bool normalize) const {
  double arg1 = (argument[0]-center[0])/sigma[0];
  double arg2 = (argument[1]-center[1])/sigma[1];
  double corr = correlation[0];
  double value = (arg1*arg1 + arg2*arg2 - 2.0*corr*arg1*arg2);
  value *= -1.0 / ( 2.0*(1.0-corr*corr) );
  value = exp(value);
  if(normalize) {
    value /=  2*pi*sigma[0]*sigma[1]*sqrt(1.0-corr*corr);
  }
  return value;
}

}
}
