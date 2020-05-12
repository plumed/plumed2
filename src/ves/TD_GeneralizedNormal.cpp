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

//+PLUMEDOC VES_TARGETDIST TD_GENERALIZED_NORMAL
/*
Target distribution given by a sum of generalized normal distributions (static).

Employ a target distribution that is given by a sum where each
term is a product of one-dimensional
[generalized normal distributions](https://en.wikipedia.org/wiki/Generalized_normal_distribution)
(version 1, also know as an exponential power distribution), defined as
\f[
p(\mathbf{s}) = \sum_{i} \, w_{i}
\prod_{k}^{d}
\frac{\beta_{k,i}}{2 \, \alpha_{k,i} \, \Gamma(1/\beta_{k,i})}
\exp\left( -\left\vert \frac{s_{k}-\mu_{k,i}}{\alpha_{k,i}} \right\vert^{\beta_{k,i}} \right)
\f]
where \f$(\mu_{1,i},\mu_{2,i},\ldots,\mu_{d,i})\f$
are the centers of the distributions,
\f$(\alpha_{1,i},\alpha_{2,i},\ldots,\alpha_{d,i})\f$ are the scale
parameters of the distributions,
\f$(\beta_{1,i},\beta_{2,i},\ldots,\beta_{d,i})\f$ are the shape
parameters of the distributions, and \f$\Gamma(x)\f$ is the
gamma function.
The weights \f$w_{i}\f$ are normalized to 1, \f$\sum_{i}w_{i}=1\f$.

Employing \f$\beta=2\f$ results in a
Gaussian (normal) distributions with mean
\f$\mu\f$ and variance \f$\alpha^2/2\f$,
\f$\beta=1\f$ gives the Laplace distribution, and
the limit \f$\beta \to \infty\f$ results in a
uniform  distribution on the interval \f$[\mu-\alpha,\mu+\alpha]\f$.

The centers \f$(\mu_{1,i},\mu_{2,i},\ldots,\mu_{d,i})\f$
are given using the numbered CENTER keywords, the scale
parameters \f$(\alpha_{1,i},\alpha_{2,i},\ldots,\alpha_{d,i})\f$
using the numbered SCALE keywords, and the shape parameters
\f$(\beta_{1,i},\beta_{2,i},\ldots,\beta_{d,i})\f$ using the
numbered SHAPE keywords.
The weights are given using the WEIGHTS keywords, if no weights are
given are all terms weighted equally.

\par Examples

A generalized normal distribution in one-dimensional
\plumedfile
td1: TD_GENERALIZED_NORMAL CENTER1=+20.0  ALPHA1=5.0  BETA1=4.0
\endplumedfile

A sum of two one-dimensional generalized normal distributions
\plumedfile
TD_GENERALIZED_NORMAL ...
 CENTER1=+20.0  ALPHA1=5.0  BETA1=4.0
 CENTER2=-20.0  ALPHA2=5.0  BETA2=3.0
 LABEL=td1
... TD_GENERALIZED_NORMAL
\endplumedfile

A sum of two two-dimensional generalized normal distributions
\plumedfile
TD_GENERALIZED_NORMAL ...
 CENTER1=-20.0,-20.0 ALPHA1=5.0,3.0 BETA1=2.0,4.0
 CENTER2=-20.0,+20.0 ALPHA2=3.0,5.0 BETA2=4.0,2.0
 WEIGHTS=2.0,1.0
 LABEL=td1
... TD_GENERALIZED_NORMAL
\endplumedfile

*/
//+ENDPLUMEDOC

class TD_GeneralizedNormal: public TargetDistribution {
  std::vector< std::vector<double> > centers_;
  std::vector< std::vector<double> > alphas_;
  std::vector< std::vector<double> > betas_;
  std::vector< std::vector<double> > normalization_;
  std::vector<double> weights_;
  unsigned int ncenters_;
  double ExponentialPowerDiagonal(const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, const std::vector<double>&) const;
public:
  static void registerKeywords(Keywords&);
  explicit TD_GeneralizedNormal(const ActionOptions& ao);
  double getValue(const std::vector<double>&) const override;
};


PLUMED_REGISTER_ACTION(TD_GeneralizedNormal,"TD_GENERALIZED_NORMAL")


void TD_GeneralizedNormal::registerKeywords(Keywords& keys) {
  TargetDistribution::registerKeywords(keys);
  keys.add("numbered","CENTER","The center of each generalized normal distribution.");
  keys.add("numbered","ALPHA","The alpha parameters for each generalized normal distribution.");
  keys.add("numbered","BETA","The beta parameters for each generalized normal distribution.");
  keys.add("optional","WEIGHTS","The weights of the generalized normal distribution. By default all are weighted equally.");
  keys.use("WELLTEMPERED_FACTOR");
  keys.use("SHIFT_TO_ZERO");
  keys.use("NORMALIZE");
}


TD_GeneralizedNormal::TD_GeneralizedNormal(const ActionOptions& ao):
  PLUMED_VES_TARGETDISTRIBUTION_INIT(ao),
  centers_(0),
  alphas_(0),
  betas_(0),
  normalization_(0),
  weights_(0),
  ncenters_(0)
{
  for(unsigned int i=1;; i++) {
    std::vector<double> tmp_center;
    if(!parseNumberedVector("CENTER",i,tmp_center) ) {break;}
    centers_.push_back(tmp_center);
  }
  for(unsigned int i=1;; i++) {
    std::vector<double> tmp_alpha;
    if(!parseNumberedVector("ALPHA",i,tmp_alpha) ) {break;}
    for(unsigned int k=0; k<tmp_alpha.size(); k++) {
      if(tmp_alpha[k]<=0.0) {plumed_merror(getName()+": the values given in ALPHA should be positive");}
    }
    alphas_.push_back(tmp_alpha);
  }
  for(unsigned int i=1;; i++) {
    std::vector<double> tmp_beta;
    if(!parseNumberedVector("BETA",i,tmp_beta) ) {break;}
    for(unsigned int k=0; k<tmp_beta.size(); k++) {
      if(tmp_beta[k]<=0.0) {plumed_merror(getName()+": the values given in BETA should be positive");}
    }
    betas_.push_back(tmp_beta);
  }
  //
  if(centers_.size()==0) {
    plumed_merror(getName()+": CENTER keywords seem to be missing. Note that numbered keywords start at CENTER1.");
  }
  //
  if(centers_.size()!=alphas_.size() || centers_.size()!=betas_.size() ) {
    plumed_merror(getName()+": there has to be an equal amount of CENTER, ALPHA, and BETA keywords");
  }
  //
  setDimension(centers_[0].size());
  ncenters_ = centers_.size();
  //
  // check centers and sigmas
  for(unsigned int i=0; i<ncenters_; i++) {
    if(centers_[i].size()!=getDimension()) {
      plumed_merror(getName()+": one of the CENTER keyword does not match the given dimension");
    }
    if(alphas_[i].size()!=getDimension()) {
      plumed_merror(getName()+": one of the ALPHA keyword does not match the given dimension");
    }
    if(betas_[i].size()!=getDimension()) {
      plumed_merror(getName()+": one of the BETA keyword does not match the given dimension");
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
  normalization_.resize(ncenters_);
  for(unsigned int i=0; i<ncenters_; i++) {
    normalization_[i].resize(getDimension());
    for(unsigned int k=0; k<getDimension(); k++) {
      normalization_[i][k] = 0.5*betas_[i][k]/(alphas_[i][k]*tgamma(1.0/betas_[i][k]));
    }
  }
  checkRead();
}


double TD_GeneralizedNormal::getValue(const std::vector<double>& argument) const {
  double value=0.0;
  for(unsigned int i=0; i<ncenters_; i++) {
    value+=weights_[i]*ExponentialPowerDiagonal(argument,centers_[i],alphas_[i],betas_[i],normalization_[i]);
  }
  return value;
}


double TD_GeneralizedNormal::ExponentialPowerDiagonal(const std::vector<double>& argument, const std::vector<double>& center, const std::vector<double>& alpha, const std::vector<double>& beta, const std::vector<double>& normalization) const {
  double value = 1.0;
  for(unsigned int k=0; k<argument.size(); k++) {
    double arg=(std::abs(argument[k]-center[k]))/alpha[k];
    arg = pow(arg,beta[k]);
    value*=normalization[k]*exp(-arg);
  }
  return value;
}



}
}
