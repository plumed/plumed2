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
#include "GridIntegrationWeights.h"

#include "core/ActionRegister.h"
#include "tools/Tools.h"

#include <iostream>



namespace PLMD {
namespace ves {

//+PLUMEDOC VES_TARGETDIST TD_VONMISES
/*
Target distribution given by a sum of Von Mises distributions (static).

Employ a target distribution that is given by a sum where each
term is a product of one-dimensional
[Von Mises distributions](https://en.wikipedia.org/wiki/Von_Mises_distribution),
\f[
p(\mathbf{s}) = \sum_{i} \, w_{i}
\prod_{k}^{d}
\frac{\exp\left(\kappa_{k,i} \, \cos (s_{k}-\mu_{k,i}) \right)}
{2\pi I_{0}(\kappa_{k,i})}
\f]
where \f$(\mu_{1,i},\mu_{2,i},\ldots,\mu_{d,i})\f$
are the centers of the distributions,
\f$(\kappa_{1,i},\kappa_{2,i},\ldots,\kappa_{d,i})\f$
are parameters that determine the extend of each distribution,
and \f$I_{0}(x)\f$ is the modified Bessel function of order 0.
The weights \f$w_{i}\f$ are normalized to 1, \f$\sum_{i}w_{i}=1\f$.

The Von Mises distribution is defined for periodic variables with a
periodicity of \f$2\pi\f$ and is analogous to the Gaussian distribution.
The parameter \f$ \sqrt{1/\kappa}\f$ is comparable to the standard deviation
\f$\sigma\f$ for the Gaussian distribution.

To use this target distribution you need to give the centers
\f$(\mu_{1,i},\mu_{2,i},\ldots,\mu_{d,i})\f$ by
using the numbered CENTER keywords and the "standard deviations"
\f$(\sqrt{1/\kappa_{1,i}},\sqrt{1/\kappa_{2,i}},\ldots,\sqrt{1/\kappa_{d,i}})\f$ using the numbered SIGMA keywords.


\par Examples

Sum of two Von Mises distribution in one dimension that have equal weights
as no weights are given.
\plumedfile
TD_VONMISES ...
 CENTER1=+2.0 SIGMA1=0.6
 CENTER2=-2.0 SIGMA2=0.7
 LABEL=td
... TD_VONMISES
\endplumedfile

Sum of two Von Mises distribution in two dimensions that have different weights.
Note that the weights are automatically normalized to 1 such that
specifying WEIGHTS=1.0,2.0 is equal to specifying WEIGHTS=0.33333,0.66667.
\plumedfile
TD_VONMISES ...
 CENTER1=+2.0,+2.0 SIGMA1=0.6,0.7
 CENTER2=-2.0,+2.0 SIGMA2=0.7,0.6
 WEIGHTS=1.0,2.0
 LABEL=td
... TD_VONMISES
\endplumedfile

*/
//+ENDPLUMEDOC

class TD_VonMises: public TargetDistribution {
  // properties of the Gaussians
  std::vector< std::vector<double> > sigmas_;
  std::vector< std::vector<double> > kappas_;
  std::vector< std::vector<double> > centers_;
  std::vector< std::vector<double> > normalization_;
  std::vector<double> weights_;
  std::vector<double> periods_;
  unsigned int ncenters_;
  double VonMisesDiagonal(const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, const std::vector<double>&) const;
  double getNormalization(const double, const double) const;
public:
  static void registerKeywords(Keywords&);
  explicit TD_VonMises(const ActionOptions& ao);
  double getValue(const std::vector<double>&) const override;
};


PLUMED_REGISTER_ACTION(TD_VonMises,"TD_VONMISES")


void TD_VonMises::registerKeywords(Keywords& keys) {
  TargetDistribution::registerKeywords(keys);
  keys.add("numbered","CENTER","The centers of the Von Mises distributions.");
  keys.add("numbered","SIGMA","The \"standard deviations\" of the Von Mises distributions.");
  keys.add("optional","WEIGHTS","The weights of the Von Mises distributions. Have to be as many as the number of centers given with the numbered CENTER keywords. If no weights are given the distributions are weighted equally. The weights are automatically normalized to 1.");
  keys.add("hidden","PERIODS","The periods for each of the dimensions. By default they are 2*pi for each dimension.");
  keys.use("WELLTEMPERED_FACTOR");
  keys.use("SHIFT_TO_ZERO");
  //keys.use("NORMALIZE");
}


TD_VonMises::TD_VonMises(const ActionOptions& ao):
  PLUMED_VES_TARGETDISTRIBUTION_INIT(ao),
  sigmas_(0),
  centers_(0),
  normalization_(0),
  weights_(0),
  periods_(0),
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
  //
  plumed_massert(centers_.size()==sigmas_.size(),"there has to be an equal amount of CENTER and SIGMA keywords");
  if(centers_.size()==0) {
    plumed_merror(getName()+": CENTER and SIGMA keywords seem to be missing. Note that numbered keywords start at CENTER1 and SIGMA1.");
  }
  //
  setDimension(centers_[0].size());
  ncenters_ = centers_.size();
  //
  // check centers and sigmas
  for(unsigned int i=0; i<ncenters_; i++) {
    if(centers_[i].size()!=getDimension()) {plumed_merror(getName()+": one of the CENTER keyword does not match the given dimension");}
    if(sigmas_[i].size()!=getDimension()) {plumed_merror(getName()+": one of the SIGMA keyword does not match the given dimension");}
  }
  //
  kappas_.resize(sigmas_.size());
  for(unsigned int i=0; i<sigmas_.size(); i++) {
    kappas_[i].resize(sigmas_[i].size());
    for(unsigned int k=0; k<kappas_[i].size(); k++) {
      kappas_[i][k] = 1.0/(sigmas_[i][k]*sigmas_[i][k]);
    }
  }
  //
  parseVector("WEIGHTS",weights_);
  if(weights_.size()==0) {weights_.assign(centers_.size(),1.0);}
  if(centers_.size()!=weights_.size()) {plumed_merror(getName() + ": there has to be as many weights given in WEIGHTS as numbered CENTER keywords");}
  //
  if(periods_.size()==0) {periods_.assign(getDimension(),2*pi);}
  parseVector("PERIODS",periods_);
  if(periods_.size()!=getDimension()) {plumed_merror(getName() + ": the number of values given in PERIODS does not match the dimension of the distribution");}
  //
  double sum_weights=0.0;
  for(unsigned int i=0; i<weights_.size(); i++) {sum_weights+=weights_[i];}
  for(unsigned int i=0; i<weights_.size(); i++) {weights_[i]/=sum_weights;}
  //
  normalization_.resize(ncenters_);
  for(unsigned int i=0; i<ncenters_; i++) {
    normalization_[i].resize(getDimension());
    for(unsigned int k=0; k<getDimension(); k++) {
      normalization_[i][k] = getNormalization(kappas_[i][k],periods_[k]);
    }
  }
  checkRead();
}


double TD_VonMises::getValue(const std::vector<double>& argument) const {
  double value=0.0;
  for(unsigned int i=0; i<ncenters_; i++) {
    value+=weights_[i]*VonMisesDiagonal(argument, centers_[i], kappas_[i],periods_,normalization_[i]);
  }
  return value;
}


double TD_VonMises::VonMisesDiagonal(const std::vector<double>& argument, const std::vector<double>& center, const std::vector<double>& kappa, const std::vector<double>& periods, const std::vector<double>& normalization) const {
  double value = 1.0;
  for(unsigned int k=0; k<argument.size(); k++) {
    double arg = kappa[k]*cos( ((2*pi)/periods[k])*(argument[k]-center[k]) );
    value*=normalization[k]*exp(arg);
  }
  return value;
}


double TD_VonMises::getNormalization(const double kappa, const double period) const {
  //
  std::vector<double> centers(1);
  centers[0] = 0.0;
  std::vector<double> kappas(1);
  kappas[0] = kappa;
  std::vector<double> periods(1);
  periods[0] = period;
  std::vector<double> norm(1);
  norm[0] = 1.0;
  //
  const unsigned int nbins = 1001;
  std::vector<double> points;
  std::vector<double> weights;
  double min = 0.0;
  double max = period;
  GridIntegrationWeights::getOneDimensionalIntegrationPointsAndWeights(points,weights,nbins,min,max);
  //
  double sum = 0.0;
  for(unsigned int l=0; l<nbins; l++) {
    std::vector<double> arg(1); arg[0]= points[l];
    sum += weights[l] * VonMisesDiagonal(arg,centers,kappas,periods,norm);
  }
  return 1.0/sum;
}


}
}
