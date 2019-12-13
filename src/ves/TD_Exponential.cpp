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

//+PLUMEDOC VES_TARGETDIST TD_EXPONENTIAL
/*
Exponential distribution (static).

Employ a target distribution given by an
[exponential distribution](https://en.wikipedia.org/wiki/Exponential_distribution)
that is defined as
\f[
p(s) =
\lambda e^{-\lambda(s-a)}
\f]
where \f$a\f$ is the minimum of the distribution that is defined on the interval \f$[a,\infty)\f$,
and \f$\lambda>0\f$ is the so-called rate parameter.

The minimum \f$a\f$ is given using the MINIMUM keyword, and the rate parameter \f$\lambda\f$ is given
using the LAMBDA keyword.

This target distribution action is only defined for one dimension, for multiple dimensions
it should be used in combination with \ref TD_PRODUCT_DISTRIBUTION action.

\par Examples

Exponential distribution with \f$a=10.0\f$ and \f$\lambda=0.5\f$
\plumedfile
td: TD_EXPONENTIAL  MINIMUM=-10.0  LAMBDA=0.5
\endplumedfile

The exponential distribution is only defined for one dimension so for multiple
dimensions we have to use it in combination with the \ref TD_PRODUCT_DISTRIBUTION action as shown in
the following example where we have a uniform distribution for argument 1 and
and an exponential distribution for argument 2
\plumedfile
td_uni: TD_UNIFORM

td_exp: TD_EXPONENTIAL  MINIMUM=-10.0  LAMBDA=0.5

td_pd: TD_PRODUCT_DISTRIBUTION DISTRIBUTIONS=td_uni,td_exp
\endplumedfile


*/
//+ENDPLUMEDOC

class TD_Exponential: public TargetDistribution {
  std::vector<double> minima_;
  std::vector<double> lambda_;
public:
  static void registerKeywords(Keywords&);
  explicit TD_Exponential(const ActionOptions& ao);
  double getValue(const std::vector<double>&) const override;
};


PLUMED_REGISTER_ACTION(TD_Exponential,"TD_EXPONENTIAL")


void TD_Exponential::registerKeywords(Keywords& keys) {
  TargetDistribution::registerKeywords(keys);
  keys.add("compulsory","MINIMUM","The minimum of the exponential distribution.");
  keys.add("compulsory","LAMBDA","The \\f$\\lambda\\f$ parameter of the exponential distribution given as positive number.");
  keys.use("WELLTEMPERED_FACTOR");
  keys.use("SHIFT_TO_ZERO");
  keys.use("NORMALIZE");
}


TD_Exponential::TD_Exponential(const ActionOptions& ao):
  PLUMED_VES_TARGETDISTRIBUTION_INIT(ao),
  minima_(0),
  lambda_(0)
{
  parseVector("MINIMUM",minima_);
  parseVector("LAMBDA",lambda_);
  for(unsigned int k=0; k<lambda_.size(); k++) {
    if(lambda_[k] < 0.0) {plumed_merror(getName()+": the value given in LAMBDA should be positive.");}
  }


  setDimension(minima_.size());
  if(getDimension()>1) {plumed_merror(getName()+": only defined for one dimension, for multiple dimensions it should be used in combination with the TD_PRODUCT_DISTRIBUTION action.");}
  if(lambda_.size()!=getDimension()) {plumed_merror(getName()+": the LAMBDA keyword does not match the given dimension in MINIMUM");}
  checkRead();
}


double TD_Exponential::getValue(const std::vector<double>& argument) const {
  double value = 1.0;
  for(unsigned int k=0; k<argument.size(); k++) {
    double arg = (argument[k]-minima_[k])*lambda_[k];
    if(arg<0.0) {plumed_merror(getName()+": the exponential distribution is not defined for values less that ones given in MINIMUM");}
    value *= lambda_[k]*exp(-arg);
  }
  return value;
}



}
}
