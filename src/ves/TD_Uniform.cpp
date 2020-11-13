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

//+PLUMEDOC VES_TARGETDIST TD_UNIFORM
/*
Uniform target distribution (static).

Using this keyword you can define a uniform target distribution which is a
product of one-dimensional distributions \f$p_{k}(s_{k})\f$ that are uniform
over a given interval \f$[a_{k},b_{k}]\f$

MISSING EQUATION TO BE FIXED

The overall distribution is then given as
\f[
p(\mathbf{s}) =
\prod^{d}_{k} p_{k}(s_{k}) =
\begin{cases}
\prod^{d}_{k} \frac{1}{(b_{k}-a_{k})}
& \mathrm{if} \ a_{k} \leq s_{k} \leq b_{k} \ \mathrm{for\ all}\ k \\
\\
0 & \mathrm{otherwise}
\end{cases}
\f]
The distribution is thus uniform inside a rectangular for two arguments
and a cube for a three arguments.

The limits of the intervals \f$ a_{k}\f$ and \f$ b_{k}\f$ are given
with the MINIMA and MAXIMA keywords, respectively. If one or both of
these keywords are missing the code should automatically detect the limits.


It is also possible to use one-dimensional distributions
that go smoothly to zero at the boundaries.
This is done by employing a function with
Gaussian switching functions at the boundaries \f$a_{k}\f$ and \f$b_{k}\f$
\f[
f_{k}(s_{k}) =
\begin{cases}
\exp\left(-\frac{(s_{k}-a_{k})^2}{2 \sigma^2_{a,k}}\right)
& \mathrm{if}\, s_{k} < a_{k} \\
\\
1 & \mathrm{if}\, a_{k} \leq s_{k} \leq b_{k} \\
\\
\exp\left(-\frac{(s_{k}-b_{k})^2}{2 \sigma^2_{b,k}}\right)
& \mathrm{if}\, s_{k} > b_{k}
\end{cases}
\f]
where the standard deviation parameters \f$\sigma_{a,k}\f$
and \f$\sigma_{b,k}\f$ determine how quickly the switching functions
goes to zero.
The overall distribution is then normalized
\f[
p(\mathbf{s}) =
\prod^{d}_{k} p_{k}(s_{k}) =
\prod^{d}_{k} \frac{f(s_{k})}{\int d s_{k} \, f(s_{k})}
\f]
To use this option you need to provide the standard deviation
parameters \f$\sigma_{a,k}\f$ and \f$\sigma_{b,k}\f$ by using the
SIGMA_MINIMA and SIGMA_MAXIMA keywords, respectively. Giving a value of
0.0 means that the boundary is sharp, which is the default behavior.






\par Examples

If one or both of the MINIMA or MAXIMA keywords are missing
the code should automatically detect the limits not given.
Therefore, if we consider a target distribution that is
defined over an interval from 0.0 to 10.0 for the first
argument and from 0.2 to 1.0 for the second argument are
the following example
\plumedfile
td: TD_UNIFORM
\endplumedfile

is equivalent to this one

\plumedfile
TD_UNIFORM ...
 MINIMA=0.0,0.2
 MAXIMA=10.0,1.0
 LABEL=td
 ... TD_UNIFORM
\endplumedfile

and this one

\plumedfile
td: TD_UNIFORM  MAXIMA=10.0,1.0
\endplumedfile

and also this one

\plumedfile
td: TD_UNIFORM MINIMA=0.0,0,2
\endplumedfile


We can also define a target distribution that goes smoothly to zero
at the boundaries of the uniform distribution. In the following
we consider an interval of 0 to 10 for the target distribution.
The following input would result in a target distribution that
would be uniform from 2 to 7 and then smoothly go to zero from
2 to 0 and from 7 to 10.
\plumedfile
TD_UNIFORM ...
 MINIMA=2.0
 MAXIMA=+7.0
 SIGMA_MINIMA=0.5
 SIGMA_MAXIMA=1.0
 LABEL=td
... TD_UNIFORM
\endplumedfile
It is also possible to employ a smooth switching function for just one
of the boundaries as shown here where the target distribution
would be uniform from 0 to 7 and then smoothly go to zero from 7 to 10.
\plumedfile
TD_UNIFORM ...
 MAXIMA=+7.0
 SIGMA_MAXIMA=1.0
 LABEL=td
... TD_UNIFORM
\endplumedfile
Furthermore, it is possible to employ a sharp boundary by
using
\plumedfile
TD_UNIFORM ...
 MAXIMA=+7.0
 SIGMA_MAXIMA=0.0
 LABEL=td
... TD_UNIFORM
\endplumedfile
or
\plumedfile
td: TD_UNIFORM MAXIMA=+7.0
\endplumedfile


*/
//+ENDPLUMEDOC

class TD_Uniform : public TargetDistribution {
  std::vector<double> minima_;
  std::vector<double> maxima_;
  std::vector<double> sigma_min_;
  std::vector<double> sigma_max_;
  double GaussianSwitchingFunc(const double, const double, const double) const;
  void setupAdditionalGrids(const std::vector<Value*>&, const std::vector<std::string>&, const std::vector<std::string>&, const std::vector<unsigned int>&) override;
public:
  static void registerKeywords( Keywords&);
  explicit TD_Uniform(const ActionOptions& ao);
  double getValue(const std::vector<double>&) const override;
};


PLUMED_REGISTER_ACTION(TD_Uniform,"TD_UNIFORM")


void TD_Uniform::registerKeywords(Keywords& keys) {
  TargetDistribution::registerKeywords(keys);
  keys.add("optional","MINIMA","The minimum of the intervals where the target distribution is taken as uniform. You should give one value for each argument.");
  keys.add("optional","MAXIMA","The maximum of the intervals where the target distribution is taken as uniform. You should give one value for each argument.");
  keys.add("optional","SIGMA_MINIMA","The standard deviation parameters of the Gaussian switching functions for the minima of the intervals. You should give one value for each argument. Value of 0.0 means that switch is done without a smooth switching function, this is the default behavior.");
  keys.add("optional","SIGMA_MAXIMA","The standard deviation parameters of the Gaussian switching functions for the maximum of the intervals. You should give one value for each argument. Value of 0.0 means that switch is done without a smooth switching function, this is the default behavior.");
}


TD_Uniform::TD_Uniform(const ActionOptions& ao):
  PLUMED_VES_TARGETDISTRIBUTION_INIT(ao),
  minima_(0),
  maxima_(0),
  sigma_min_(0),
  sigma_max_(0)
{
  parseVector("MINIMA",minima_);
  parseVector("MAXIMA",maxima_);

  parseVector("SIGMA_MINIMA",sigma_min_);
  parseVector("SIGMA_MAXIMA",sigma_max_);
  if(minima_.size()==0 && sigma_min_.size()>0) {plumed_merror(getName()+": you cannot give SIGMA_MINIMA if MINIMA is not given");}
  if(maxima_.size()==0 && sigma_max_.size()>0) {plumed_merror(getName()+": you cannot give SIGMA_MAXIMA if MAXIMA is not given");}

  if(minima_.size()>0 && maxima_.size()>0) {
    // both MINIMA and MAXIMA given, do all checks
    if(minima_.size()!=maxima_.size()) {plumed_merror(getName()+": MINIMA and MAXIMA do not have the same number of values.");}
    setDimension(minima_.size());
    for(unsigned int k=0; k<getDimension(); k++) {
      if(minima_[k]>maxima_[k]) {
        plumed_merror(getName()+": error in MINIMA and MAXIMA keywords, one of the MINIMA values is larger than the corresponding MAXIMA values");
      }
    }
  }
  else if(minima_.size()>0 && maxima_.size()==0) {
    // only MINIMA given, MAXIMA assigned later on.
    setDimension(minima_.size());
  }
  else if(maxima_.size()>0 && minima_.size()==0) {
    // only MAXIMA given, MINIMA assigned later on.
    setDimension(maxima_.size());
  }
  else if(maxima_.size()==0 && minima_.size()==0) {
    // neither MAXIMA nor MINIMA givenm, both assigned later on.
    setDimension(0);
  }

  if(sigma_min_.size()==0) {sigma_min_.assign(getDimension(),0.0);}
  if(sigma_max_.size()==0) {sigma_max_.assign(getDimension(),0.0);}
  if(sigma_min_.size()!=getDimension()) {plumed_merror(getName()+": SIGMA_MINIMA has the wrong number of values");}
  if(sigma_max_.size()!=getDimension()) {plumed_merror(getName()+": SIGMA_MAXIMA has the wrong number of values");}
  //
  setForcedNormalization();
  checkRead();
}


void TD_Uniform::setupAdditionalGrids(const std::vector<Value*>& arguments, const std::vector<std::string>& min, const std::vector<std::string>& max, const std::vector<unsigned int>& nbins) {

  if(minima_.size()==0) {
    minima_.assign(getDimension(),0.0);
    for(unsigned int k=0; k<getDimension(); k++) {Tools::convert(min[k],minima_[k]);}
  }

  if(maxima_.size()==0) {
    maxima_.assign(getDimension(),0.0);
    for(unsigned int k=0; k<getDimension(); k++) {Tools::convert(max[k],maxima_[k]);}
  }

}


double TD_Uniform::getValue(const std::vector<double>& argument) const {
  //
  double value = 1.0;
  for(unsigned int k=0; k<getDimension(); k++) {
    double tmp;
    if(argument[k] < minima_[k]) {
      tmp = GaussianSwitchingFunc(argument[k],minima_[k],sigma_min_[k]);
    }
    else if(argument[k] > maxima_[k]) {
      tmp = GaussianSwitchingFunc(argument[k],maxima_[k],sigma_max_[k]);
    }
    else {
      tmp = 1.0;
    }
    value *= tmp;
  }
  return value;
}

inline
double TD_Uniform::GaussianSwitchingFunc(const double argument, const double center, const double sigma) const {
  if(sigma>0.0) {
    double arg=(argument-center)/sigma;
    return exp(-0.5*arg*arg);
  }
  else {
    return 0.0;
  }
}






}
}
