/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2017 The VES code team
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

#include "BasisFunctions.h"

#include "core/ActionRegister.h"

namespace PLMD {
namespace ves {

//+PLUMEDOC VES_BASISF BF_GAUSSIAN
/*
Gaussian basis functions.

\attention
__These basis functions should not be used in conventional biasing simulations__.
Instead you should use orthogonal basis functions like Legendre or
Chebyshev polynomials. They are only included for usage in \ref ves_md_linearexpansion
and some special cases.

Basis functions given by Gaussian functions with shifted means defined on a
bounded interval. You need to provide the interval \f$[a,b]\f$ on which the
basis functions are to be used. The order of the expansion \f$N\f$ determines
the position of the means of the Gaussian functions by separating the
intervall into \f$N\f$ parts.
It is also possible to specify the width (i.e. variance) of the Gaussians using
the WIDTH keyword. By default it is set to the sub-intervall length.

The optimization then adjusts the heigths of the individual Gaussians.

Add stuff here...

*/
//+ENDPLUMEDOC

class BF_Gaussian : public BasisFunctions {
  // variance of the Gaussians
  double variance_;
  // positions of the means
  std::vector<double> mean_;
  // normfactor of the Gaussians
  double normfactor_;
  virtual void setupLabels();
public:
  static void registerKeywords( Keywords&);
  explicit BF_Gaussian(const ActionOptions&);
  double getValue(const double, const unsigned int, double&, bool&) const;
  void getAllValues(const double, double&, bool&, std::vector<double>&, std::vector<double>&) const;
};


PLUMED_REGISTER_ACTION(BF_Gaussian,"BF_GAUSSIAN")


void BF_Gaussian::registerKeywords(Keywords& keys) {
  BasisFunctions::registerKeywords(keys);
  keys.add("optional","VARIANCE","The variance (i.e. width) of the Gaussian functions. By default it is equal to the sub-intervall size.");
  keys.remove("NUMERICAL_INTEGRALS");
}

BF_Gaussian::BF_Gaussian(const ActionOptions&ao):
  PLUMED_VES_BASISFUNCTIONS_INIT(ao),
  variance_(intervalRange() / getOrder())
{
  setNumberOfBasisFunctions(getOrder());
  setIntrinsicInterval(intervalMin(),intervalMax());
  if(parse("VARIANCE",variance_)) {addKeywordToList("VARIANCE",variance_);}
  mean_.reserve(getNumberOfBasisFunctions());
  for(unsigned int i=0; i < getNumberOfBasisFunctions; i++)
    mean_.push_back(intervalMin()+(i+0.5)*(intervalRange()/getNumberOfBasisFunctions()));
  normfactor_ = 1/sqrt(2*pi*variance_);
  setNonPeriodic();
  setNonOrthogonal();
  setIntervalBounded();
  setType("gaussian_functions");
  setDescription("Gaussian Functions");
  setupBF();
  log.printf("   variance: %f\n",variance_);
  checkRead();
}


void BF_Gaussian::getAllValues(const double arg, double& argT, bool& inside_range, std::vector<double>& values, std::vector<double>& derivs) const {
  inside_range=true;
  argT=checkIfArgumentInsideInterval(arg,inside_range);
  //
  for(unsigned int i=0; i < getNumberOfBasisFunctions(); i++) {
    // double io = static_cast<double>(i);
    // values[i] = pow(argT,io);
    // derivs[i] = io*pow(argT,io-1.0);
    values[i] = normfactor_ * exp(-pow(argT-mean_[i],2.0)/2*variance_);
    derivs[i] = -values[i] * (argT-mean_[i])/variance_;
  }
  if(!inside_range) {for(unsigned int i=0; i<derivs.size(); i++) {derivs[i]=0.0;}}
}


// label according to positions?
void BF_Gaussian::setupLabels() {
  for(unsigned int i=0; i < getNumberOfBasisFunctions(); i++) {
    std::string is; Tools::convert(mean_[i],is);
    setLabel(i,"s="+is);
  }
}

}
}
