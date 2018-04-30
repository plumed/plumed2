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
the number of equally sized sub-intervalls to be used. On the borders of each
of these sub-intervalls the mean of a Gaussian basis function is placed.

It is also possible to specify the width (i.e. standart deviation) of the
Gaussians using the WIDTH keyword. By default it is set to the sub-intervall
length.

The optimization procedure then adjusts the heigths of the individual Gaussians.

Add stuff here...

*/
//+ENDPLUMEDOC

class BF_Gaussian : public BasisFunctions {
  // width of the Gaussians
  double width_;
  // positions of the means
  std::vector<double> mean_;
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
  keys.add("optional","WIDTH","The width (i.e. standart deviation) of the Gaussian functions. By default it is equal to the sub-intervall size.");
  keys.remove("NUMERICAL_INTEGRALS");
}

BF_Gaussian::BF_Gaussian(const ActionOptions&ao):
  PLUMED_VES_BASISFUNCTIONS_INIT(ao),
  width_((intervalMax()-intervalMin()) / getOrder())
{
  setNumberOfBasisFunctions(getOrder()+2);
  setIntrinsicInterval(intervalMin(),intervalMax());
  parse("WIDTH",width_);
  if(width_ <= 0.0) {plumed_merror("WIDTH should be larger than 0");}
  if(width_ != (intervalMax()-intervalMin())/getOrder()) {addKeywordToList("WIDTH",width_);}
  mean_.reserve(getNumberOfBasisFunctions());
  mean_.push_back(1.0);
  for(unsigned int i=1; i < getNumberOfBasisFunctions(); i++) {
    mean_.push_back(intervalMin()+(1-i)*((intervalMax()-intervalMin())/getOrder()));
  }
  setNonPeriodic();
  setNonOrthogonal();
  setIntervalBounded();
  setType("gaussian_functions");
  setDescription("Gaussian Functions with shifted means that are being optimized in their height");
  setupBF();
  log.printf("   width: %f\n",width_);
  checkRead();
}


void BF_Gaussian::getAllValues(const double arg, double& argT, bool& inside_range, std::vector<double>& values, std::vector<double>& derivs) const {
  inside_range=true;
  argT=checkIfArgumentInsideInterval(arg,inside_range);
  values[0]=1.0;
  derivs[0]=0.0;
  for(unsigned int i=1; i < getNumberOfBasisFunctions(); i++) {
    values[i] = exp(-0.5*pow((argT-mean_[i])/width_,2.0));
    derivs[i] = -values[i] * (argT-mean_[i])/pow(width_,2.0);
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
