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
__These basis functions are still experimental and should not be used in
conventional biasing simulations__.
Instead you should use orthogonal basis functions like Legendre or Chebyshev
polynomials.

Basis functions given by Gaussian distributions with shifted centers defined on a
bounded interval.
You need to provide the interval \f$[a,b]\f$ on which the bias is to be
expanded.
The ORDER keyword of the basis set \f$N\f$ determines the number of equally sized
sub-intervalls to be used.
On the borders of each of these sub-intervalls the mean \f$\mu\f$ of a Gaussian
basis function is placed:
\f{align}{
  \mu_i = a + (i-1) \frac{b-a}{N}
\f}

The total number of basis functions is \f$N+4\f$ as the constant
\f$f_{0}(x)=1\f$, as well as two additional Gaussians at the Boundaries are also included.

The basis functions are given by
\f{align}{
  f_0(x) &= 1 \\
  f_i(x) &= \exp\left(-\frac{{\left(x-\mu_i\right)}^2}{2\sigma^2}\right)
\f}

It is possible to specify the width \f$\sigma\f$ (i.e. the standard deviation)
of the Gaussians using the WIDTH keyword.
By default it is set to the sub-intervall length.

The optimization procedure then adjusts the heigths of the individual Gaussians.

\par Examples

The bias is expanded with Gaussian functions in the intervall from 0.0 to
10.0 using order 20.
This results in 24 basis functions.

\plumedfile
bfG: BF_GAUSSIANS MINIMUM=0.0 MAXIMUM=10.0 ORDER=20
\endplumedfile

Because it was not specified, the width of the Gaussians is by default
set to the sub-intervall length, i.e.\ \f$\sigma=0.5\f$.
To e.g. enhance the overlap between neighbouring basis functions, it can be
specified explicitely:

\plumedfile
bfG: BF_GAUSSIANS MINIMUM=0.0 MAXIMUM=10.0 ORDER=20 WIDTH=0.7
\endplumedfile

*/
//+ENDPLUMEDOC

class BF_Gaussians : public BasisFunctions {
  /// one over width of the Gaussians
  double inv_sigma_;
  /// positions of the centers
  std::vector<double> centers_;
  void setupLabels() override;
public:
  static void registerKeywords( Keywords&);
  explicit BF_Gaussians(const ActionOptions&);
  void getAllValues(const double, double&, bool&, std::vector<double>&, std::vector<double>&) const override;
};


PLUMED_REGISTER_ACTION(BF_Gaussians,"BF_GAUSSIANS")


void BF_Gaussians::registerKeywords(Keywords& keys) {
  BasisFunctions::registerKeywords(keys);
  keys.add("optional","WIDTH","The width (i.e. standart deviation) of the Gaussian functions. By default it is equal to the sub-intervall size.");
  keys.addFlag("PERIODIC", false, "Use periodic version of basis set.");
  keys.remove("NUMERICAL_INTEGRALS");
}

BF_Gaussians::BF_Gaussians(const ActionOptions&ao):
  PLUMED_VES_BASISFUNCTIONS_INIT(ao)
{
  setIntrinsicInterval(intervalMin(),intervalMax());

  double width = (intervalMax()-intervalMin()) / getOrder();
  parse("WIDTH",width);
  if(width <= 0.0) {plumed_merror("WIDTH should be larger than 0");}
  if(width != (intervalMax()-intervalMin())/getOrder()) {addKeywordToList("WIDTH",width);}
  inv_sigma_ = 1/(width);

  bool periodic = false;
  parseFlag("PERIODIC",periodic);
  if (periodic) {addKeywordToList("PERIODIC",periodic);}

  // 1 constant, getOrder() on interval, 1 (left) + 2 (right) at boundaries if not periodic
  unsigned int num_BFs = periodic ? getOrder()+1U : getOrder()+4U;
  setNumberOfBasisFunctions(num_BFs);

  centers_.push_back(0.0); // constant one
  for(unsigned int i=1; i < getNumberOfBasisFunctions(); i++) {
    centers_.push_back(intervalMin()+(static_cast<int>(i) - 1 - static_cast<int>(!periodic))*(intervalMax()-intervalMin())/getOrder());
  }
  periodic ? setPeriodic() : setNonPeriodic();
  setIntervalBounded();
  setType("gaussian_functions");
  setDescription("Gaussian functions with shifted centers that are being optimized in their height");
  setupBF();
  log.printf("   width: %f\n",width);
  checkRead();
}


void BF_Gaussians::getAllValues(const double arg, double& argT, bool& inside_range, std::vector<double>& values, std::vector<double>& derivs) const {
  inside_range=true;
  argT=checkIfArgumentInsideInterval(arg,inside_range);
  values[0]=1.0;
  derivs[0]=0.0;
  for(unsigned int i=1; i < getNumberOfBasisFunctions(); i++) {
    double dist = argT - centers_[i];
    if(arePeriodic()) { // wrap around similar to MetaD
      dist /= intervalRange();
      dist = Tools::pbc(dist);
      dist *= intervalRange();
    }
    values[i] = exp(-0.5*pow(dist*inv_sigma_,2.0));
    derivs[i] = -values[i] * (dist)*pow(inv_sigma_,2.0);
  }
  if(!inside_range) {for (auto& d: derivs){d=0.0;}}
}


// label according to position of mean
void BF_Gaussians::setupLabels() {
  setLabel(0,"const");
  for(unsigned int i=1; i < getNumberOfBasisFunctions(); i++) {
    std::string is; Tools::convert(centers_[i],is);
    setLabel(i,"m="+is);
  }
}

}
}
