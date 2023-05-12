/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2021 The VES code team
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
#include "core/PlumedMain.h"


namespace PLMD {
namespace ves {

//+PLUMEDOC VES_BASISF BF_CUBIC_B_SPLINES
/*
Cubic B spline basis functions.

\attention
__These basis functions do not form orthogonal bases. We recommend using wavelets (\ref BF_WAVELETS) instead that do for orthogonal bases__.

A basis using cubic B spline functions according to \cite habermann_multidimensional_2007. See \cite ValssonPampel_Wavelets_2022 for full details.

The mathematical expression of the individual splines is given by
\f{align*}{
  h\left(x\right) =
  \begin{cases}
    \left(2 - \lvert x \rvert\right)^3, & 1 \leq \lvert x \rvert \leq 2\\
    4 - 6\lvert x \rvert^2 + 3 \lvert x \rvert^3,\qquad & \lvert x \rvert \leq 1\\
    0, & \text{elsewhere}.
  \end{cases}
\f}

The full basis consists of equidistant splines at positions \f$\mu_i\f$ which are optimized in their height:
\f{align*}{
  f_i\left(x\right) = h\left(\frac{x-\mu_i}{\sigma}\right)
\f}

Note that the distance between individual splines cannot be chosen freely but is equal to the width: \f$\mu_{i+1} = \mu_{i} + \sigma\f$.


The ORDER keyword of the basis set determines the number of equally sized sub-intervalls to be used.
On the borders of each of these sub-intervalls the mean \f$\mu_i\f$ of a spline function is placed.

The total number of basis functions is \f$\text{ORDER}+4\f$ as the constant \f$f_{0}(x)=1\f$, as well as the two splines with means just outside the interval are also included.

As an example two adjacent basis functions can be seen below.
The full basis consists of shifted splines in the full specified interval.

\image html ves_basisf-splines.png

When the splines are used for a periodic CV (with the PERIODIC keyword),
the sub-intervals are chosen in the same way, but only \f$\text{ORDER}+1\f$ functions
are required to fill it (the ones at the boundary coincide and the ones outside
can be omitted).

To avoid 'blind' optimization of the basis functions outside the currently sampled area, it is often beneficial to use the OPTIMIZATION_THRESHOLD keyword of the VES_LINEAR_EXPANSION (set it to a small value, e.g. 1e-6)

\par Examples
The bias is expanded with cubic B splines in the intervall from 0.0 to 10.0 specifying an order of 20.
This results in 24 basis functions.

\plumedfile
bf: BF_CUBIC_B_SPLINES MINIMUM=0.0 MAXIMUM=10.0 ORDER=20
\endplumedfile

*/
//+ENDPLUMEDOC

class BF_CubicBsplines : public BasisFunctions {
  double spacing_;
  double inv_spacing_;
  double inv_normfactor_;
  double spline(const double, double&) const;
public:
  static void registerKeywords( Keywords&);
  explicit BF_CubicBsplines(const ActionOptions&);
  void getAllValues(const double, double&, bool&, std::vector<double>&, std::vector<double>&) const;
};


PLUMED_REGISTER_ACTION(BF_CubicBsplines,"BF_CUBIC_B_SPLINES")

// See DOI 10.1007/s10614-007-9092-4 for more information;


void BF_CubicBsplines::registerKeywords(Keywords& keys) {
  BasisFunctions::registerKeywords(keys);
  keys.add("optional","NORMALIZATION","the normalization factor that is used to normalize the basis functions by dividing the values. By default it is 2.");
  keys.addFlag("PERIODIC", false, "Use periodic version of basis set.");
  keys.remove("NUMERICAL_INTEGRALS");
}

BF_CubicBsplines::BF_CubicBsplines(const ActionOptions&ao):
  PLUMED_VES_BASISFUNCTIONS_INIT(ao)
{
  log.printf("  Cubic B spline basis functions, see and cite ");
  log << plumed.cite("Pampel and Valsson, J. Chem. Theory Comput. 18, 4127-4141 (2022) - DOI:10.1021/acs.jctc.2c00197");

  setIntrinsicInterval(intervalMin(),intervalMax());
  spacing_=(intervalMax()-intervalMin())/static_cast<double>(getOrder());
  inv_spacing_ = 1.0/spacing_;

  bool periodic = false;
  parseFlag("PERIODIC",periodic);
  if (periodic) {addKeywordToList("PERIODIC",periodic);}

  // 1 constant, getOrder() on interval, 1 (left) + 2 (right) at boundaries if not periodic
  unsigned int num_BFs = periodic ? getOrder()+1U : getOrder()+4U;
  setNumberOfBasisFunctions(num_BFs);

  double normfactor_=2.0;
  parse("NORMALIZATION",normfactor_);
  if(normfactor_!=2.0) {addKeywordToList("NORMALIZATION",normfactor_);}
  inv_normfactor_=1.0/normfactor_;

  periodic ? setPeriodic() : setNonPeriodic();
  setIntervalBounded();
  setType("splines_2nd-order");
  setDescription("Cubic B-splines (2nd order splines)");
  setLabelPrefix("S");
  setupBF();
  log.printf("   normalization factor: %f\n",normfactor_);
  checkRead();
}


void BF_CubicBsplines::getAllValues(const double arg, double& argT, bool& inside_range, std::vector<double>& values, std::vector<double>& derivs) const {
  // plumed_assert(values.size()==numberOfBasisFunctions());
  // plumed_assert(derivs.size()==numberOfBasisFunctions());
  inside_range=true;
  argT=checkIfArgumentInsideInterval(arg,inside_range);
  //
  values[0]=1.0;
  derivs[0]=0.0;
  //
  for(unsigned int i=1; i < getNumberOfBasisFunctions(); i++) {
    double argx = ((argT-intervalMin())*inv_spacing_) - (static_cast<double>(i) - 2.0);
    if(arePeriodic()) { // periodic range of argx is [-intervalRange/spacing,+intervalRange/spacing]
      argx /= intervalRange()*inv_spacing_;
      argx = Tools::pbc(argx);
      argx *= (intervalRange()*inv_spacing_);
    }
    values[i] = spline(argx, derivs[i]);
    derivs[i] *= inv_spacing_;
  }
  if(!inside_range) {for(unsigned int i=0; i<derivs.size(); i++) {derivs[i]=0.0;}}
}


double BF_CubicBsplines::spline(const double arg, double& deriv) const {
  double value=0.0;
  double x=arg;
  // derivative of abs(x);
  double dx = 1.0;
  if(x < 0) {
    x=-x;
    dx = -1.0;
  }
  //
  if(x > 2) {
    value=0.0;
    deriv=0.0;
  }
  else if(x >= 1) {
    value = ((2.0-x)*(2.0-x)*(2.0-x));
    deriv = dx*(-3.0*(2.0-x)*(2.0-x));
    // value=((2.0-x)*(2.0-x)*(2.0-x))/6.0;
    // deriv=-x*x*(2.0-x)*(2.0-x);
  }
  else {
    value = 4.0-6.0*x*x+3.0*x*x*x;
    deriv = dx*(-12.0*x+9.0*x*x);
    // value=x*x*x*0.5-x*x+2.0/3.0;
    // deriv=(3.0/2.0)*x*x-2.0*x;
  }
  value *= inv_normfactor_;
  deriv *= inv_normfactor_;
  return value;
}


}
}
