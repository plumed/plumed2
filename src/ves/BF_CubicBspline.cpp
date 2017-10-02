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

//+PLUMEDOC VES_BASISF BF_CUBIC_B_SPLINES
/*
Cubic B spline basis functions

\par Examples

\par Test

*/
//+ENDPLUMEDOC

class BF_CubicBspline : public BasisFunctions {
  double spacing_;
  double inv_spacing_;
  double inv_normfactor_;
  double spline(const double, double&) const;
public:
  static void registerKeywords( Keywords&);
  explicit BF_CubicBspline(const ActionOptions&);
  void getAllValues(const double, double&, bool&, std::vector<double>&, std::vector<double>&) const;
};


PLUMED_REGISTER_ACTION(BF_CubicBspline,"BF_CUBIC_B_SPLINES")

// See DOI 10.1007/s10614-007-9092-4 for more information;


void BF_CubicBspline::registerKeywords(Keywords& keys) {
  BasisFunctions::registerKeywords(keys);
  keys.add("optional","NORMALIZATION","the normalization factor that is used to normalize the basis functions by dividing the values. By default it is 2.");
  keys.remove("NUMERICAL_INTEGRALS");
}

BF_CubicBspline::BF_CubicBspline(const ActionOptions&ao):
  PLUMED_VES_BASISFUNCTIONS_INIT(ao)
{
  setNumberOfBasisFunctions((getOrder()+3)+1);
  setIntrinsicInterval(intervalMin(),intervalMax());
  spacing_=(intervalMax()-intervalMin())/static_cast<double>(getOrder());
  inv_spacing_ = 1.0/spacing_;
  double normfactor_=2.0;
  parse("NORMALIZATION",normfactor_);
  if(normfactor_!=2.0) {addKeywordToList("NORMALIZATION",normfactor_);}
  inv_normfactor_=1.0/normfactor_;
  setNonPeriodic();
  setIntervalBounded();
  setType("splines_2nd-order");
  setDescription("Cubic B-splines (2nd order splines)");
  setLabelPrefix("S");
  setupBF();
  log.printf("   normalization factor: %f\n",normfactor_);
  checkRead();
}


void BF_CubicBspline::getAllValues(const double arg, double& argT, bool& inside_range, std::vector<double>& values, std::vector<double>& derivs) const {
  // plumed_assert(values.size()==numberOfBasisFunctions());
  // plumed_assert(derivs.size()==numberOfBasisFunctions());
  inside_range=true;
  argT=checkIfArgumentInsideInterval(arg,inside_range);
  //
  values[0]=1.0;
  derivs[0]=0.0;
  //
  for(unsigned int i=1; i < getNumberOfBasisFunctions(); i++) {
    double argx = ((argT-intervalMin())*inv_spacing_) - (static_cast<double>(i)-2.0);
    values[i]  = spline(argx, derivs[i]);
    derivs[i]*=inv_spacing_;
  }
  if(!inside_range) {for(unsigned int i=0; i<derivs.size(); i++) {derivs[i]=0.0;}}
}


double BF_CubicBspline::spline(const double arg, double& deriv) const {
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
