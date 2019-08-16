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

#include "BasisFunctions.h"

#include "core/ActionRegister.h"


namespace PLMD {
namespace ves {

//+PLUMEDOC VES_BASISF BF_POWERS
/*
Polynomial power basis functions.

\attention
__These basis functions should not be used in conventional biasing simulations__.
Instead you should use orthogonal basis functions like Legendre or
Chebyshev polynomials. They are only included for usage in \ref ves_md_linearexpansion
and some special cases.

Basis functions given by polynomial powers defined on a bounded interval.
You need to provide the interval \f$[a,b]\f$
on which the basis functions are to be used, and the order of the
expansion \f$N\f$ (i.e. the highest power used).
The total number of basis functions is \f$N+1\f$ as the constant \f$f_{0}(x)=1\f$
is also included.
These basis functions should not be used for periodic CVs.

The basis functions are given by
\f{align}{
f_{0}(x)    &= 1 \\
f_{1}(x)    &= x \\
f_{2}(x)    &= x^2 \\
& \vdots \\
f_{n}(x)    &= x^n \\
& \vdots \\
f_{N}(x)    &= x^N \\
\f}

Note that these basis functions are __not__ orthogonal. In fact the integral
over the uniform target distribution blows up as the interval is increased.
Therefore they should not be used in conventional biasing simulations.
However, they can be useful for usage with \ref ves_md_linearexpansion.

\par Examples

Here we employ a polynomial power expansion of order 5
over the interval -2.0 to 2.0.
This results in a total number of 6 basis functions.
The label used to identify  the basis function action can then be
referenced later on in the input file.
\plumedfile
BF_POWERS MINIMUM=-2.0 MAXIMUM=2.0 ORDER=5 LABEL=bf_pow
\endplumedfile


*/
//+ENDPLUMEDOC

class BF_Powers : public BasisFunctions {
  double inv_normfactor_;
  void setupLabels() override;
public:
  static void registerKeywords( Keywords&);
  explicit BF_Powers(const ActionOptions&);
  void getAllValues(const double, double&, bool&, std::vector<double>&, std::vector<double>&) const override;
};


PLUMED_REGISTER_ACTION(BF_Powers,"BF_POWERS")


void BF_Powers::registerKeywords(Keywords& keys) {
  BasisFunctions::registerKeywords(keys);
  keys.add("optional","NORMALIZATION","The normalization factor that is used to normalize the basis functions. By default it is 1.0.");
  keys.remove("NUMERICAL_INTEGRALS");
}

BF_Powers::BF_Powers(const ActionOptions&ao):
  PLUMED_VES_BASISFUNCTIONS_INIT(ao)
{
  setNumberOfBasisFunctions(getOrder()+1);
  setIntrinsicInterval(intervalMin(),intervalMax());
  double normfactor_=1.0;
  parse("NORMALIZATION",normfactor_);
  if(normfactor_!=1.0) {addKeywordToList("NORMALIZATION",normfactor_);}
  inv_normfactor_=1.0/normfactor_;
  setNonPeriodic();
  setIntervalBounded();
  setType("polynom_powers");
  setDescription("Polynomial Powers");
  setupBF();
  log.printf("   normalization factor: %f\n",normfactor_);
  checkRead();
}


void BF_Powers::getAllValues(const double arg, double& argT, bool& inside_range, std::vector<double>& values, std::vector<double>& derivs) const {
  inside_range=true;
  argT=checkIfArgumentInsideInterval(arg,inside_range);
  //
  values[0]=1.0;
  derivs[0]=0.0;
  //
  for(unsigned int i=1; i < getNumberOfBasisFunctions(); i++) {
    // double io = static_cast<double>(i);
    // values[i] = pow(argT,io);
    // derivs[i] = io*pow(argT,io-1.0);
    values[i] = argT*values[i-1];
    derivs[i]=values[i-1]+argT*derivs[i-1];
  }
  if(!inside_range) {for(unsigned int i=0; i<derivs.size(); i++) {derivs[i]=0.0;}}
}


void BF_Powers::setupLabels() {
  setLabel(0,"1");
  for(unsigned int i=1; i < getOrder()+1; i++) {
    std::string is; Tools::convert(i,is);
    setLabel(i,"s^"+is);
  }
}

}
}
