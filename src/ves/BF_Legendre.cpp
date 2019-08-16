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

//+PLUMEDOC VES_BASISF BF_LEGENDRE
/*
Legendre polynomials basis functions.

Use as basis functions [Legendre polynomials](https://en.wikipedia.org/wiki/Legendre_polynomials)
\f$P_{n}(x)\f$ defined on a bounded interval.
You need to provide the interval \f$[a,b]\f$
on which the basis functions are to be used, and the order of the
expansion \f$N\f$ (i.e. the highest order polynomial used).
The total number of basis functions is \f$N+1\f$ as the constant \f$P_{0}(x)=1\f$
is also included.
These basis functions should not be used for periodic CVs.

Intrinsically the Legendre polynomials are defined on the interval \f$[-1,1]\f$.
A variable \f$t\f$ in the interval \f$[a,b]\f$ is transformed to a variable \f$x\f$
in the intrinsic interval \f$[-1,1]\f$ by using the transform function
\f[
x(t) = \frac{t-(a+b)/2}
{(b-a)/2}
\f]

The Legendre polynomials are given by the recurrence relation
\f{align}{
P_{0}(x)    &= 1 \\
P_{1}(x)    &= x \\
P_{n+1}(x)  &= \frac{2n+1}{n+1} \, x \, P_{n}(x) -  \frac{n}{n+1} \, P_{n-1}(x)
\f}

The first 6 polynomials are shown below
\image html ves_basisf-legendre.png

The Legendre polynomial are orthogonal over the interval \f$[-1,1]\f$
\f[
\int_{-1}^{1} dx \, P_{n}(x)\, P_{m}(x)  =  \frac{2}{2n+1} \delta_{n,m}
\f]
By using the SCALED keyword the polynomials are scaled by a factor of
\f$ \sqrt{\frac{2n+1}{2}}\f$ such that they are orthonormal to 1.


From the above equation it follows that integral of the basis functions
over the uniform target distribution \f$p_{\mathrm{u}}(x)\f$ are given by
\f[
\int_{-1}^{1} dx \, P_{n}(x) p_{\mathrm{u}}(x) =  \delta_{n,0},
\f]
and thus always zero except for the constant \f$P_{0}(x)=1\f$.


For further mathematical properties of the Legendre polynomials see for example
the [Wikipedia page](https://en.wikipedia.org/wiki/Legendre_polynomials).

\par Examples

Here we employ a Legendre expansion of order 20 over the interval -4.0 to 8.0.
This results in a total number of 21 basis functions.
The label used to identify  the basis function action can then be
referenced later on in the input file.
\plumedfile
bf_leg: BF_LEGENDRE MINIMUM=-4.0 MAXIMUM=8.0 ORDER=20
\endplumedfile

\par Examples

*/
//+ENDPLUMEDOC

class BF_Legendre : public BasisFunctions {
  bool scaled_;
  void setupUniformIntegrals() override;
public:
  static void registerKeywords(Keywords&);
  explicit BF_Legendre(const ActionOptions&);
  void getAllValues(const double, double&, bool&, std::vector<double>&, std::vector<double>&) const override;
};


PLUMED_REGISTER_ACTION(BF_Legendre,"BF_LEGENDRE")


void BF_Legendre::registerKeywords(Keywords& keys) {
  BasisFunctions::registerKeywords(keys);
  keys.addFlag("SCALED",false,"Scale the polynomials such that they are orthonormal to 1.");
}

BF_Legendre::BF_Legendre(const ActionOptions&ao):
  PLUMED_VES_BASISFUNCTIONS_INIT(ao),
  scaled_(false)
{
  parseFlag("SCALED",scaled_); addKeywordToList("SCALED",scaled_);
  setNumberOfBasisFunctions(getOrder()+1);
  setIntrinsicInterval("-1.0","+1.0");
  setNonPeriodic();
  setIntervalBounded();
  setType("Legendre");
  setDescription("Legendre polynomials");
  setLabelPrefix("L");
  setupBF();
  checkRead();
}


void BF_Legendre::getAllValues(const double arg, double& argT, bool& inside_range, std::vector<double>& values, std::vector<double>& derivs) const {
  // plumed_assert(values.size()==numberOfBasisFunctions());
  // plumed_assert(derivs.size()==numberOfBasisFunctions());
  inside_range=true;
  argT=translateArgument(arg, inside_range);
  std::vector<double> derivsT(derivs.size());
  //
  values[0]=1.0;
  derivsT[0]=0.0;
  derivs[0]=0.0;
  values[1]=argT;
  derivsT[1]=1.0;
  derivs[1]=intervalDerivf();
  for(unsigned int i=1; i < getOrder(); i++) {
    double io = static_cast<double>(i);
    values[i+1]  = ((2.0*io+1.0)/(io+1.0))*argT*values[i] - (io/(io+1.0))*values[i-1];
    derivsT[i+1] = ((2.0*io+1.0)/(io+1.0))*(values[i]+argT*derivsT[i])-(io/(io+1.0))*derivsT[i-1];
    derivs[i+1]  = intervalDerivf()*derivsT[i+1];
  }
  if(scaled_) {
    // L0 is also scaled!
    for(unsigned int i=0; i<values.size(); i++) {
      double io = static_cast<double>(i);
      double sf = sqrt(io+0.5);
      values[i] *= sf;
      derivs[i] *= sf;
    }
  }
  if(!inside_range) {for(unsigned int i=0; i<derivs.size(); i++) {derivs[i]=0.0;}}
}


void BF_Legendre::setupUniformIntegrals() {
  setAllUniformIntegralsToZero();
  double L0_int = 1.0;
  if(scaled_) {L0_int = sqrt(0.5);}
  setUniformIntegral(0,L0_int);
}


}
}
