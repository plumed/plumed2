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


namespace PLMD {
namespace ves {

//+PLUMEDOC VES_BASISF BF_CHEBYSHEV
/*
Chebyshev polynomial basis functions.

Use as basis functions [Chebyshev polynomials](https://en.wikipedia.org/wiki/Chebyshev_polynomials)
of the first kind $T_{n}(x)$ defined on a bounded interval.
You need to provide the interval $[a,b]$
on which the basis functions are to be used, and the order of the
expansion $N$ (i.e. the highest order polynomial used).
The total number of basis functions is $N+1$ as the constant $T_{0}(x)=1$
is also included.
These basis functions should not be used for periodic CVs.

Intrinsically the Chebyshev polynomials are defined on the interval $[-1,1]$.
A variable $t$ in the interval $[a,b]$ is transformed to a variable $x$
in the intrinsic interval $[-1,1]$ by using the transform function

$$
x(t) = \frac{t-(a+b)/2}
{(b-a)/2}
$$

The Chebyshev polynomials are given by the recurrence relation

$$
\begin{aligned}
T_{0}(x)    &= 1 \\
T_{1}(x)    &= x \\
T_{n+1}(x)  &= 2 \, x \, T_{n}(x) - T_{n-1}(x)
\end{aligned}
$$

The first 6 polynomials are shown below

![Graph showing first 6 Chebyshev polynomials](figures/ves_basisf-chebyshev.png)

The Chebyshev polynomial are orthogonal over the interval $[-1,1]$
with respect to the weight $\frac{1}{\sqrt{1-x^2}}$

$$
\int_{-1}^{1} dx \, T_{n}(x)\, T_{m}(x) \, \frac{1}{\sqrt{1-x^2}} =
\begin{cases}
0 & n \neq m \\
\pi & n = m = 0 \\
\pi/2 & n = m \neq 0
\end{cases}
$$

For further mathematical properties of the Chebyshev polynomials see for example
the [Wikipedia page](https://en.wikipedia.org/wiki/Chebyshev_polynomials).

## Examples

Here we employ a Chebyshev expansion of order 20 over the interval 0.0 to 10.0.
This results in a total number of 21 basis functions.
The label used to identify  the basis function action can then be
referenced later on in the input file.

```plumed
bfC: BF_CHEBYSHEV MINIMUM=0.0 MAXIMUM=10.0 ORDER=20
```

*/
//+ENDPLUMEDOC


class BF_Chebyshev : public BasisFunctions {
  void setupUniformIntegrals() override;
public:
  static void registerKeywords(Keywords&);
  explicit BF_Chebyshev(const ActionOptions&);
  void getAllValues(const double, double&, bool&, std::vector<double>&, std::vector<double>&) const override;
};


PLUMED_REGISTER_ACTION(BF_Chebyshev,"BF_CHEBYSHEV")


void BF_Chebyshev::registerKeywords(Keywords& keys) {
  BasisFunctions::registerKeywords(keys);
}

BF_Chebyshev::BF_Chebyshev(const ActionOptions&ao):
  PLUMED_VES_BASISFUNCTIONS_INIT(ao) {
  setNumberOfBasisFunctions(getOrder()+1);
  setIntrinsicInterval("-1.0","+1.0");
  setNonPeriodic();
  setIntervalBounded();
  setType("chebyshev-1st-kind");
  setDescription("Chebyshev polynomials of the first kind");
  setLabelPrefix("T");
  setupBF();
  checkRead();
}


void BF_Chebyshev::getAllValues(const double arg, double& argT, bool& inside_range, std::vector<double>& values, std::vector<double>& derivs) const {
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
    values[i+1]  = 2.0*argT*values[i]-values[i-1];
    derivsT[i+1] = 2.0*values[i]+2.0*argT*derivsT[i]-derivsT[i-1];
    derivs[i+1]  = intervalDerivf()*derivsT[i+1];
  }
  if(!inside_range) {
    for(unsigned int i=0; i<derivs.size(); i++) {
      derivs[i]=0.0;
    }
  }
}


void BF_Chebyshev::setupUniformIntegrals() {
  for(unsigned int i=0; i<numberOfBasisFunctions(); i++) {
    double io = i;
    double value = 0.0;
    if(i % 2 == 0) {
      value = -2.0/( pow(io,2.0)-1.0)*0.5;
    }
    setUniformIntegral(i,value);
  }
}


}
}
