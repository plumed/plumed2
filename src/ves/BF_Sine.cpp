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

//+PLUMEDOC VES_BASISF BF_SINE
/*
Fourier sine basis functions.

Use as basis functions Fourier sine series defined on a periodic interval.
You need to provide the periodic interval $[a,b]$
on which the basis functions are to be used, and the order of the
expansion $N$ (i.e. the highest Fourier sine mode used).
The total number of basis functions is $N+1$ as
the constant $f_{0}(x)=1$ is also included.
These basis functions should only be used for periodic CVs.
They can be useful if the periodic function being expanded is an
odd function, i.e. $F(-x)=-F(x)$.

The Fourier sine basis functions are given by

$$
\begin{aligned}
f_{0}(x)    &= 1 \\
f_{1}(x)    &= sin(\frac{2\pi }{P} x) \\
f_{2}(x)    &= sin(2 \cdot \frac{2\pi}{P} x) \\
f_{3}(x)    &= sin(3 \cdot \frac{2\pi}{P} x) \\
& \vdots \\
f_{n}(x) &= sin(n \cdot \frac{2\pi}{P} x) \\
& \vdots \\
f_{N}(x)   &= sin(N \cdot \frac{2\pi}{P} x) \\
\end{aligned}
$$

where $P=(b-a)$ is the periodicity of the interval.
They are orthogonal over the interval $[a,b]$

$$
\int_{a}^{b} dx \, f_{n}(x)\, f_{m}(x)  =
\begin{cases}
0 & n \neq m \\
(b-a) & n = m = 0 \\
(b-a)/2 & n = m \neq 0
\end{cases}.
$$

## Examples

Here we employ a Fourier sine expansion of order 10 over the periodic interval
$-\pi$ to $+\pi$.
This results in a total number of 11 basis functions.
The label used to identify  the basis function action can then be
referenced later on in the input file.

```plumed
BF_SINE MINIMUM=-pi MAXIMUM=+pi ORDER=10 LABEL=bfS
```

*/
//+ENDPLUMEDOC

class BF_Sine : public BasisFunctions {
  void setupLabels() override;
  void setupUniformIntegrals() override;
public:
  static void registerKeywords(Keywords&);
  explicit BF_Sine(const ActionOptions&);
  void getAllValues(const double, double&, bool&, std::vector<double>&, std::vector<double>&) const override;
};


PLUMED_REGISTER_ACTION(BF_Sine,"BF_SINE")


void BF_Sine::registerKeywords(Keywords& keys) {
  BasisFunctions::registerKeywords(keys);
}


BF_Sine::BF_Sine(const ActionOptions&ao):
  PLUMED_VES_BASISFUNCTIONS_INIT(ao) {
  setNumberOfBasisFunctions(getOrder()+1);
  setIntrinsicInterval("-pi","+pi");
  setPeriodic();
  setIntervalBounded();
  setType("trigonometric_sin");
  setDescription("Sine");
  setupBF();
  checkRead();
}


void BF_Sine::getAllValues(const double arg, double& argT, bool& inside_range, std::vector<double>& values, std::vector<double>& derivs) const {
  // plumed_assert(values.size()==numberOfBasisFunctions());
  // plumed_assert(derivs.size()==numberOfBasisFunctions());
  inside_range=true;
  argT=translateArgument(arg, inside_range);
  values[0]=1.0;
  derivs[0]=0.0;
  for(unsigned int i=1; i < getOrder()+1; i++) {
    double io = i;
    double cos_tmp = cos(io*argT);
    double sin_tmp = sin(io*argT);
    values[i] = sin_tmp;
    derivs[i] = io*cos_tmp*intervalDerivf();
  }
  if(!inside_range) {
    for(unsigned int i=0; i<derivs.size(); i++) {
      derivs[i]=0.0;
    }
  }
}


void BF_Sine::setupLabels() {
  setLabel(0,"1");
  for(unsigned int i=1; i < getOrder()+1; i++) {
    std::string is;
    Tools::convert(i,is);
    setLabel(i,"sin("+is+"*s)");
  }
}


void BF_Sine::setupUniformIntegrals() {
  setAllUniformIntegralsToZero();
  setUniformIntegral(0,1.0);
}


}
}
