/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2020 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "function/FunctionTemplateBase.h"
#include "function/FunctionShortcut.h"
#include "function/FunctionOfMatrix.h"
#include "core/ActionRegister.h"

#include <string>
#include <cmath>

namespace PLMD {
namespace symfunc {

//+PLUMEDOC MCOLVAR FCCUBIC_FUNC
/*
Measure how similar the environment around atoms is to that found in a FCC structure.

\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC MCOLVAR FCCUBIC_FUNC_MATRIX
/*
Measure how similar the environment around atoms is to that found in a FCC structure.

\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC MCOLVAR FCCUBIC
/*
Measure how similar the environment around atoms is to that found in a FCC structure.

This shortcut is an example of a [COORDINATION_SHELL_AVERAGE](COORDINATION_SHELL_AVERAGE.md),
which we can use to measure how similar the environment around atom $i$ is to an fcc structure.
The function that is used to make this determination is as follows:

$$
s_i = \frac{ \sum_{i \ne j} \sigma(r_{ij}) \left\{ a\left[ \frac{(x_{ij}y_{ij})^4 + (x_{ij}z_{ij})^4 + (y_{ij}z_{ij})^4}{r_{ij}^8} - \frac{\alpha (x_{ij}y_{ij}z_{ij})^4}{r_{ij}^{12}} \right] + b \right\} }{ \sum_{i \ne j} \sigma(r_{ij}) }
$$


In this expression $x_{ij}$, $y_{ij}$ and $z_{ij}$ are the $x$, $y$ and $z$ components of the vector connecting atom $i$ to
atom $j$ and $r_{ij}$ is the magnitude of this vector.  $\sigma(r_{ij})$ is a switching function that acts on the distance between
atom $i$ and atom $j$ and its inclusion in the numerator and the denominator of the above expression as well as the fact that we are summing
over all of the other atoms in the system ensures that we are calculating an average
of the function of $x_{ij}$, $y_{ij}$ and $z_{ij}$ for the atoms in the first coordination sphere around atom $i$.  Lastly, $\alpha$
is a parameter that can be set by the user, which by default is equal to three.  The values of $a$ and $b$ are calculated from $\alpha$ using:

$$
a = \frac{ 80080}{ 2717 + 16 \alpha} \qquad \textrm{and} \qquad b = \frac{ 16(\alpha - 143) }{2717 + 16\alpha}
$$

This action was been used in all the articles in the bibliography. We thus wrote an explict
action to calculate it in PLUMED instead of using a shortcut as we did for [SIMPLECUBIC](SIMPLECUBIC.md) so that we could get
good computational performance.

The following input calculates the FCCUBIC parameter for the 64 atoms in the system
and then calculates and prints the average value for this quantity.

```plumed
d: FCCUBIC SPECIES=1-64 SWITCH={RATIONAL D_0=3.0 R_0=1.5} MEAN
PRINT ARG=d.* FILE=colv
```

Notice that you can you can rotate the bond vectors before computing the
function in the above expression as is discussed in the documentation for [COORDINATION_SHELL_FUNCTION](COORDINATION_SHELL_FUNCTION.md)

*/
//+ENDPLUMEDOC


class Fccubic : public function::FunctionTemplateBase {
private:
  double alpha, a1, b1;
public:
  void registerKeywords( Keywords& keys ) override;
  void read( ActionWithArguments* action ) override;
  void calc( const ActionWithArguments* action, const std::vector<double>& args, std::vector<double>& vals, Matrix<double>& derivatives ) const override;
};

typedef function::FunctionShortcut<Fccubic> FccubicShortcut;
PLUMED_REGISTER_ACTION(FccubicShortcut,"FCCUBIC_FUNC")
typedef function::FunctionOfMatrix<Fccubic> MatrixFccubic;
PLUMED_REGISTER_ACTION(MatrixFccubic,"FCCUBIC_FUNC_MATRIX")

void Fccubic::registerKeywords( Keywords& keys ) {
  keys.add("compulsory","ALPHA","3.0","The alpha parameter of the angular function");
  keys.setValueDescription("matrix","a function that measures the similarity with an fcc environment");
}

void Fccubic::read( ActionWithArguments* action ) {
  // Scaling factors such that '1' corresponds to fcc lattice
  // and '0' corresponds to isotropic (liquid)
  parse(action,"ALPHA",alpha);
  a1 = 80080. / (2717. + 16*alpha);
  b1 = 16.*(alpha-143)/(2717+16*alpha);
  action->log.printf("  setting alpha paramter equal to %f \n",alpha);
}

void Fccubic::calc( const ActionWithArguments* action, const std::vector<double>& args, std::vector<double>& vals, Matrix<double>& derivatives ) const {
  double x2 = args[0]*args[0];
  double x4 = x2*x2;

  double y2 = args[1]*args[1];
  double y4 = y2*y2;

  double z2 = args[2]*args[2];
  double z4 = z2*z2;

  double d2 = x2 + y2 + z2;
  double r8 = pow( d2, 4 );
  double r12 = pow( d2, 6 );

  double tmp = ((x4*y4)+(x4*z4)+(y4*z4))/r8-alpha*x4*y4*z4/r12;

  double t0 = (x2*y4+x2*z4)/r8-alpha*x2*y4*z4/r12;
  double t1 = (y2*x4+y2*z4)/r8-alpha*y2*x4*z4/r12;
  double t2 = (z2*x4+z2*y4)/r8-alpha*z2*x4*y4/r12;
  double t3 = (2*tmp-alpha*x4*y4*z4/r12)/d2;

  derivatives(0,0)=4*a1*args[0]*(t0-t3);
  derivatives(0,1)=4*a1*args[1]*(t1-t3);
  derivatives(0,2)=4*a1*args[2]*(t2-t3);

  // Set the value and the derivatives
  vals[0] = (a1*tmp+b1);
}

}
}

