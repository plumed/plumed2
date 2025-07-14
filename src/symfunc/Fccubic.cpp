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
#include "function/FunctionSetup.h"
#include "function/FunctionShortcut.h"
#include "function/FunctionOfMatrix.h"
#include "core/ActionRegister.h"

#include <string>
#include <cmath>

namespace PLMD {
namespace symfunc {

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
good computational performance.  Notice also that you can you can rotate the bond vectors before computing the
function in the above expression by using the PHI, THETA and PSI keywords as is discussed in the documentation for [COORDINATION_SHELL_FUNCTION](COORDINATION_SHELL_FUNCTION.md).

The following input calculates the FCCUBIC parameter for the 64 atoms in the system
and then calculates and prints the average value for this quantity.

```plumed
d: FCCUBIC SPECIES=1-64 D_0=3.0 R_0=1.5 NN=6 MM=12
dm: MEAN ARG=d PERIODIC=NO
PRINT ARG=dm FILE=colv
```

In the input above we use a rational [switching function](LESS_THAN.md) with the parameters above. We would recommend using SWITCH syntax
rather than the syntax above when giving the parameters for the switching function as you can then use any of the switching functions described
in the documentation for [LESS_THAN](LESS_THAN.md).  More importantly, however, using this syntax allows you to set the D_MAX parameter for the
switching function as demonstrated below:

```plumed
d: FCCUBIC SPECIES=1-64 SWITCH={RATIONAL D_0=3.0 R_0=1.5 D_MAX=6.0}
dm: MEAN ARG=d PERIODIC=NO
PRINT ARG=dm FILE=colv
```

Setting the `D_MAX` can substantially improve PLUMED performance as it turns on the linked list algorithm that is discussed in the optimisation details part
of the documentation for [CONTACT_MATRIX](CONTACT_MATRIX.md).

## Working with two types of atom

If you would like to calculate whether the atoms in GROUPB are arranged around the atoms in GROUPA as they in in an FCC structure you use an input like the one
shown below:

```plumed
d: FCCUBIC SPECIESA=1-64 SPECIESB=65-200 SWITCH={RATIONAL D_0=3.0 R_0=1.5 D_MAX=6.0}
lt: MORE_THAN ARG=d SWITCH={RATIONAL R_0=0.5}
s: SUM ARG=lt PERIODIC=NO
PRINT ARG=s FILE=colv
```

This input calculates how many of the 64 atoms that were input to the SPECIESA keyword have an fcc value that is greater than 0.5.

## The MASK keyword

You can use the MASK keyword with this action in the same way that it is used with [COORDINATIONNUMBER](COORDINATIONNUMBER.md).  This keyword thus expects a vector in
input, which tells FCCUBIC that it is safe not to calculate the FCCUBIC parameter for some of the atoms.  As illustrated below, this is useful if you are using functionality
from the [volumes module](module_volumes.md) to calculate the average value of the FCCUBIC parameter for only those atoms that lie in a certain part of the simulation box.

```plumed
# Fixed virtual atom which serves as the probe volume's center (pos. in nm)
center: FIXEDATOM AT=2.5,2.5,2.5
# Vector in which element i is one if atom i is in sphere of interest and zero otherwise
sphere: INSPHERE ATOMS=1-400 CENTER=center RADIUS={GAUSSIAN D_0=0.5 R_0=0.01 D_MAX=0.52}
# Calculate the fccubic parameter of the atoms
cc: FCCUBIC ...
  SPECIES=1-400 MASK=sphere
  SWITCH={RATIONAL D_0=3.0 R_0=1.5 D_MAX=6.0}
...
# Multiply fccubic parameters numbers by sphere vector
prod: CUSTOM ARG=cc,sphere FUNC=x*y PERIODIC=NO
# Sum of coordination numbers for atoms that are in the sphere of interest
numer: SUM ARG=prod PERIODIC=NO
# Number of atoms that are in sphere of interest
denom: SUM ARG=sphere PERIODIC=NO
# Average coordination number for atoms in sphere of interest
av: CUSTOM ARG=numer,denom FUNC=x/y PERIODIC=NO
# And print out final CV to a file
PRINT ARG=av FILE=colvar STRIDE=1
```

This input calculate the average value of the FCCUBIC parameter for only those atoms that are within a spherical region that is centered on the point
$(2.5,2.5,2.5)$.

## Deprecated syntax

More information on the deprecated keywords that are given below is available in the documentation for the [DISTANCES](DISTANCES.md) command.

*/
//+ENDPLUMEDOC

//+PLUMEDOC MCOLVAR FCCUBIC_FUNC
/*
Measure how similar the environment around atoms is to that found in a FCC structure.

This is the function that is used in the [FCCUBIC](FCCUBIC.md) shortcut. You can see an example that shows how it is used
if you expand the FCCUBIC shortcut in the following example input:

```plumed
d: FCCUBIC SPECIES=1-64 SWITCH={RATIONAL D_0=3.0 R_0=1.5} MEAN
PRINT ARG=d.* FILE=colv
```

Notice that the decomposition of the function that is illustrated above allows you to do things like the calculation in the following input:

```plumed
# Calculate the distances between atoms
d_mat: DISTANCE_MATRIX GROUP=1-64 CUTOFF=4.5 COMPONENTS
# Find the six nearest atom to each of the coordinates
nn: NEIGHBORS ARG=d_mat.w NLOWEST=6
# Evalulate the FCC function
d_vfunc: FCCUBIC_FUNC ARG=d_mat.x,d_mat.y,d_mat.z MASK=nn ALPHA=3.0
# Take the product of the weights with the fcc function evaluations
d_wvfunc: CUSTOM ARG=d_vfunc,nn FUNC=x*y PERIODIC=NO
# Calculate the sum of fcc cubic function values for each atom
d_ones: ONES SIZE=64
d: MATRIX_VECTOR_PRODUCT ARG=d_wvfunc,d_ones
# Calculate the number of neighbours
d_denom: MATRIX_VECTOR_PRODUCT ARG=nn,d_ones
# Calculate the average value of the fcc cubic function per bonds
d_n: CUSTOM ARG=d,d_denom FUNC=x/y PERIODIC=NO
d_mean: MEAN ARG=d_n PERIODIC=NO
PRINT ARG=d_mean FILE=colv
```

In this input we evaluate the [FCCUBIC](FCCUBIC.md) function for the six nearest neighbours to each atom rather than the atoms that are within a certain cutoff.
By using the MASK keyword in the input to FCCUBIC_FUNC we ensure that the [FCCUBIC](FCCUBIC.md) function is only evaluated for the six nearest neighbors.

*/
//+ENDPLUMEDOC

class Fccubic {
public:
  double alpha, a1, b1;
  static void registerKeywords( Keywords& keys );
  static void read( Fccubic& func, ActionWithArguments* action, function::FunctionOptions& funcout );
  static void calc( const Fccubic& func, bool noderiv, View<const double> args, function::FunctionOutput& funcout );
  Fccubic& operator=(const Fccubic& m) {
    alpha = m.alpha;
    a1 = m.a1;
    b1 = m.b1;
    return *this;
  }
};

typedef function::FunctionShortcut<Fccubic> FccubicShortcut;
PLUMED_REGISTER_ACTION(FccubicShortcut,"FCCUBIC_FUNC")
typedef function::FunctionOfMatrix<Fccubic> MatrixFccubic;
PLUMED_REGISTER_ACTION(MatrixFccubic,"FCCUBIC_FUNC_MATRIX")

void Fccubic::registerKeywords( Keywords& keys ) {
  keys.add("compulsory","ALPHA","3.0","The alpha parameter of the angular function");
  keys.setValueDescription("matrix","a function that measures the similarity with an fcc environment");
}

void Fccubic::read( Fccubic& func, ActionWithArguments* action, function::FunctionOptions& funcout ) {
  // Scaling factors such that '1' corresponds to fcc lattice
  // and '0' corresponds to isotropic (liquid)
  action->parse("ALPHA",func.alpha);
  func.a1 = 80080. / (2717. + 16*func.alpha);
  func.b1 = 16.*(func.alpha-143)/(2717+16*func.alpha);
  action->log.printf("  setting alpha paramter equal to %f \n",func.alpha);
}

void Fccubic::calc( const Fccubic& func, bool noderiv, const View<const double> args, function::FunctionOutput& funcout ) {
  double x2 = args[0]*args[0];
  double x4 = x2*x2;

  double y2 = args[1]*args[1];
  double y4 = y2*y2;

  double z2 = args[2]*args[2];
  double z4 = z2*z2;

  double d2 = x2 + y2 + z2;
  if( d2 < epsilon ) {
    d2 = 1;
  }
  double r8 = pow( d2, 4 );
  double r12 = pow( d2, 6 );

  double tmp = ((x4*y4)+(x4*z4)+(y4*z4))/r8-func.alpha*x4*y4*z4/r12;

  double t0 = (x2*y4+x2*z4)/r8-func.alpha*x2*y4*z4/r12;
  double t1 = (y2*x4+y2*z4)/r8-func.alpha*y2*x4*z4/r12;
  double t2 = (z2*x4+z2*y4)/r8-func.alpha*z2*x4*y4/r12;
  double t3 = (2*tmp-func.alpha*x4*y4*z4/r12)/d2;

  if( !noderiv ) {
    funcout.derivs[0][0]=4*func.a1*args[0]*(t0-t3);
    funcout.derivs[0][1]=4*func.a1*args[1]*(t1-t3);
    funcout.derivs[0][2]=4*func.a1*args[2]*(t2-t3);
  }

  // Set the value and the derivatives
  funcout.values[0] = (func.a1*tmp+func.b1);
}

}
}

