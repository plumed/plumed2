/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2018 The plumed team
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
#include "core/ActionRegister.h"
#include "core/ActionShortcut.h"
#include "core/ActionWithValue.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"

//+PLUMEDOC FUNCTION NORMALIZED_EUCLIDEAN_DISTANCE
/*
Calculate the normalised euclidean distance between two points in CV space

If we have two $n$-dimensional vectors $u$ and $v$ and an $n$ dimensional vector of inverse covariance values, $a$,
we can calculate the normalised Euclidean distance between the two points as

$$
d = \sqrt{ \sum_{i=1}^n a_i (u_i - v_i)^2 }
$$

which can be expressed in matrix form as:

$$
d^2 = (u-v)^T a \odot (u-v)
$$

where $\odot$ here is used to indicate the [Hadamard product](https://en.wikipedia.org/wiki/Hadamard_product_(matrices))
of the two vectors.

The inputs below shows an example where this is used to calculate the Normalized Euclidean distance
between the instaneous values of some torsional angles and some reference values
for these torsion.  The inverse covriance values are provided in the constant value with label `m`.
In this first example the input values are vectors:

```plumed
m: CONSTANT VALUES=0.1,0.2,0.3
c: CONSTANT VALUES=1,2,3
d: DISTANCE ATOMS1=1,2 ATOMS2=3,4 ATOMS3=5,6
dd: NORMALIZED_EUCLIDEAN_DISTANCE ARG1=c ARG2=d METRIC=m
PRINT ARG=dd FILE=colvar
```

while this second example does the same thing but uses scalars in input.

```plumed
m: CONSTANT VALUES=0.1,0.2,0.3
c1: CONSTANT VALUE=1
d1: DISTANCE ATOMS=1,2
c2: CONSTANT VALUE=2
d2: DISTANCE ATOMS=3,4
c3: CONSTANT VALUE=3
d3: DISTANCE ATOMS=5,6
dd: NORMALIZED_EUCLIDEAN_DISTANCE ARG1=c1,c2,c3 ARG2=d1,d2,d3 METRIC=m
PRINT ARG=dd FILE=colvar
```

Lastly, note that if you want to calculate the square of the distance rather than the distance you can use
the `SQUARED` flag as shown below:

```plumed
m: CONSTANT VALUES=0.1,0.2,0.3
c: CONSTANT VALUES=1,2,3
d: DISTANCE ATOMS1=1,2 ATOMS2=3,4 ATOMS3=5,6
dd: NORMALIZED_EUCLIDEAN_DISTANCE ARG1=c ARG2=d METRIC=m SQUARED
PRINT ARG=dd FILE=colvar
```

Calculating the square of the distance is slightly cheapter than computing the distance as you avoid taking the square root.

## Calculating multiple distances

Suppose that we now have $m$ reference configurations we can define the following $m$ distances
from these reference configurations:

$$
d_j^2 = (u-v_j)^T a \odot (u-v_j)
$$

Lets suppose that we put the $m$, $n$-dimensional $(u-v_j)$ vectors in this expression into a
$n\times m$ matrix, $A$, by using the [DISPLACEMENT](DISPLACEMENT.md) command.  It is then
straightforward to show that the $d_j^2$ values in the above expression are the diagonal
elements of the matrix product $A^T K \cdot A$, where $K$ is an $n \times m$ matrix that contains
$m$ copies of the inverse covariance matrix $a$ in its columns.

We can use this idea to calculate multiple NORMALIZED_EUCLIDEAN_DISTANCE values in the following inputs.
This first example calculates the three distances between the instaneoues values of two torsions
and three reference configurations.

```plumed
m: CONSTANT VALUES=0.1,0.2
ref_psi: CONSTANT VALUES=2.25,1.3,-1.5
ref_phi: CONSTANT VALUES=-1.91,-0.6,2.4

psi: TORSION ATOMS=1,2,3,4
phi: TORSION ATOMS=13,14,15,16

dd: NORMALIZED_EUCLIDEAN_DISTANCE ARG2=psi,phi ARG1=ref_psi,ref_phi METRIC=m
PRINT ARG=dd FILE=colvar
```

This second example calculates the three distances between a single reference value for the two
torsions and three instances of this pair of torsions.

```plumed
m: CONSTANT VALUES=0.1,0.2
ref_psi: CONSTANT VALUES=2.25
ref_phi: CONSTANT VALUES=-1.91

psi: TORSION ATOMS1=1,2,3,4 ATOMS2=5,6,7,8 ATOMS3=9,10,11,12
phi: TORSION ATOMS1=13,14,15,16 ATOMS2=17,18,19,20 ATOMS3=21,22,23,24

dd: NORMALIZED_EUCLIDEAN_DISTANCE ARG1=psi,phi ARG2=ref_psi,ref_phi METRIC=m
PRINT ARG=dd FILE=colvar
```

This final example then computes three distances between three pairs of torsional angles and threee
reference values for these three values.

```plumed
m: CONSTANT VALUES=0.1,0.2
ref_psi: CONSTANT VALUES=2.25,1.3,-1.5
ref_phi: CONSTANT VALUES=-1.91,-0.6,2.4

psi: TORSION ATOMS1=1,2,3,4 ATOMS2=5,6,7,8 ATOMS3=9,10,11,12
phi: TORSION ATOMS1=13,14,15,16 ATOMS2=17,18,19,20 ATOMS3=21,22,23,24

dd: NORMALIZED_EUCLIDEAN_DISTANCE ARG1=psi,phi ARG2=ref_psi,ref_phi METRIC=m
PRINT ARG=dd FILE=colvar
```

!!! note "scalars must be specified in ARG2"

    If you use a mixture of vectors are scalars when specifying the input to to this action the
    vectors should be passed using the ARG1 keyword and the scalars must be passed in the ARG2 keyword
    as is done in the example inputs above.
*/
//+ENDPLUMEDOC

namespace PLMD {
namespace refdist {

class NormalizedEuclideanDistance : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  explicit NormalizedEuclideanDistance(const ActionOptions&ao);
};

PLUMED_REGISTER_ACTION(NormalizedEuclideanDistance,"NORMALIZED_EUCLIDEAN_DISTANCE")

void NormalizedEuclideanDistance::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords(keys);
  keys.add("compulsory","ARG1","The poin that we are calculating the distance from");
  keys.add("compulsory","ARG2","The point that we are calculating the distance to");
  keys.add("compulsory","METRIC","The inverse covariance matrix that should be used when calculating the distance");
  keys.addFlag("SQUARED",false,"The squared distance should be calculated");
  keys.setValueDescription("scalar/vector","the normalized euclidean distances between the input vectors");
  keys.needsAction("DISPLACEMENT");
  keys.needsAction("CUSTOM");
  keys.needsAction("OUTER_PRODUCT");
  keys.needsAction("TRANSPOSE");
  keys.needsAction("MATRIX_PRODUCT_DIAGONAL");
  keys.needsAction("ONES");
}

NormalizedEuclideanDistance::NormalizedEuclideanDistance( const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  std::string arg1, arg2, metstr;
  parse("ARG1",arg1);
  parse("ARG2",arg2);
  parse("METRIC",metstr);
  // Vectors are in rows here
  readInputLine( getShortcutLabel() + "_diff: DISPLACEMENT ARG1=" + arg1 + " ARG2=" + arg2 );
  // Vectors are in columns here
  readInputLine( getShortcutLabel() + "_diffT: TRANSPOSE ARG=" + getShortcutLabel() + "_diff");
  // Get the action that computes the differences
  ActionWithValue* av = plumed.getActionSet().selectWithLabel<ActionWithValue*>( getShortcutLabel() + "_diffT");
  plumed_assert( av );
  // If this is a matrix we need create a matrix to multiply by
  if( av->copyOutput(0)->getRank()==2 ) {
    // Create some ones
    std::string nones;
    Tools::convert( av->copyOutput(0)->getShape()[1], nones );
    readInputLine( getShortcutLabel() + "_ones: ONES SIZE=" + nones);
    // Now do some multiplication to create a matrix that can be multiplied by our "inverse variance" vector
    if( av->copyOutput(0)->getShape()[0]==1 ) {
      readInputLine( getShortcutLabel() + "_" + metstr + "T: CUSTOM ARG=" + getShortcutLabel() + "_ones," + metstr + " FUNC=x*y PERIODIC=NO");
      readInputLine( getShortcutLabel() + "_" + metstr + ": TRANSPOSE ARG=" + getShortcutLabel() + "_" + metstr + "T");
    } else {
      readInputLine( getShortcutLabel() + "_" + metstr + ": OUTER_PRODUCT ARG=" + metstr + "," + getShortcutLabel() + "_ones");
    }
    metstr = getShortcutLabel() + "_" + metstr;
  }
  // Now do the multiplication
  readInputLine( getShortcutLabel() + "_sdiff: CUSTOM ARG=" + getShortcutLabel() +"_diffT," + metstr + " FUNC=x*y PERIODIC=NO");
  bool squared;
  parseFlag("SQUARED",squared);
  std::string olab = getShortcutLabel();
  if( !squared ) {
    olab += "_2";
  }
  readInputLine( olab + ": MATRIX_PRODUCT_DIAGONAL ARG=" + getShortcutLabel() +"_diff," + getShortcutLabel() + "_sdiff");
  if( !squared ) {
    readInputLine( getShortcutLabel() + ": CUSTOM ARG=" + getShortcutLabel() + "_2 FUNC=sqrt(x) PERIODIC=NO");
  }
}

}
}
