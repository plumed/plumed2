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
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/ActionShortcut.h"
#include "core/ActionWithValue.h"

//+PLUMEDOC FUNCTION MAHALANOBIS_DISTANCE
/*
Calculate the Mahalanobis distance between two points in CV space

If we have two $n$-dimensional vectors $u$ and $v$ we can calculate the
[Mahalanobis distance](https://en.wikipedia.org/wiki/Mahalanobis_distance) between the two points as

$$
d = \sqrt{ \sum_{i=1}^n \sum_{j=1}^n m_{ij} (u_i - v_i)(u_j - v_j) }
$$

which can be expressed in matrix form as:

$$
d^2 = (u-v)^T M (u-v)
$$

The inputs below shows an example where this is used to calculate the Mahalanobis distance
between the instaneous values of some torsional angles and some reference values
for these distances.  The inverse covriance values are provided in the constant value with label `m`.
In this first example the input values are vectors:

```plumed
m: CONSTANT VALUES=2.45960237E-0001,-1.30615381E-0001,-1.30615381E-0001,2.40239117E-0001 NROWS=2 NCOLS=2
c: CONSTANT VALUES=1,2
d: DISTANCE ATOMS1=1,2 ATOMS2=3,4
dd: MAHALANOBIS_DISTANCE ARG1=c ARG2=d METRIC=m
PRINT ARG=dd FILE=colvar
```

while this second example does the same thing but uses scalars in input.

```plumed
m: CONSTANT VALUES=2.45960237E-0001,-1.30615381E-0001,-1.30615381E-0001,2.40239117E-0001 NROWS=2 NCOLS=2
c1: CONSTANT VALUE=1
d1: DISTANCE ATOMS=1,2
c2: CONSTANT VALUE=2
d2: DISTANCE ATOMS=3,4
dd: MAHALANOBIS_DISTANCE ARG1=c1,c2 ARG2=d1,d2 METRIC=m
PRINT ARG=dd FILE=colvar
```

Lastly, note that if you want to calculate the square of the distance rather than the distance you can use
the `SQUARED` flag as shown below:

```plumed
m: CONSTANT VALUES=2.45960237E-0001,-1.30615381E-0001,-1.30615381E-0001,2.40239117E-0001 NROWS=2 NCOLS=2
c: CONSTANT VALUES=1,2
d: DISTANCE ATOMS1=1,2 ATOMS2=3,4
dd: MAHALANOBIS_DISTANCE ARG1=c ARG2=d METRIC=m SQUARED
PRINT ARG=dd FILE=colvar
```

Calculating the square of the distance is slightly cheapter than computing the distance as you avoid taking the square root.

## Dealing with periodic variables

When you are calculating a distance from a reference point you need to be careful when the input variables
are periodic. If you are calculating the distance using the [EUCLIDEAN_DISTANCE](EUCLIDEAN_DISTANCE.md) and
[NORMALIZED_EUCLIDEAN_DISTANCE](NORMALIZED_EUCLIDEAN_DISTANCE.md) commands this is not a problem. The problems are
specific to the Mahalanobis distance command and have been resolved in the papers that are cited below by defining
the following alternatative to the Mahalanobis distance:

$$
d^2 = 2\sum_{i=1}^n m_{ii} \left[ 1 - \cos\left( \frac{2\pi(u_i-v_i)}{P_i} \right) \right] + \sum_{i\ne j} m_{ij} \sin\left( \frac{2\pi(u_i-v_i)}{P_i} \right) \sin\left( \frac{2\pi(u_j-v_j)}{P_j} \right)
$$

In this expression, $P_i$ indicates the periodicity of the domain for variable $i$. If you would like to compute this
distance with PLUMED you use the `VON_MISSES` shown below:

```plumed
m: CONSTANT VALUES=2.45960237E-0001,-1.30615381E-0001,-1.30615381E-0001,2.40239117E-0001 NROWS=2 NCOLS=2
c: CONSTANT VALUES=1,2
d: TORSION ATOMS1=1,2,3,4 ATOMS2=5,6,7,8
dd: MAHALANOBIS_DISTANCE ARG1=c ARG2=d METRIC=m VON_MISSES
PRINT ARG=dd FILE=colvar
```

## Calculating multiple distances

Suppose that we now have $m$ reference configurations we can define the following $m$ distances
from these reference configurations:

$$
d_j^2 = (u-v_j)^T M (u-v_j)
$$

Lets suppose that we put the $m$, $n$-dimensional $(u-v_j)$ vectors in this expression into a
$n\times m$ matrix, $A$, by using the [DISPLACEMENT](DISPLACEMENT.md) command.  It is then
straightforward to show that the $d_j^2$ values in the above expression are the diagonal
elements of the matrix product $A^T M A$.

We can use this idea to calculate multiple MAHALANOBIS_DISTANCE values in the following inputs.
This first example calculates the three distances between the instaneoues values of two torsions
and three reference configurations.

```plumed
m: CONSTANT VALUES=2.45960237E-0001,-1.30615381E-0001,-1.30615381E-0001,2.40239117E-0001 NROWS=2 NCOLS=2
ref_psi: CONSTANT VALUES=2.25,1.3,-1.5
ref_phi: CONSTANT VALUES=-1.91,-0.6,2.4

psi: TORSION ATOMS=1,2,3,4
phi: TORSION ATOMS=13,14,15,16

dd: MAHALANOBIS_DISTANCE ARG2=psi,phi ARG1=ref_psi,ref_phi METRIC=m
PRINT ARG=dd FILE=colvar
```

This section example calculates the three distances between a single reference value for the two
torsions and three instances of this pair of torsions.

```plumed
m: CONSTANT VALUES=2.45960237E-0001,-1.30615381E-0001,-1.30615381E-0001,2.40239117E-0001 NROWS=2 NCOLS=2
ref_psi: CONSTANT VALUES=2.25
ref_phi: CONSTANT VALUES=-1.91

psi: TORSION ATOMS1=1,2,3,4 ATOMS2=5,6,7,8 ATOMS3=9,10,11,12
phi: TORSION ATOMS1=13,14,15,16 ATOMS2=17,18,19,20 ATOMS3=21,22,23,24

dd: MAHALANOBIS_DISTANCE ARG1=psi,phi ARG2=ref_psi,ref_phi METRIC=m
PRINT ARG=dd FILE=colvar
```

This final example then computes three distances between three pairs of torsional angles and threee
reference values for these three values.

```plumed
m: CONSTANT VALUES=2.45960237E-0001,-1.30615381E-0001,-1.30615381E-0001,2.40239117E-0001 NROWS=2 NCOLS=2
ref_psi: CONSTANT VALUES=2.25,1.3,-1.5
ref_phi: CONSTANT VALUES=-1.91,-0.6,2.4

psi: TORSION ATOMS1=1,2,3,4 ATOMS2=5,6,7,8 ATOMS3=9,10,11,12
phi: TORSION ATOMS1=13,14,15,16 ATOMS2=17,18,19,20 ATOMS3=21,22,23,24

dd: MAHALANOBIS_DISTANCE ARG1=psi,phi ARG2=ref_psi,ref_phi METRIC=m
PRINT ARG=dd FILE=colvar
```

Notice, finally, that you can also calculate multiple distances if you use the `VON_MISSES` option:

```plumed
m: CONSTANT VALUES=2.45960237E-0001,-1.30615381E-0001,-1.30615381E-0001,2.40239117E-0001 NROWS=2 NCOLS=2
ref_psi: CONSTANT VALUES=2.25,1.3,-1.5
ref_phi: CONSTANT VALUES=-1.91,-0.6,2.4

psi: TORSION ATOMS1=1,2,3,4 ATOMS2=5,6,7,8 ATOMS3=9,10,11,12
phi: TORSION ATOMS1=13,14,15,16 ATOMS2=17,18,19,20 ATOMS3=21,22,23,24

dd: MAHALANOBIS_DISTANCE ARG1=psi,phi ARG2=ref_psi,ref_phi METRIC=m VON_MISSES
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

class MahalanobisDistance : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  explicit MahalanobisDistance(const ActionOptions&ao);
};

PLUMED_REGISTER_ACTION(MahalanobisDistance,"MAHALANOBIS_DISTANCE")

void MahalanobisDistance::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords(keys);
  keys.add("compulsory","ARG1","The point that we are calculating the distance from");
  keys.add("compulsory","ARG2","The point that we are calculating the distance to");
  keys.add("compulsory","METRIC","The inverse covariance matrix that should be used when calculating the distance");
  keys.addFlag("SQUARED",false,"The squared distance should be calculated");
  keys.addFlag("VON_MISSES",false,"Compute the mahalanobis distance in a way that is more sympathetic to the periodic boundary conditions");
  keys.setValueDescription("scalar/vector","the Mahalanobis distances between the input vectors");
  keys.needsAction("DISPLACEMENT");
  keys.needsAction("CUSTOM");
  keys.needsAction("OUTER_PRODUCT");
  keys.needsAction("TRANSPOSE");
  keys.needsAction("MATRIX_PRODUCT_DIAGONAL");
  keys.needsAction("CONSTANT");
  keys.needsAction("MATRIX_VECTOR_PRODUCT");
  keys.needsAction("MATRIX_PRODUCT");
  keys.needsAction("COMBINE");
  keys.addDOI("10.1073/pnas.1011511107");
  keys.addDOI("10.1021/acs.jctc.7b00993");
}

MahalanobisDistance::MahalanobisDistance( const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  std::string arg1, arg2, metstr;
  parse("ARG1",arg1);
  parse("ARG2",arg2);
  parse("METRIC",metstr);
  // Check on input metric
  ActionWithValue* mav=plumed.getActionSet().selectWithLabel<ActionWithValue*>( metstr );
  if( !mav ) {
    error("could not find action named " + metstr + " to use for metric");
  }
  if( mav->copyOutput(0)->getRank()!=2 ) {
    error("metric has incorrect rank");
  }

  readInputLine( getShortcutLabel() + "_diff: DISPLACEMENT ARG1=" + arg1 + " ARG2=" + arg2 );
  readInputLine( getShortcutLabel() + "_diffT: TRANSPOSE ARG=" + getShortcutLabel() + "_diff");
  bool von_miss, squared;
  parseFlag("VON_MISSES",von_miss);
  parseFlag("SQUARED",squared);
  if( von_miss ) {
    unsigned nrows = mav->copyOutput(0)->getShape()[0];
    if( mav->copyOutput(0)->getShape()[1]!=nrows ) {
      error("metric is not symmetric");
    }
    // Create a matrix that can be used to compute the off diagonal elements
    std::string valstr, nrstr;
    Tools::convert( mav->copyOutput(0)->get(0), valstr );
    Tools::convert( nrows, nrstr );
    std::string diagmet = getShortcutLabel() + "_diagmet: CONSTANT VALUES=" + valstr;
    std::string offdiagmet = getShortcutLabel() + "_offdiagmet: CONSTANT NROWS=" + nrstr + " NCOLS=" + nrstr + " VALUES=0";
    for(unsigned i=0; i<nrows; ++i) {
      for(unsigned j=0; j<nrows; ++j) {
        Tools::convert( mav->copyOutput(0)->get(i*nrows+j), valstr );
        if( i==j && i>0 ) {
          offdiagmet += ",0";
          diagmet += "," + valstr;
        } else if( i!=j ) {
          offdiagmet += "," + valstr;
        }
      }
    }
    readInputLine( diagmet );
    readInputLine( offdiagmet );
    // Compute distances scaled by periods
    ActionWithValue* av=plumed.getActionSet().selectWithLabel<ActionWithValue*>( getShortcutLabel() + "_diff" );
    plumed_assert( av );
    if( !av->copyOutput(0)->isPeriodic() ) {
      error("VON_MISSES only works with periodic variables");
    }
    std::string min, max;
    av->copyOutput(0)->getDomain(min,max);
    readInputLine( getShortcutLabel() + "_scaled: CUSTOM ARG=" + getShortcutLabel() + "_diffT FUNC=2*pi*x/(" + max +"-" + min + ") PERIODIC=NO");
    // We start calculating off-diagonal elements by computing the sines of the scaled differences (this is a column vector)
    readInputLine( getShortcutLabel() + "_sinediffT: CUSTOM ARG=" + getShortcutLabel() + "_scaled FUNC=sin(x) PERIODIC=NO");
    // Transpose sines to get a row vector
    readInputLine( getShortcutLabel() + "_sinediff: TRANSPOSE ARG=" + getShortcutLabel() + "_sinediffT");
    // Compute the off diagonal elements
    ActionWithValue* avs=plumed.getActionSet().selectWithLabel<ActionWithValue*>( getShortcutLabel() + "_sinediffT" );
    plumed_assert( avs && avs->getNumberOfComponents()==1 );
    if( (avs->copyOutput(0))->getRank()==1 ) {
      readInputLine( getShortcutLabel() + "_matvec: MATRIX_VECTOR_PRODUCT ARG=" + metstr + "," + getShortcutLabel() +"_sinediffT");
    } else {
      readInputLine( getShortcutLabel() + "_matvec: MATRIX_PRODUCT ARG=" + getShortcutLabel() + "_offdiagmet," + getShortcutLabel() +"_sinediffT");
    }
    readInputLine( getShortcutLabel() + "_offdiag: MATRIX_PRODUCT_DIAGONAL ARG=" + getShortcutLabel() + "_sinediff," + getShortcutLabel() +"_matvec");
    // Sort out the metric for the diagonal elements
    std::string metstr2 = getShortcutLabel() + "_diagmet";
    // If this is a matrix we need create a matrix to multiply by
    if( av->copyOutput(0)->getShape()[0]>1 ) {
      // Create some ones
      std::string ones=" VALUES=1";
      for(unsigned i=1; i<av->copyOutput(0)->getShape()[0]; ++i ) {
        ones += ",1";
      }
      readInputLine( getShortcutLabel() + "_ones: CONSTANT " + ones );
      // Now do some multiplication to create a matrix that can be multiplied by our "inverse variance" vector
      readInputLine( getShortcutLabel() + "_" + metstr + ": OUTER_PRODUCT ARG=" + metstr2 + "," + getShortcutLabel() + "_ones");
      metstr2 = getShortcutLabel() + "_" + metstr;
    }
    // Compute the diagonal elements
    readInputLine( getShortcutLabel() + "_prod: CUSTOM ARG=" + getShortcutLabel() + "_scaled," + metstr2 + " FUNC=2*(1-cos(x))*y PERIODIC=NO");
    std::string ncstr;
    Tools::convert( nrows, ncstr );
    Tools::convert( av->copyOutput(0)->getShape()[0], nrstr );
    std::string ones=" VALUES=1";
    for(unsigned i=1; i<av->copyOutput(0)->getNumberOfValues(); ++i) {
      ones += ",1";
    }
    readInputLine( getShortcutLabel() + "_matones: CONSTANT NROWS=" + nrstr + " NCOLS=" + ncstr + ones );
    readInputLine( getShortcutLabel() + "_diag: MATRIX_PRODUCT_DIAGONAL ARG=" + getShortcutLabel() + "_matones," + getShortcutLabel() + "_prod");
    // Sum everything
    if( !squared ) {
      readInputLine( getShortcutLabel() + "_2: COMBINE ARG=" + getShortcutLabel() + "_offdiag," + getShortcutLabel() + "_diag PERIODIC=NO");
    } else {
      readInputLine( getShortcutLabel() + ": COMBINE ARG=" + getShortcutLabel() + "_offdiag," + getShortcutLabel() + "_diag PERIODIC=NO");
    }
  } else {
    ActionWithValue* av=plumed.getActionSet().selectWithLabel<ActionWithValue*>( getShortcutLabel() + "_diffT" );
    plumed_assert( av && av->getNumberOfComponents()==1 );
    if( (av->copyOutput(0))->getRank()==1 ) {
      readInputLine( getShortcutLabel() + "_matvec: MATRIX_VECTOR_PRODUCT ARG=" + metstr + "," + getShortcutLabel() +"_diffT");
    } else {
      readInputLine( getShortcutLabel() + "_matvec: MATRIX_PRODUCT ARG=" + metstr + "," + getShortcutLabel() +"_diffT");
    }
    std::string olab = getShortcutLabel();
    if( !squared ) {
      olab += "_2";
    }
    readInputLine( olab + ": MATRIX_PRODUCT_DIAGONAL ARG=" + getShortcutLabel() + "_diff," + getShortcutLabel() +"_matvec");
  }
  if( !squared ) {
    readInputLine( getShortcutLabel() + ": CUSTOM ARG=" + getShortcutLabel() + "_2 FUNC=sqrt(x) PERIODIC=NO");
  }
}

}
}
