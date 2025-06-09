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

//+PLUMEDOC MCOLVAR DISPLACEMENT
/*
Calculate the displacement vector between the pair of input vectors

This shortcut can be used to calculate the vector of displacements between two input vectors as shown below.

```plumed
c: CONSTANT VALUES=1,2,3
d: DISTANCE ATOMS1=1,2 ATOMS2=3,4 ATOMS3=5,6
e: DISPLACEMENT ARG1=c ARG2=d
PRINT ARG=e FILE=colvar
```

The output here, `e`, is a $1 \times 3$ matrix for reasons that will become clear later in this documentation that is computed as:

$$
\mathbf{e} = \mathbf{c} - \mathbf{d}
$$

Notice that we can obtain the same result by specifying the input vectors here as two sets of three scalars as shown
below:

```plumed
c1: CONSTANT VALUE=1
d1: DISTANCE ATOMS=1,2
c2: CONSTANT VALUE=2
d2: DISTANCE ATOMS=3,4
c3: CONSTANT VALUE=3
d3: DISTANCE ATOMS=5,6
e: DISPLACEMENT ARG1=c1,c2,c3 ARG2=d1,d2,d3
PRINT ARG=e FILE=colvar
```

The DISPLACEMENT command that has been introduced in the above inputs is primarily used within the [EUCLIDEAN_DISTANCE](EUCLIDEAN_DISTANCE.md),
[NORMALIZED_EUCLIDEAN_DISTANCE](NORMALIZED_EUCLIDEAN_DISTANCE.md) and [MAHALANOBIS_DISTANCE](MAHALANOBIS_DISTANCE.md) shortcuts.  If the $1 \times N$ matrix
of displacements that that we obtainfrom these commands is, $E$, these three actions calculate

$$
d = E M E^T
$$

The $N \times N$ matrix $M$ here is the identity if you are using [EUCLIDEAN_DISTANCE](EUCLIDEAN_DISTANCE.md), a diagonal matrix if you are using
[NORMALIZED_EUCLIDEAN_DISTANCE](NORMALIZED_EUCLIDEAN_DISTANCE.md) and a full matrix if you are computing the [MAHALANOBIS_DISTANCE](MAHALANOBIS_DISTANCE.md).

## Calculating multiple displacement vectors

The reason the output of DISPLACEMENT is a $1 \times 3$ matrix here becomes clearer once we consider the following input:

```plumed
ref_psi: CONSTANT VALUES=2.25
ref_phi: CONSTANT VALUES=-1.91

psi: TORSION ATOMS1=1,2,3,4 ATOMS2=5,6,7,8 ATOMS3=9,10,11,12
phi: TORSION ATOMS1=13,14,15,16 ATOMS2=17,18,19,20 ATOMS3=21,22,23,24

dd: DISPLACEMENT ARG1=psi,phi ARG2=ref_psi,ref_phi
PRINT ARG=dd FILE=colvar
```

The output from the input above is a $3\times 2$ matrix.  The rows of this matrix run over the 3 different
torsion values that have been specified in the `psi` and `phi` commands.  The first column of the matrix
contains the differences between each of the instantaneous `psi` aingles and the reference value for this
angle, while the second columns contains the differences between the `phi` angles and the reference.

In other words, we can calculate multiple displacement vectors at once as each row of the final output matrix will
contain a vector of displacements between two vectors.  Notice that we can use a similar input to calculate the
differences between the instantaneous value of a pair of torsions and 3 reference values as shown below:

```plumed
ref_psi: CONSTANT VALUES=2.25,1.3,-1.5
ref_phi: CONSTANT VALUES=-1.91,-0.6,2.4

psi: TORSION ATOMS=1,2,3,4
phi: TORSION ATOMS=13,14,15,16

dd: DISPLACEMENT ARG2=psi,phi ARG1=ref_psi,ref_phi
PRINT ARG=dd FILE=colvar
```

!!! note "scalars must be specified in ARG2"

    If you use a mixture of vectors are scalars when specifying the input to to this action the
    vectors should be passed using the ARG1 keyword and the scalars must be passed in the ARG2 keyword
    as is done in the example inputs above. Obviously, this limitation sometimes means that you have to add in
    an additional [CUSTOM](CUSTOM.md) action that multiplies the output vectors by -1 to get the vectors pointing
    in the direction you desire.

The output here will again be a $3\times 2$ matrix with each of the three rows holding a vector of displacements
between the 2 instananeous values and one of the three sets of reference values.

Lastly, we can use two sets of vectors in the input to DISPLACEMENT as shown below:

```plumed
ref_psi: CONSTANT VALUES=2.25,1.3,-1.5
ref_phi: CONSTANT VALUES=-1.91,-0.6,2.4

psi: TORSION ATOMS1=1,2,3,4 ATOMS2=5,6,7,8 ATOMS3=9,10,11,12
phi: TORSION ATOMS1=13,14,15,16 ATOMS2=17,18,19,20 ATOMS3=21,22,23,24

dd: DISPLACEMENT ARG1=psi,phi ARG2=ref_psi,ref_phi
PRINT ARG=dd FILE=colvar
```

The output here is still a $3 \times 2$ matrix. Now, however, each of the three instantaneous angles we have calculated
has its own set of reference values. A different pair of instaneous and reference values is used to calculate each element
of the resulting matrix.

DISPLACEMENT actions that compute $M\times N$ matrices, $D$, are used within the [EUCLIDEAN_DISTANCE](EUCLIDEAN_DISTANCE.md),
[NORMALIZED_EUCLIDEAN_DISTANCE](NORMALIZED_EUCLIDEAN_DISTANCE.md) and [MAHALANOBIS_DISTANCE](MAHALANOBIS_DISTANCE.md) shortcuts.
Doing so is useful as if you take the diagonal elements of a product of matrices that is similar to the product of vectors and matrices that we introduced
earlier:

$$
d = D M D^T
$$

you can calculate $M$ values for the [EUCLIDEAN_DISTANCE](EUCLIDEAN_DISTANCE.md), [NORMALIZED_EUCLIDEAN_DISTANCE](NORMALIZED_EUCLIDEAN_DISTANCE.md)
and [MAHALANOBIS_DISTANCE](MAHALANOBIS_DISTANCE.md).

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace refdist {

class Displacement : public ActionShortcut {
public:
  static std::string fixArgumentDot( const std::string& argin );
  static void registerKeywords( Keywords& keys );
  Value* getValueWithLabel( const std::string& name ) const ;
  explicit Displacement(const ActionOptions&ao);
};

PLUMED_REGISTER_ACTION(Displacement,"DISPLACEMENT")

void Displacement::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords(keys);
  keys.add("compulsory","ARG1","The point that we are calculating the distance from");
  keys.add("compulsory","ARG2","The point that we are calculating the distance to");
  keys.setValueDescription("vector/matrix","the differences between the input arguments");
  keys.needsAction("DIFFERENCE");
  keys.needsAction("TRANSPOSE");
  keys.needsAction("VSTACK");
}

std::string Displacement::fixArgumentDot( const std::string& argin ) {
  std::string argout = argin;
  std::size_t dot=argin.find(".");
  if( dot!=std::string::npos ) {
    argout = argin.substr(0,dot) + "_" + argin.substr(dot+1);
  }
  return argout;
}

Displacement::Displacement( const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  // Read in argument names
  std::vector<std::string> arg1f, arg2f;
  parseVector("ARG1",arg1f);
  parseVector("ARG2",arg2f);
  // Check if one of the input arguments is a reference cluster
  if( arg1f.size()!=arg2f.size() ) {
    error("number of arguments specified to ARG1 should be same as number for ARG2");
  }

  Value* val1=getValueWithLabel( arg1f[0] );
  if( arg1f.size()==1 && val1->getRank()!=0 ) {
    Value* val2=getValueWithLabel( arg2f[0] );
    if( val1->getNumberOfValues()==val2->getNumberOfValues() ) {
      readInputLine( getShortcutLabel() + "_" + fixArgumentDot(arg1f[0]) + "_diff: DIFFERENCE ARG=" + arg1f[0] + "," + arg2f[0] );
      readInputLine( getShortcutLabel() + ": TRANSPOSE ARG=" + getShortcutLabel() + "_" + fixArgumentDot(arg1f[0]) + "_diff");
    } else {
      readInputLine( getShortcutLabel() + ": DIFFERENCE ARG=" + arg1f[0] + "," + arg2f[0] );
    }
  } else {
    for(unsigned i=0; i<arg1f.size(); ++i) {
      readInputLine( getShortcutLabel() + "_" + fixArgumentDot(arg1f[i]) + "_diff: DIFFERENCE ARG=" + arg1f[i] + "," + arg2f[i] );
    }
    std::string argdat = "ARG=" + getShortcutLabel() + "_" + fixArgumentDot(arg1f[0]) + "_diff";
    for(unsigned i=1; i<arg1f.size(); ++i) {
      argdat += "," +  getShortcutLabel() + "_" + fixArgumentDot(arg1f[i]) + "_diff";
    }
    readInputLine( getShortcutLabel() + ": VSTACK " + argdat );
  }
}

Value* Displacement::getValueWithLabel( const std::string& name ) const {
  std::size_t dot=name.find(".");
  std::string sname = name.substr(0,dot);
  ActionWithValue* vv=plumed.getActionSet().selectWithLabel<ActionWithValue*>( sname );
  if( !vv ) {
    error("cannot find value with name " + name );
  }
  if( dot==std::string::npos ) {
    return vv->copyOutput(0);
  }
  if( !vv->exists(name) ) {
    error("cannot find value with name " + name );
  }
  return vv->copyOutput( name );
}

}
}
