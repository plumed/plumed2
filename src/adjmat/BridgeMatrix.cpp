/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2020 The plumed team
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
#include "AdjacencyMatrixBase.h"
#include "tools/SwitchingFunction.h"
#include "core/ActionRegister.h"

#include <string>
#include <cmath>

namespace PLMD {
namespace adjmat {

//+PLUMEDOC MCOLVAR BRIDGE_MATRIX
/*
Calculate the number of atoms that bridge two parts of a structure

This adjacency matrix is used to implement the [BRIDGE](BRIDGE.md) shortcut. The action outputs a adjacency matrix
in the same way as [CONTACT_MATRIX](CONTACT_MATRIX.md).  However, the  $j,k$ element of the adjacency matrix is calculated
using:

$$
M_{jk} = \sum_i s_A(r_{ij})s_B(r_{ik})
$$

In this expression, the sum runs over all the atoms that were specified using the `BRIDGING_ATOMS` keyword, $s_A$ and
$s_B$ are switching functions, and $r_{ij}$ and $r_{ik}$ are the distances between atom $i$ and $j$ and between atoms
$i$ and $k$.  Less formally, this formula ensures that $j,k$ element of the output matrix is one if there is a bridging
atom between atom $j$ and $k$.

In the following example input atoms 100-200 can serve as bridging atoms between the atoms in GROUPA and GROUPB and the
two switching functions $s_A$ and $s_B$ in the formula above are identical.

```plumed
w1: BRIDGE_MATRIX ...
   BRIDGING_ATOMS=100-200
   GROUPA=1-10 GROUPB=11-20
   SWITCH={RATIONAL R_0=0.2}
...
```

If you use a single GROUP keyword as in the input below as a single SWITCH keyword the output matrix is symmetric.

```plumed
w2: BRIDGE_MATRIX ...
   BRIDGING_ATOMS=100-200 GROUP=1-10
   SWITCH={RATIONAL R_0=0.2}
...
```

However, if the two switching functions are not identical, as in the following example, then the output matrix is __not__ symmetric
even if GROUP is used rather than GROUPA/GROUPB.

```plumed
w2: BRIDGE_MATRIX ...
   BRIDGING_ATOMS=100-200 GROUP=1-10
   SWITCHA={RATIONAL R_0=0.2}
   SWITCHB={RATIONAL R_0=0.4}
...
```

Notice that in all the inputs above the $r_{ij}$ and $r_{ik}$ values that enter the formula above are calculated in a way that takes the
periodic boundary conditions into account.  If you want to ignore the periodic boundary conditions you can use the NOPBC flag as shown below.

```plumed
w2: BRIDGE_MATRIX ...
   BRIDGING_ATOMS=100-200 GROUP=1-10
   SWITCH={RATIONAL R_0=0.2}
   NOPBC
...
```

## COMPONENTS flag

If you add the flag COMPONENTS to the input as shown below:

```plumed
c4: BRIDGE_MATRIX ...
   BRIDGING_ATOMS=100-200 GROUP=1-10
   SWITCH={RATIONAL R_0=0.2}
   COMPONENTS
...
```

then four matrices with the labels `c4.w`, `c4.x`, `c4.y` and `c4.z` are output by the action. The matrix with the label `c4.w` is the adjacency matrix
that would be output if you had not added the COMPONENTS flag. The $i,j$ component of the matrices `c4.x`, `c4.y` and `c4.z` contain the $x$, $y$ and $z$
components of the vector connecting atoms $j$ and $k$. Importantly, however, the components of these vectors are only stored in `c4.x`, `c4.y` and `c4.z`
if the elements of `c4.w` are non-zero. Using the COMPONENTS flag in this way ensures that you can use BRIDGE_MATRIX in tandem with many of the functionalities
that are part of the [symfunc module](module_symfunc.md).

## The MASK keyword

You use the MASK keyword with BRIDGE_MATRIX in the same way that is used in [CONTACT_MATRIX](CONTACT_MATRIX.md).  This keyword thus expects a vector in input,
which tells BRIDGE_MATRIX that it is safe to not calculate certain rows of the output matrix.  An example where this keyword is used is shown below:

```plumed
# The atoms that are of interest
ow: GROUP ATOMS=1-1650
# Fixed virtual atom which serves as the probe volume's center (pos. in nm)
center: FIXEDATOM AT=2.5,2.5,2.5
# Vector in which element i is one if atom i is in sphere of interest and zero otherwise
sphere: INSPHERE ATOMS=ow CENTER=center RADIUS={GAUSSIAN D_0=0.5 R_0=0.01 D_MAX=0.52}
# Calculates cooordination numbers
cmap: BRIDGE_MATRIX ...
   GROUP=ow BRIDGING_ATOMS=1650-3000
   SWITCH={GAUSSIAN D_0=0.32 R_0=0.01 D_MAX=0.34} MASK=sphere
...
ones: ONES SIZE=1650
cc: MATRIX_VECTOR_PRODUCT ARG=cmap,ones
# Multiply coordination numbers by sphere vector
prod: CUSTOM ARG=cc,sphere FUNC=x*y PERIODIC=NO
# Sum of coordination numbers for atoms that are in the sphere of interest
numer: SUM ARG=prod PERIODIC=NO
# Number of atoms that are in sphere of interest
denom: SUM ARG=sphere PERIODIC=NO
# Average coordination number for atoms in sphere of interest
av: CUSTOM ARG=prod,sphere FUNC=x/y PERIODIC=NO
# And print out final CV to a file
PRINT ARG=av FILE=colvar STRIDE=1
```

This input calculates the average for a type of coordination number for those atoms that are within a spherical region that is centered on the point
$(2.5,2.5,2.5)$.  The coordination number for the atoms that is evaluated here counts an atom $k$ in the coordination sphere of atom $j$ if there is a
bridging atom within 0.34 nm of atoms $j$ and $k$.

*/
//+ENDPLUMEDOC

class BridgeMatrix {
public:
  SwitchingFunction sf1;
  SwitchingFunction sf2;
  static void registerKeywords( Keywords& keys );
  void parseInput( AdjacencyMatrixBase<BridgeMatrix>* action );
  static void calculateWeight( const BridgeMatrix& data,
                               const AdjacencyMatrixInput& input,
                               MatrixOutput output );
};

typedef AdjacencyMatrixBase<BridgeMatrix> bmap;
PLUMED_REGISTER_ACTION(bmap,"BRIDGE_MATRIX")

void BridgeMatrix::registerKeywords( Keywords& keys ) {
  keys.add("atoms","BRIDGING_ATOMS","The list of atoms that can form the bridge between the two interesting parts "
           "of the structure.");
  keys.add("optional","SWITCH","The parameters of the two switching functions in the above formula");
  keys.linkActionInDocs("SWITCH","LESS_THAN");
  keys.add("optional","SWITCHA","The switching function on the distance between bridging atoms and the atoms in "
           "group A");
  keys.linkActionInDocs("SWITCHA","LESS_THAN");
  keys.add("optional","SWITCHB","The switching function on the distance between the bridging atoms and the atoms in "
           "group B");
  keys.linkActionInDocs("SWITCHB","LESS_THAN");
}

void BridgeMatrix::parseInput( AdjacencyMatrixBase<BridgeMatrix>* action ) {
  bool oneswitch;
  std::string errors;
  std::string sf1input;
  action->parse("SWITCH",sf1input);
  if( sf1input.length()>0 ) {
    sf1.set(sf1input,errors);
    oneswitch=true;
    if( errors.length()!=0 ) {
      action->error("problem reading SWITCH keyword : " + errors );
    }
    sf2.set(sf1input,errors);
    if( errors.length()!=0 ) {
      action->error("problem reading SWITCH keyword : " + errors );
    }
  } else {
    action->parse("SWITCHA",sf1input);
    if(sf1input.length()>0) {
      sf1.set(sf1input,errors);
      oneswitch=false;
      if( errors.length()!=0 ) {
        action->error("problem reading SWITCHA keyword : " + errors );
      }
      std::string sf2input=sf1input;
      action->parse("SWITCHB",sf2input);
      if(sf2input.length()==0) {
        action->error("found SWITCHA keyword without SWITCHB");
      }
      sf2.set(sf2input,errors);
      if( errors.length()!=0 ) {
        action->error("problem reading SWITCHB keyword : " + errors );
      }
    } else {
      action->error("missing definition of switching functions");
    }
  }
  action->log.printf("  distance between bridging atoms and atoms in GROUPA must be less than %s\n",sf1.description().c_str());
  action->log.printf("  distance between bridging atoms and atoms in GROUPB must be less than %s\n",sf2.description().c_str());

  // Setup link cells
  action->setLinkCellCutoff( oneswitch, sf1.get_dmax() + sf2.get_dmax() );
}

void BridgeMatrix::calculateWeight( const BridgeMatrix& data,
                                    const AdjacencyMatrixInput& input,
                                    MatrixOutput output ) {
  output.val[0] = 0;
  if( input.pos.modulo2()<epsilon ) {
    return;
  }
  for(unsigned i=0; i<input.natoms; ++i) {
    Vector dij(input.extra_positions[i][0],
               input.extra_positions[i][1],
               input.extra_positions[i][2]);
    double dijm = dij.modulo2();
    double dw1, w1=data.sf1.calculateSqr( dijm, dw1 );
    if( dijm<epsilon ) {
      w1=0.0;
      dw1=0.0;
    }
    Vector dik=input.pbc->distance( dij, input.pos );
    double dikm=dik.modulo2();
    double dw2, w2=data.sf2.calculateSqr( dikm, dw2 );
    if( dikm<epsilon ) {
      w2=0.0;
      dw2=0.0;
    }

    output.val[0] += w1*w2;
    // And finish the calculation
    output.deriv[0] += -w2*dw1*dij[0];
    output.deriv[1] += -w2*dw1*dij[1];
    output.deriv[2] += -w2*dw1*dij[2];
    output.deriv[3] += w1*dw2*dik[0];
    output.deriv[4] += w1*dw2*dik[1];
    output.deriv[5] += w1*dw2*dik[2];
    output.deriv[6+i*3+0] = -w1*dw2*dik[0] + w2*dw1*dij[0];
    output.deriv[6+i*3+1] = -w1*dw2*dik[1] + w2*dw1*dij[1];
    output.deriv[6+i*3+2] = -w1*dw2*dik[2] + w2*dw1*dij[2];
    Tensor vir = w1*(-dw2)*Tensor(dik,dik)+w2*(-dw1)*Tensor(dij,dij);
    output.deriv[6 + 3*input.natoms + 0] += vir[0][0];
    output.deriv[6 + 3*input.natoms + 1] += vir[0][1];
    output.deriv[6 + 3*input.natoms + 2] += vir[0][2];
    output.deriv[6 + 3*input.natoms + 3] += vir[1][0];
    output.deriv[6 + 3*input.natoms + 4] += vir[1][1];
    output.deriv[6 + 3*input.natoms + 5] += vir[1][2];
    output.deriv[6 + 3*input.natoms + 6] += vir[2][0];
    output.deriv[6 + 3*input.natoms + 7] += vir[2][1];
    output.deriv[6 + 3*input.natoms + 8] += vir[2][2];
  }
}

}
}
