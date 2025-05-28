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
w1: BRIDGE_MATRIX BRIDGING_ATOMS=100-200 GROUPA=1-10 GROUPB=11-20 SWITCH={RATIONAL R_0=0.2}
```  

*/
//+ENDPLUMEDOC

class BridgeMatrix {
public:
  std::string sf1input, sf2input;
  SwitchingFunction sf1,sf2;
  static void registerKeywords( Keywords& keys );
  void parseInput( AdjacencyMatrixBase<BridgeMatrix>* action );
  BridgeMatrix & operator=( const BridgeMatrix& m ) {
    sf1input=m.sf1input;
    sf2input=m.sf2input;
    std::string errors;
    sf1.set(sf1input,errors);
    sf2.set(sf2input,errors);
    return *this;
  }
  static void calculateWeight( const BridgeMatrix& data, const AdjacencyMatrixInput& input, MatrixOutput& output );
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
  action->parse("SWITCH",sf1input);
  if( sf1input.length()>0 ) {
    sf1.set(sf1input,errors);
    oneswitch=true;
    if( errors.length()!=0 ) {
      action->error("problem reading SWITCH keyword : " + errors );
    }
    sf2input=sf1input;
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

void BridgeMatrix::calculateWeight( const BridgeMatrix& data, const AdjacencyMatrixInput& input, MatrixOutput& output ) {
  output.val[0] = 0;
  if( input.pos.modulo2()<epsilon ) {
    return;
  }
  for(unsigned i=0; i<input.natoms; ++i) {
    Vector dij(input.extra_positions[i][0],input.extra_positions[i][1],input.extra_positions[i][2]);
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
