/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2020 The plumed team
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
#include "tools/Angle.h"
#include "core/ActionRegister.h"

#include <string>
#include <cmath>

//+PLUMEDOC MATRIX HBOND_MATRIX
/*
Adjacency matrix in which two atoms are adjacent if there is a hydrogen bond between them.

A useful tool for developing complex collective variables is the notion of the
so called adjacency matrix.  An adjacency matrix is an $N \times N$ matrix in which the $i$th, $j$th element tells you whether
or not the $i$th and $j$th atoms/molecules from a set of $N$ atoms/molecules are adjacent or not.  As detailed in the documentation
for [CONTACT_MATRIX](CONTACT_MATRIX.md) there are then a range of further analyses that you can perform on these matrices.

For this action the elements of the adjacency matrix are calculated using:

$$
a_{ij} = \sigma_{oo}( |\mathbf{r}_{ij}| ) \sum_{k=1}^N \sigma_{oh}( |\mathbf{r}_{ik}| ) \sigma_{\theta}( \theta_{kij} )
$$

This expression was derived by thinking about how to detect if there is a hydrogen bond between atoms $i$ and $j$.  The notion is that
if the hydrogen bond is present atoms $i$ and $j$ should be within a certain cutoff distance.  In addition, there should be a hydrogen
within a certain cutoff distance of atom $i$ and this hydrogen should lie on or close to the vector connecting atoms $i$ and $j$.
As such $\sigma_{oo}(r_{ij})$ is a switching function that acts on the modulus of the vector connecting atom $i$ to atom
$j$.  The sum over $k$ then runs over all the hydrogen atoms that are specified using using HYDROGEN keyword.  $\sigma_{oh}(r_{ik})$
is a switching function that acts on the modulus of the vector connecting atom $i$ to atom $k$ and $\sigma_{\theta}(\theta_{kij})$
is a switching function that acts on the angle between the vector connecting atoms $i$ and $j$ and the vector connecting atoms $i$ and
$k$.

It is important to note that hydrogen bonds, unlike regular bonds, are asymmetric. In other words, the hydrogen atom does not sit at the
mid point between the two other atoms in this three-center bond.  As a result of this adjacency matrices calculated using HBOND_MATRIX are not
symmetric like those calculated by [CONTACT_MATRIX](CONTACT_MATRIX.md).

The following input can be used to analyze the number of hydrogen bonds each of the oxygen atoms in a box of water participates in.  Each
water molecule can participate in a hydrogen bond in one of two ways.  It can either donate one of its hydrogen atom to the neighboring oxygen or
it can accept a bond between the hydrogen of a neighboring water molecule and its own oxygen.  The input below allows you to output information
on the number of hydrogen bonds each of the water molecules donates and accepts.  This information is output in two xyz files which each contain
five columns of data.  The first four of these columns are a label for the atom and the x, y and z position of the oxygen.  The last column is then
the number of accepted/donated hydrogen bonds.

```plumed
mat: HBOND_MATRIX ATOMS=1-192:3 HYDROGENS=2-192:3,3-192:3 SWITCH={RATIONAL R_0=3.20} HSWITCH={RATIONAL R_0=2.30} ASWITCH={RATIONAL R_0=0.167pi}
ones: ONES SIZE=64
rsums: MATRIX_VECTOR_PRODUCT ARG=mat,ones
DUMPATOMS ATOMS=1-192:3 ARG=rsums FILE=donors.xyz
```

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace adjmat {

class HbondMatrix {
public:
  SwitchingFunction distanceOOSwitch;
  SwitchingFunction distanceOHSwitch;
  SwitchingFunction angleSwitch;
  static void registerKeywords( Keywords& keys );
  void parseInput( AdjacencyMatrixBase<HbondMatrix>* action );
  static void calculateWeight( const HbondMatrix& data,
                               const AdjacencyMatrixInput& input,
                               MatrixOutput output );
};

typedef AdjacencyMatrixBase<HbondMatrix> hmap;
PLUMED_REGISTER_ACTION(hmap,"HBOND_MATRIX")

void HbondMatrix::registerKeywords( Keywords& keys ) {
  keys.add("atoms-2","DONORS","The list of atoms which can donate a hydrogen bond");
  keys.add("atoms-2","ACCEPTORS","The list of atoms which can accept a hydrogen bond");
  keys.add("atoms","HYDROGENS","The list of atoms that can form the bridge between the two interesting parts "
           "of the structure.");
  keys.add("numbered","SWITCH","The switchingfunction that specifies how close a pair of atoms must be together for there to be a hydrogen bond between them");
  keys.linkActionInDocs("SWITCH","LESS_THAN");
  keys.add("numbered","HSWITCH","The switchingfunction that specifies how close the hydrogen must be to the donor atom of the hydrogen bond for it to be "
           "considered a hydrogen bond");
  keys.linkActionInDocs("HSWITCH","LESS_THAN");
  keys.add("numbered","ASWITCH","A switchingfunction that is used to specify what the angle between the vector connecting the donor atom to the acceptor atom and "
           "the vector connecting the donor atom to the hydrogen must be in order for it considered to be a hydrogen bond");
  keys.linkActionInDocs("ASWITCH","LESS_THAN");
}

void HbondMatrix::parseInput( AdjacencyMatrixBase<HbondMatrix>* action ) {
  std::string errors;
  std::string OOinput;
  action->parse("SWITCH",OOinput);
  if( OOinput.length()==0 ) {
    action->error("could not find SWITCH keyword");
  }
  distanceOOSwitch.set(OOinput,errors);
  if( errors.length()!=0 ) {
    action->error("problem reading SWITCH keyword : " + errors );
  }

  std::string OHinput;
  action->parse("HSWITCH",OHinput);
  if( OHinput.length()==0 ) {
    action->error("could not find HSWITCH keyword");
  }
  distanceOHSwitch.set(OHinput,errors);
  if( errors.length()!=0 ) {
    action->error("problem reading HSWITCH keyword : " + errors );
  }

  std::string anginput;
  action->parse("ASWITCH",anginput);
  if( anginput.length()==0 ) {
    action->error("could not find SWITCH keyword");
  }
  angleSwitch.set(anginput,errors);
  if( errors.length()!=0 ) {
    action->error("problem reading SWITCH keyword : " + errors );
  }

  // Setup link cells
  action->setLinkCellCutoff( false, distanceOOSwitch.get_dmax() );
}

void HbondMatrix::calculateWeight( const HbondMatrix& data,
                                   const AdjacencyMatrixInput& input,
                                   MatrixOutput output ) {
  Vector ood = input.pos;
  double ood_l = ood.modulo2(); // acceptor - donor
  if( ood_l<epsilon) {
    return;
  }
  double ood_df;
  double ood_sw=data.distanceOOSwitch.calculateSqr( ood_l, ood_df );

  for(unsigned i=0; i<input.natoms; ++i) {
    Vector ohd(input.extra_positions[i][0],
               input.extra_positions[i][1],
               input.extra_positions[i][2]);
    double ohd_l=ohd.modulo2();
    double ohd_df;
    double ohd_sw=data.distanceOHSwitch.calculateSqr( ohd_l, ohd_df );

    Angle a;
    Vector ood_adf;
    Vector ohd_adf;
    double angle=a.compute( ood, ohd, ood_adf, ohd_adf );
    double angle_df;
    double angle_sw=data.angleSwitch.calculate( angle, angle_df );
    output.val[0] += ood_sw*ohd_sw*angle_sw;

    if( !input.noderiv ) {
      Vector d1 = angle_sw*ohd_sw*(-ood_df)*ood
                  + angle_sw*ood_sw*(-ohd_df)*ohd
                  + ood_sw*ohd_sw*angle_df*angle*(-ood_adf-ohd_adf);
      output.deriv[0] += d1[0];
      output.deriv[1] += d1[1];
      output.deriv[2] += d1[2];
      Vector d2 = angle_sw*ohd_sw*(+ood_df)*ood
                  + ood_sw*ohd_sw*angle_df*angle*ood_adf;
      output.deriv[3] += d2[0];
      output.deriv[4] += d2[1];
      output.deriv[5] += d2[2];
      Vector d3 = angle_sw*ood_sw*(+ohd_df)*ohd
                  + ood_sw*ohd_sw*angle_df*angle*ohd_adf;
      output.deriv[6+i*3+0] = d3[0];
      output.deriv[6+i*3+1] = d3[1];
      output.deriv[6+i*3+2] = d3[2];
      Tensor vir = angle_sw*ohd_sw*(-ood_df)*Tensor(ood,ood)
                   + angle_sw*ood_sw*(-ohd_df)*Tensor(ohd,ohd)
                   - ood_sw*ohd_sw*angle_df*angle*(Tensor(ood,ood_adf)
                       + Tensor(ohd,ohd_adf));
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
}
