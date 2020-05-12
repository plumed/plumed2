/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2020 The plumed team
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
#include "ActionWithVirtualAtom.h"
#include "ActionRegister.h"
#include "tools/Vector.h"
#include "tools/Exception.h"

using namespace std;

namespace PLMD {
namespace vatom {

//+PLUMEDOC VATOM GHOST
/*
Calculate the absolute position of a ghost atom with fixed coordinates in the local reference frame formed by three atoms.

The computed ghost atom is stored as a virtual atom that can be accessed in
 an atom list through the the label for the GHOST action that creates it.

\par Examples

The following input instructs plumed to print the distance between the
ghost atom and the center of mass for atoms 15,20:
\plumedfile
c1: GHOST ATOMS=1,5,10 COORDINATES=10.0,10.0,10.0
c2: COM ATOMS=15,20
d1: DISTANCE ATOMS=c1,c2
PRINT ARG=d1
\endplumedfile

*/
//+ENDPLUMEDOC


class Ghost:
  public ActionWithVirtualAtom
{
  vector<double> coord;
public:
  explicit Ghost(const ActionOptions&ao);
  void calculate() override;
  static void registerKeywords( Keywords& keys );
};

PLUMED_REGISTER_ACTION(Ghost,"GHOST")

void Ghost::registerKeywords(Keywords& keys) {
  ActionWithVirtualAtom::registerKeywords(keys);
  keys.add("atoms","COORDINATES","coordinates of the ghost atom in the local reference frame");
}

Ghost::Ghost(const ActionOptions&ao):
  Action(ao),
  ActionWithVirtualAtom(ao)
{
  vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  if(atoms.size()!=3) error("ATOMS should contain a list of three atoms");

  parseVector("COORDINATES",coord);
  if(coord.size()!=3) error("COORDINATES should be a list of three real numbers");

  checkRead();
  log.printf("  of atoms");
  for(unsigned i=0; i<atoms.size(); ++i) log.printf(" %d",atoms[i].serial());
  log.printf("\n");
  requestAtoms(atoms);
}

void Ghost::calculate() {
  Vector pos;
  vector<Tensor> deriv(getNumberOfAtoms());
  vector<Vector> n;

// first versor
  Vector n01 = delta(getPosition(0), getPosition(1));
  n.push_back(n01/n01.modulo());

// auxiliary vector
  Vector n02 = delta(getPosition(0), getPosition(2));

// second versor
  Vector n03 = crossProduct(n[0],n02);
  double n03_norm = n03.modulo();
  n.push_back(n03/n03_norm);

// third versor
  n.push_back(crossProduct(n[0],n[1]));

// origin of the reference system
  pos = getPosition(0);

  for(unsigned i=0; i<3; ++i) {
    pos += coord[i] * n[i];
  }

  setPosition(pos);
  setMass(1.0);
  setCharge(0.0);

// some useful tensors for derivatives
  Tensor dn0d0  = (-Tensor::identity()+Tensor(n[0],n[0]))/n01.modulo();
  Tensor dn0d1  = (+Tensor::identity()-Tensor(n[0],n[0]))/n01.modulo();
  Tensor dn02d0 = -Tensor::identity();
  Tensor dn02d2 =  Tensor::identity();

// derivative of n1 = n0 x n02
  Tensor dn1d0, dn1d1, dn1d2;
  Vector aux0, aux1, aux2;

  for(unsigned j=0; j<3; ++j) {
// derivative of n0 x n02 with respect to point 0, coordinate j
    Vector tmp00  = Vector( dn0d0(j,0),  dn0d0(j,1),  dn0d0(j,2));
    Vector tmp020 = Vector(dn02d0(j,0), dn02d0(j,1), dn02d0(j,2));
    Vector tmp0   = crossProduct(tmp00,n02) + crossProduct(n[0],tmp020);
    aux0[j]       = dotProduct(tmp0,n[1]);
// derivative of n0 x n02 with respect to point 1, coordinate j
    Vector tmp01  = Vector( dn0d1(j,0),  dn0d1(j,1),  dn0d1(j,2));
    Vector tmp1   = crossProduct(tmp01,n02);
    aux1[j]       = dotProduct(tmp1,n[1]);
// derivative of n0 x n02 with respect to point 2, coordinate j
    Vector tmp022 = Vector(dn02d2(j,0), dn02d2(j,1), dn02d2(j,2));
    Vector tmp2   = crossProduct(n[0],tmp022);
    aux2[j]       = dotProduct(tmp2,n[1]);
// derivative of n1 = (n0 x n02) / || (n0 x n02) ||
    for(unsigned i=0; i<3; ++i) {
      dn1d0(j,i) = ( tmp0[i] - aux0[j] * n[1][i] ) / n03_norm;
      dn1d1(j,i) = ( tmp1[i] - aux1[j] * n[1][i] ) / n03_norm;
      dn1d2(j,i) = ( tmp2[i] - aux2[j] * n[1][i] ) / n03_norm;
    }
  }

// Derivative of the last versor n2 = n0 x n1 =  ( n0( n0 n02 ) - n02 ) / || n0 x n02 ||
// Scalar product and derivatives
  double n0_n02 = dotProduct(n[0],n02);
  Vector dn0_n02d0, dn0_n02d1, dn0_n02d2;

  for(unsigned j=0; j<3; ++j) {
    for(unsigned i=0; i<3; ++i) {
      dn0_n02d0[j] += dn0d0(j,i)*n02[i] + n[0][i]*dn02d0(j,i);
      dn0_n02d1[j] += dn0d1(j,i)*n02[i];
      dn0_n02d2[j] +=                     n[0][i]*dn02d2(j,i);
    }
  }

  Tensor dn2d0, dn2d1, dn2d2;
  for(unsigned j=0; j<3; ++j) {
    for(unsigned i=0; i<3; ++i) {
      dn2d0(j,i) = ( dn0d0(j,i) * n0_n02 + n[0][i] * dn0_n02d0[j] - dn02d0(j,i) - ( n[0][i] * n0_n02 - n02[i] ) * aux0[j] / n03_norm ) / n03_norm;
      dn2d1(j,i) = ( dn0d1(j,i) * n0_n02 + n[0][i] * dn0_n02d1[j]               - ( n[0][i] * n0_n02 - n02[i] ) * aux1[j] / n03_norm ) / n03_norm;
      dn2d2(j,i) = (                       n[0][i] * dn0_n02d2[j] - dn02d2(j,i) - ( n[0][i] * n0_n02 - n02[i] ) * aux2[j] / n03_norm ) / n03_norm;
    }
  }

// Finally, the derivative tensor
  deriv[0] = Tensor::identity() + coord[0]*dn0d0 + coord[1]*dn1d0 + coord[2]*dn2d0;
  deriv[1] =                      coord[0]*dn0d1 + coord[1]*dn1d1 + coord[2]*dn2d1;
  deriv[2] =                                       coord[1]*dn1d2 + coord[2]*dn2d2;

  setAtomsDerivatives(deriv);

// Virial contribution
  setBoxDerivativesNoPbc();
}

}
}
