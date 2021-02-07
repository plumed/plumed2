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
#include "core/ActionAtomistic.h"
#include "core/ActionPilot.h"
#include "core/ActionRegister.h"
#include "tools/Vector.h"
#include "tools/Matrix.h"
#include "tools/AtomNumber.h"
#include "tools/Tools.h"
#include "core/Atoms.h"
#include "tools/Pbc.h"

namespace PLMD {
namespace generic {

//+PLUMEDOC GENERIC RESET_CELL
/*
This action is used to rotate the full cell

This can be used to modify the periodic box. Notice that
this is done at fixed scaled coordinates,
so that also atomic coordinates for the entire system are affected.
To see what effect try
the \ref DUMPATOMS directive to output the atomic positions.

Also notice that PLUMED propagate forces correctly so that you can add a bias on a CV computed
after rotation. See also \ref FIT_TO_TEMPLATE

Currently, only TYPE=TRIANGULAR is implemented, which allows one to reset
the cell to a lower triangular one. Namely, a proper rotation is found that allows
rotating the box so that the first lattice vector is in the form (ax,0,0),
the second lattice vector is in the form (bx,by,0), and the third lattice vector is
arbitrary.

\attention
The implementation of this action is available but should be considered in testing phase. Please report any
strange behavior.

\attention
This directive modifies the stored position at the precise moment
it is executed. This means that only collective variables
which are below it in the input script will see the corrected positions.
Unless you
know exactly what you are doing, leave the default stride (1), so that
this action is performed at every MD step.

\par Examples

Reset cell to be triangular after a rototranslational fit
\plumedfile
DUMPATOMS FILE=dump-original.xyz ATOMS=1-20
FIT_TO_TEMPLATE STRIDE=1 REFERENCE=ref.pdb TYPE=OPTIMAL
DUMPATOMS FILE=dump-fit.xyz ATOMS=1-20
RESET_CELL TYPE=TRIANGULAR
DUMPATOMS FILE=dump-reset.xyz ATOMS=1-20
\endplumedfile

The reference file for the FIT_TO_TEMPLATE is just a normal pdb file with the format shown below:

\auxfile{ref.pdb}
ATOM      8  HT3 ALA     2      -1.480  -1.560   1.212  1.00  1.00      DIA  H
ATOM      9  CAY ALA     2      -0.096   2.144  -0.669  1.00  1.00      DIA  C
ATOM     10  HY1 ALA     2       0.871   2.385  -0.588  1.00  1.00      DIA  H
ATOM     12  HY3 ALA     2      -0.520   2.679  -1.400  1.00  1.00      DIA  H
ATOM     14  OY  ALA     2      -1.139   0.931  -0.973  1.00  1.00      DIA  O
END
\endauxfile

*/
//+ENDPLUMEDOC


class ResetCell:
  public ActionPilot,
  public ActionAtomistic
{
  std::string type;
  Tensor rotation,newbox;

public:
  explicit ResetCell(const ActionOptions&ao);
  static void registerKeywords( Keywords& keys );
  void calculate() override;
  void apply() override;
};

PLUMED_REGISTER_ACTION(ResetCell,"RESET_CELL")

void ResetCell::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys );
  keys.add("compulsory","STRIDE","1","the frequency with which molecules are reassembled.  Unless you are completely certain about what you are doing leave this set equal to 1!");
  keys.add("compulsory","TYPE","TRIANGULAR","the manner in which the cell is reset");
}

ResetCell::ResetCell(const ActionOptions&ao):
  Action(ao),
  ActionPilot(ao),
  ActionAtomistic(ao)
{
  type.assign("TRIANGULAR");
  parse("TYPE",type);

  log<<"  type: "<<type<<"\n";
  if(type!="TRIANGULAR") error("undefined type "+type);

  checkRead();
}


void ResetCell::calculate() {

  Pbc & pbc(modifyGlobalPbc());

  Tensor box=pbc.getBox();

// moduli of lattice vectors
  double a=modulo(box.getRow(0));
  double b=modulo(box.getRow(1));
  double c=modulo(box.getRow(2));
// cos-angle between lattice vectors
  double ab=dotProduct(box.getRow(0),box.getRow(1))/(a*b);
  double ac=dotProduct(box.getRow(0),box.getRow(2))/(a*c);
  double bc=dotProduct(box.getRow(1),box.getRow(2))/(b*c);

// generate a new set of lattice vectors as a lower triangular matrix
  newbox[0][0]=a;
  newbox[1][0]=b*ab;
  newbox[1][1]=std::sqrt(b*b-newbox[1][0]*newbox[1][0]);
  newbox[2][0]=c*ac;
  newbox[2][1]=c*(bc-ac*ab)/std::sqrt(1-ab*ab);
  newbox[2][2]=std::sqrt(c*c-newbox[2][0]*newbox[2][0]-newbox[2][1]*newbox[2][1]);

  if(determinant(newbox)*determinant(box)<0) newbox[2][2]=-newbox[2][2];

// rotation matrix from old to new coordinates
  rotation=transpose(matmul(inverse(box),newbox));

// rotate all coordinates
  for(unsigned i=0; i<getTotAtoms(); i++) {
    Vector & ato (modifyGlobalPosition(AtomNumber::index(i)));
    ato=matmul(rotation,ato);
  }
// rotate box
  pbc.setBox(newbox);
}

void ResetCell::apply() {
// rotate back forces
  for(unsigned i=0; i<getTotAtoms(); i++) {
    Vector & f(modifyGlobalForce(AtomNumber::index(i)));
    f=matmul(transpose(rotation),f);
  }

  Tensor& virial(modifyGlobalVirial());
// I have no mathematical derivation for this.
// The reasoning is the following.
// virial= h^T * dU/dh, where h is the box matrix and dU/dh its derivatives.
// The final virial should be rotationally invariant, that is symmetric.
// in the rotated frame, dU/dh elements [0][1], [0][2], and [1][2] should
// be changed so as to enforce rotational invariance. Thus we here have to
// make the virial matrix symmetric.
// Since h^T is upper triangular, it can be shown that any change in these elements
// will only affect the corresponding elements of the virial matrix.
// Thus, the only possibility is to set the corresponding elements
// of the virial matrix equal to their symmetric ones.
// GB
  virial[0][1]=virial[1][0];
  virial[0][2]=virial[2][0];
  virial[1][2]=virial[2][1];
// rotate back virial
  virial=matmul(transpose(rotation),matmul(virial,rotation));



}

}
}
