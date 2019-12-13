/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2019 The plumed team
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

//+PLUMEDOC VATOM FIXEDATOM
/*
Add a virtual atom in a fixed position.

This action creates a virtual atom at a fixed position.
The coordinates can be specified in Cartesian components (by default)
or in scaled coordinates (SCALED_COMPONENTS).
It is also possible to assign a predefined charge or mass to the atom.

\attention
Similar to \ref POSITION this variable is not invariant for translation
of the system. Adding a force on it can create serious troubles.

Notice that the distance between to atoms created
using FIXEDATOM is invariant for translation.
Additionally, if one first align atoms to a reference using \ref FIT_TO_TEMPLATE,
then it is safe to add further fixed atoms without breaking translational invariance.

\par Examples

The following input instructs plumed to compute the angle between
distance of atoms 15 and 20 and the z axis and keeping it close to zero.
\plumedfile
a: FIXEDATOM AT=0,0,0
b: FIXEDATOM AT=0,0,1
an: ANGLE ATOMS=a,b,15,20
RESTRAINT ARG=an AT=0.0 KAPPA=100.0
\endplumedfile

The following input instructs plumed to align a protein to a template
and to then compute the distance between one of the atoms in the protein and the point
(10,20,30).
\plumedfile
FIT_TO_TEMPLATE STRIDE=1 REFERENCE=ref.pdb TYPE=SIMPLE
a: FIXEDATOM AT=10,20,30
d: DISTANCE ATOMS=a,20
PRINT ARG=d FILE=colvar
\endplumedfile

The reference structure to align to is provided in a pdb file called ref.pdb as shown below:

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


class FixedAtom:
  public ActionWithVirtualAtom
{
  Vector coord;
  double mass,charge;
  bool scaled_components;
public:
  explicit FixedAtom(const ActionOptions&ao);
  void calculate() override;
  static void registerKeywords( Keywords& keys );
};

PLUMED_REGISTER_ACTION(FixedAtom,"FIXEDATOM")

void FixedAtom::registerKeywords(Keywords& keys) {
  ActionWithVirtualAtom::registerKeywords(keys);
  keys.add("compulsory","AT","coordinates of the virtual atom");
  keys.add("compulsory","SET_MASS","1","mass of the virtual atom");
  keys.add("compulsory","SET_CHARGE","0","charge of the virtual atom");
  keys.addFlag("SCALED_COMPONENTS",false,"use scaled components");
}

FixedAtom::FixedAtom(const ActionOptions&ao):
  Action(ao),
  ActionWithVirtualAtom(ao)
{
  vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  if(atoms.size()!=0) error("ATOMS should be empty");

  parseFlag("SCALED_COMPONENTS",scaled_components);

  vector<double> at;
  parseVector("AT",at);
  if(at.size()!=3) error("AT should be a list of three real numbers");

  parse("SET_MASS",mass);
  parse("SET_CHARGE",charge);

  coord[0]=at[0];
  coord[1]=at[1];
  coord[2]=at[2];

  checkRead();
  log<<"  AT position "<<coord[0]<<" "<<coord[1]<<" "<<coord[2]<<"\n";
  if(scaled_components) log<<"  position is in scaled components\n";
}

void FixedAtom::calculate() {
  vector<Tensor> deriv(getNumberOfAtoms());
  if(scaled_components) {
    setPosition(getPbc().scaledToReal(coord));
  } else {
    setPosition(coord);
  }
  setMass(mass);
  setCharge(charge);
  setAtomsDerivatives(deriv);
// Virial contribution
  if(!scaled_components) setBoxDerivativesNoPbc();
// notice that with scaled components there is no additional virial contribution
}

}
}
