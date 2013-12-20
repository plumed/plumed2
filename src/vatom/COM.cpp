/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

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
#include "core/PlumedMain.h"
#include "core/Atoms.h"

using namespace std;

namespace PLMD{
namespace vatom{

//+PLUMEDOC VATOM COM
/*
Calculate the center of mass for a group of atoms.

The computed
center of mass is stored as a virtual atom that can be accessed in
an atom list through the label for the COM action that creates it.

For arbitrary weights (e.g. geometric center) see \ref CENTER.

\par Examples

The following input instructs plumed to print the distance between the
center of mass for atoms 1,2,3,4,5,6,7 and that for atoms 15,20:
\verbatim
COM ATOMS=1-7         LABEL=c1
COM ATOMS=15,20       LABEL=c2
DISTANCE ATOMS=c1,c2  LABEL=d1
PRINT ARG=d1
\endverbatim
(See also \ref DISTANCE and \ref PRINT).

*/
//+ENDPLUMEDOC


class COM:
  public ActionWithVirtualAtom
{
public:
  COM(const ActionOptions&ao);
  void calculate();
  static void registerKeywords( Keywords& keys );
};

PLUMED_REGISTER_ACTION(COM,"COM")

void COM::registerKeywords(Keywords& keys){
  ActionWithVirtualAtom::registerKeywords(keys);
}

COM::COM(const ActionOptions&ao):
  Action(ao),
  ActionWithVirtualAtom(ao)
{
  vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  if(atoms.size()==0) error("at least one atom should be specified");
  checkRead();
  log.printf("  of atoms");
  for(unsigned i=0;i<atoms.size();++i) log.printf(" %d",atoms[i].serial());
  log.printf("\n");
  requestAtoms(atoms);
}

void COM::calculate(){
  Vector pos;
  double mass(0.0);
  vector<Tensor> deriv(getNumberOfAtoms());
  for(unsigned i=0;i<getNumberOfAtoms();i++) mass+=getMass(i);
  if( plumed.getAtoms().chargesWereSet() ){
     double charge(0.0);
     for(unsigned i=0;i<getNumberOfAtoms();i++) charge+=getCharge(i);
     setCharge(charge);
  } else {
     setCharge(0.0);
  }
  for(unsigned i=0;i<getNumberOfAtoms();i++){
    pos+=(getMass(i)/mass)*getPosition(i);
    deriv[i]=(getMass(i)/mass)*Tensor::identity();
  }
  setPosition(pos);
  setMass(mass);
  setAtomsDerivatives(deriv);
}

}
}
