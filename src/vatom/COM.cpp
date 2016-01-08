/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2016 The plumed team
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

When running with periodic boundary conditions, the atoms should be 
in the proper periodic image. This is done automatically since PLUMED 2.2,
by considering the ordered list of atoms and rebuilding PBCs with a procedure
that is equivalent to that done in \ref WHOLEMOLECULES . Notice that
rebuilding is local to this action. This is different from \ref WHOLEMOLECULES
which actually modifies the coordinates stored in PLUMED. 

In case you want to recover the old behavior you should use the NOPBC flag.
In that case you need to take care that atoms are in the correct
periodic image.

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
  bool nopbc;
public:
  explicit COM(const ActionOptions&ao);
  void calculate();
  static void registerKeywords( Keywords& keys );
};

PLUMED_REGISTER_ACTION(COM,"COM")

void COM::registerKeywords(Keywords& keys){
  ActionWithVirtualAtom::registerKeywords(keys);
  keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions when calculating distances");
}

COM::COM(const ActionOptions&ao):
  Action(ao),
  ActionWithVirtualAtom(ao),
  nopbc(false)
{
  vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  if(atoms.size()==0) error("at least one atom should be specified");
  parseFlag("NOPBC",nopbc);
  checkRead();
  log.printf("  of atoms");
  for(unsigned i=0;i<atoms.size();++i) log.printf(" %d",atoms[i].serial());
  log.printf("\n");
  if(!nopbc){
    log<<"  PBC will be ignored\n";
  } else {
    log<<"  broken molecules will be rebuilt assuming atoms are in the proper order\n";
  }
  requestAtoms(atoms);
}

void COM::calculate(){
  Vector pos;
  if(!nopbc) makeWhole();
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
