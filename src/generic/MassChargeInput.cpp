/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2023 The plumed team
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
#include "core/ActionShortcut.h"
#include "core/ActionRegister.h"
#include "core/ActionToPutData.h"
#include "core/ActionAnyorder.h"
#include "core/ActionSetup.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "tools/PDB.h"

namespace PLMD {
namespace generic {

//+PLUMEDOC COLVAR READ_MASS_CHARGE
/*
Set the masses and charges from an input PDB file.

\par Examples

*/
//+ENDPLUMEDOC

class MassChargeInput : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit MassChargeInput(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(MassChargeInput,"READ_MASS_CHARGE")

void MassChargeInput::registerKeywords(Keywords& keys) {
  ActionShortcut::registerKeywords( keys );
  keys.add("compulsory","FILE","a pdb file that contains the masses and charges of the atoms in the beta and occupancy columns");
  keys.addOutputComponent("mass","default","the masses of the atoms in the system");
  keys.addOutputComponent("charges","default","the masses of the atoms in the system");
  keys.needsAction("CONSTANT");
}

MassChargeInput::MassChargeInput(const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao)
{
  const ActionSet& actionset(plumed.getActionSet());
  for(const auto & p : actionset) {
    // check that all the preceding actions are ActionSetup
    if( !dynamic_cast<ActionSetup*>(p.get()) && !dynamic_cast<ActionForInterface*>(p.get()) && !dynamic_cast<ActionAnyorder*>(p.get()) ) {
      error("Action " + getLabel() + " is a setup action, and should be only preceded by other setup actions or by actions that can be used in any order.");
    } else if( (p.get())->getName()=="READ_MASS_CHARGE" ) error("should only be one READ_MASS_CHARGE action in the input file");
  }
  std::string input; parse("FILE",input); PDB pdb;
  if( !pdb.read(input, false, 1.0 ) ) error("error reading pdb file containing masses and charges");
  // Check for correct number of atoms
  unsigned natoms=0; std::vector<ActionToPutData*> inputs=plumed.getActionSet().select<ActionToPutData*>();
  for(const auto & pp : inputs ) {
    if( pp->getRole()=="x" ) natoms = (pp->copyOutput(0))->getShape()[0];
  }
  if( natoms!=pdb.size() ) error("mismatch between number of atoms passed from MD code and number of atoms in PDB file");

  // Now get masses and charges
  std::string nnn, charges, masses;
  Tools::convert( pdb.getBeta()[0], charges );
  Tools::convert( pdb.getOccupancy()[0], masses );
  for(unsigned i=1; i<pdb.size(); ++i) {
    Tools::convert( pdb.getBeta()[i], nnn ); charges += "," + nnn;
    Tools::convert( pdb.getOccupancy()[i], nnn ); masses += "," + nnn;
  }
  // And create constant actions to hold masses and charges
  readInputLine( getShortcutLabel() + "_mass: CONSTANT VALUES=" + masses );
  readInputLine( getShortcutLabel() + "_charges: CONSTANT VALUES=" + charges );
}

}
}
