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
#include "tools/IFile.h"
#include "tools/PDB.h"

namespace PLMD {
namespace generic {

//+PLUMEDOC COLVAR READMASSCHARGE
/*
Set the masses and charges from an input PDB file.

You can use this command to read in the masses and charges for the atoms as shown in the example input below:

```plumed
#SETTINGS INPUTFILES=regtest/basic/rt-readmasscharge2/mcinpt
mq: READMASSCHARGE FILE=regtest/basic/rt-readmasscharge2/mcinpt
```

In the above example the masses and charges are in a file called `mcinpt` that contains three columns of numbers.
The first column of numbers is a numerical index, the second column is the masses and the third is the charges.

You can read masses and charges from a PDB file instead by using an input like the one shown below:

```plumed
#SETTINGS NATOMS=108 INPUTFILES=regtest/basic/rt-readmasscharge/test.pdb
mq: READMASSCHARGE PDBFILE=regtest/basic/rt-readmasscharge/test.pdb
```

In this chages masses are read from the occupancy column and charges are read from the beta column.

__Using this command to read in masses and charges is useful when you are postprocessing a trajectory using driver.
If you are using PLUMED as a plugin to a molecular dynamics code the masses and charges will often be passed from the
MD code to PLUMED.__


*/
//+ENDPLUMEDOC

class MassChargeInput : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit MassChargeInput(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(MassChargeInput,"READMASSCHARGE")

void MassChargeInput::registerKeywords(Keywords& keys) {
  ActionShortcut::registerKeywords( keys );
  keys.add("optional","FILE","input file that contains the masses and charges that should be used");
  keys.add("compulsory","PDBFILE","a pdb file that contains the masses and charges of the atoms in the beta and occupancy columns");
  keys.addOutputComponent("mass","default","vector","the masses of the atoms in the system");
  keys.addOutputComponent("charges","default","vector","the masses of the atoms in the system");
  keys.needsAction("CONSTANT");
}

MassChargeInput::MassChargeInput(const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  const ActionSet& actionset(plumed.getActionSet());
  for(const auto & p : actionset) {
    // check that all the preceding actions are ActionSetup
    if( !dynamic_cast<ActionSetup*>(p.get()) && !dynamic_cast<ActionForInterface*>(p.get()) && !dynamic_cast<ActionAnyorder*>(p.get()) ) {
      error("Action " + getLabel() + " is a setup action, and should be only preceded by other setup actions or by actions that can be used in any order.");
    } else if( (p.get())->getName()=="READMASSCHARGE" ) {
      error("should only be one READMASSCHARGE action in the input file");
    }
  }
  // Check for correct number of atoms
  unsigned natoms=0;
  std::vector<ActionToPutData*> inputs=plumed.getActionSet().select<ActionToPutData*>();
  for(const auto & pp : inputs ) {
    if( pp->getRole()=="x" ) {
      natoms = (pp->copyOutput(0))->getShape()[0];
    }
  }
  std::string input;
  parse("FILE",input);
  std::vector<double> masses( natoms ), charges( natoms );
  if( input.length()>0 ) {
    log.printf("   reading masses and charges from file named %s \n", input.c_str() );
    IFile ifile;
    ifile.open( input );
    int index;
    double mass;
    double charge;
    while(ifile.scanField("index",index).scanField("mass",mass).scanField("charge",charge).scanField()) {
      if( static_cast<unsigned>(index)>=natoms ) {
        error("indices of atoms in input file are too large");
      }
      masses[index]=mass;
      charges[index]=charge;
    }
    ifile.close();
  } else {
    std::string pdbinpt;
    parse("PDBFILE",pdbinpt);
    PDB pdb;
    log.printf("  reading masses and charges from pdb file named %s \n", pdbinpt.c_str() );
    if( !pdb.read(pdbinpt, false, 1.0 ) ) {
      error("error reading pdb file containing masses and charges");
    }
    if( natoms!=pdb.size() ) {
      error("mismatch between number of atoms passed from MD code and number of atoms in PDB file");
    }
    masses = pdb.getOccupancy();
    charges = pdb.getBeta();
  }

  // Now get masses and charges
  std::string nnn, qstr, mstr;
  Tools::convert( masses[0], mstr );
  Tools::convert( charges[0], qstr );
  for(unsigned i=1; i<natoms; ++i) {
    Tools::convert( masses[i], nnn );
    mstr += "," + nnn;
    Tools::convert( charges[i], nnn );
    qstr += "," + nnn;
  }
  // And create constant actions to hold masses and charges
  readInputLine( getShortcutLabel() + "_mass: CONSTANT VALUES=" + mstr );
  readInputLine( getShortcutLabel() + "_charges: CONSTANT VALUES=" + qstr );
}

}
}
