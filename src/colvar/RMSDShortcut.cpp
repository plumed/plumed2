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
#include "core/ActionShortcut.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "tools/PDB.h"
#include "RMSD.h"

namespace PLMD {
namespace colvar {

class RMSDShortcut : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit RMSDShortcut(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(RMSDShortcut,"RMSD")

void RMSDShortcut::registerKeywords(Keywords& keys){
  ActionShortcut::registerKeywords( keys );
  keys.add("compulsory","REFERENCE","a file in pdb format containing the reference structure and the atoms involved in the CV");
  keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions when calculating distances");
  keys.addFlag("NUMERICAL_DERIVATIVES", false, "calculate the derivatives for these quantities numerically");
  RMSD::registerRMSD( keys );
}

RMSDShortcut::RMSDShortcut(const ActionOptions& ao):
Action(ao),
ActionShortcut(ao)
{
  std::string reference; parse("REFERENCE",reference);
  // Read the reference pdb file
  PDB pdb; if( !pdb.read(reference,plumed.getAtoms().usingNaturalUnits(),0.1/plumed.getAtoms().getUnits().getLength()) ) plumed_merror("missing file " + reference );
  // Create the input reference position
  readInputLine( getShortcutLabel() + "_ref: READ_CONFIG CHECKATOMS REFERENCE=" + reference );
  // Create the input for the RMSD_CALC object
  std::vector<AtomNumber> anum( pdb.getAtomNumbers() ); std::string num; Tools::convert( anum[0].serial(), num ); 
  std::string rmsd_line = getShortcutLabel() + ": RMSD_CALC REFERENCE_ATOMS=" + getShortcutLabel() + "_ref ATOMS=" + num;
  for(unsigned i=1; i<anum.size(); ++i) { Tools::convert( anum[i].serial(), num ); rmsd_line += "," + num; }
  // Now align
  std::vector<double> align( pdb.getOccupancy() ); Tools::convert( align[0], num ); rmsd_line += " ALIGN=" + num; 
  for(unsigned i=1; i<align.size(); ++i) { Tools::convert( align[i], num ); rmsd_line += "," + num; }
  // And displace
  std::vector<double> displace( pdb.getBeta() ); Tools::convert( displace[0], num ); rmsd_line += " DISPLACE=" + num; 
  for(unsigned i=1; i<displace.size(); ++i) { Tools::convert( displace[i], num ); rmsd_line += "," + num; }
  // And create the RMSD object
  bool numder; parseFlag("NUMERICAL_DERIVATIVES",numder); if(numder) rmsd_line += " NUMERICAL_DERIVATIVES";
  bool nopbc; parseFlag("NOPBC",nopbc); if(nopbc) rmsd_line += " NOPBC";
  bool squared; parseFlag("SQUARED",squared); if(squared) rmsd_line += " SQUARED";
  bool disp; parseFlag("DISPLACEMENT",disp); if(disp) rmsd_line += " DISPLACEMENT";
  std::string tt; parse("TYPE",tt); readInputLine( rmsd_line + " TYPE=" + tt );
}

}
}
