/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2018 The plumed team
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
#include "MultiColvarShortcuts.h"

//+PLUMEDOC MCOLVAR COORD_ANGLES
/*
Calculate all the angles between bonds in the first coordination spheres of a set of atoms

\par Examples


*/
//+ENDPLUMEDOC

namespace PLMD {
namespace multicolvar {

class CoordAngles : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  static void pruneShortcuts(Keywords& keys);
  explicit CoordAngles(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(CoordAngles,"COORD_ANGLES")

void CoordAngles::registerKeywords(Keywords& keys) {
  ActionShortcut::registerKeywords( keys );
  keys.add("atoms","CATOMS","all the angles between the bonds that radiate out from these central atom are computed");
  keys.add("atoms","GROUP","a list of angls between pairs of bonds connecting one of the atoms specified using the CATOM command and two of the atoms specified here are computed");
  keys.add("compulsory","SWITCH","the switching function specifies that only those bonds that have a length that is less than a certain threshold are considered");
  MultiColvarShortcuts::shortcutKeywords( keys ); pruneShortcuts( keys );
  keys.needsAction("DISTANCE"); keys.needsAction("COMBINE"); keys.needsAction("CUSTOM");
  keys.needsAction("VSTACK"); keys.needsAction("TRANSPOSE");
  keys.needsAction("OUTER_PRODUCT"); keys.needsAction("MATRIX_PRODUCT");
}

void CoordAngles::pruneShortcuts(Keywords& keys) {
  keys.remove("ALT_MIN"); keys.remove("MIN"); keys.remove("MAX"); keys.remove("HIGHEST"); keys.remove("LOWEST");
}

CoordAngles::CoordAngles(const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao)
{
  // Parse the central atoms
  std::vector<std::string> catoms; parseVector("CATOMS",catoms); Tools::interpretRanges(catoms);
  // Parse the coordination sphere
  std::vector<std::string> group; parseVector("GROUP",group); Tools::interpretRanges(group);
  // Create the list of atoms
  std::string atlist; unsigned k=1;
  for(unsigned i=0; i<catoms.size(); ++i) {
    for(unsigned j=0; j<group.size(); ++j) { std::string num; Tools::convert( k, num ); atlist += " ATOMS" + num + "=" + catoms[i] + "," + group[j]; k++; }
  }
  // Calculate the distances
  readInputLine( getShortcutLabel() + "_dd: DISTANCE" + atlist );
  // Transform with the switching function
  std::string switch_input; parse("SWITCH",switch_input);
  readInputLine( getShortcutLabel() + "_sw: LESS_THAN ARG=" + getShortcutLabel() + "_dd SWITCH={" + switch_input +"}");
  // Now get the normalised vectors
  readInputLine( getShortcutLabel() + "_comp: DISTANCE" + atlist + " COMPONENTS");
  readInputLine( getShortcutLabel() + "_norm2: COMBINE ARG=" + getShortcutLabel() + "_comp.x" + "," + getShortcutLabel() + "_comp.y," + getShortcutLabel() + "_comp.z POWERS=2,2,2 PERIODIC=NO");
  readInputLine( getShortcutLabel() + "_norm: CUSTOM ARG=" + getShortcutLabel() + "_norm2 FUNC=sqrt(x) PERIODIC=NO");
  readInputLine( getShortcutLabel() + "_norm_x: CUSTOM ARG=" + getShortcutLabel() + "_comp.x," + getShortcutLabel() + "_norm FUNC=x/y PERIODIC=NO");
  readInputLine( getShortcutLabel() + "_norm_y: CUSTOM ARG=" + getShortcutLabel() + "_comp.y," + getShortcutLabel() + "_norm FUNC=x/y PERIODIC=NO");
  readInputLine( getShortcutLabel() + "_norm_z: CUSTOM ARG=" + getShortcutLabel() + "_comp.z," + getShortcutLabel() + "_norm FUNC=x/y PERIODIC=NO");
  readInputLine( getShortcutLabel() + "_stack: VSTACK ARG=" + getShortcutLabel() + "_norm_x" + "," + getShortcutLabel() + "_norm_y," + getShortcutLabel() + "_norm_z");
  readInputLine( getShortcutLabel() + "_stackT: TRANSPOSE ARG=" + getShortcutLabel() + "_stack");
  // Create the matrix of weights
  readInputLine( getShortcutLabel() + "_swd: OUTER_PRODUCT ELEMENTS_ON_DIAGONAL_ARE_ZERO ARG=" + getShortcutLabel() + "_sw," + getShortcutLabel() + "_sw");
  // Avoid double counting
  readInputLine( getShortcutLabel() + "_wmat: CUSTOM ARG=" + getShortcutLabel() + "_swd FUNC=0.5*x PERIODIC=NO");
  // And the matrix of dot products and the angles
  readInputLine( getShortcutLabel() + "_dpmat: MATRIX_PRODUCT ARG=" + getShortcutLabel() + "_stack," + getShortcutLabel() + "_stackT");
  readInputLine( getShortcutLabel() + "_angles: CUSTOM ARG=" + getShortcutLabel() + "_dpmat FUNC=acos(x) PERIODIC=NO");
  // Read the input
  Keywords keys; MultiColvarShortcuts::shortcutKeywords( keys ); pruneShortcuts( keys ); bool do_mean; parseFlag("MEAN",do_mean);
  std::map<std::string,std::string> keymap; readShortcutKeywords( keys, keymap ); if( do_mean ) keymap.insert(std::pair<std::string,std::string>("SUM",""));
  MultiColvarShortcuts::expandFunctions( getShortcutLabel(), getShortcutLabel() + "_angles", getShortcutLabel() + "_wmat", keymap, this );
  if( do_mean ) {
    readInputLine( getShortcutLabel() + "_denom: SUM ARG=" + getShortcutLabel() + "_wmat PERIODIC=NO");
    readInputLine( getShortcutLabel() + "_mean: CUSTOM ARG=" + getShortcutLabel() + "_sum," + getShortcutLabel() + "_denom FUNC=x/y PERIODIC=NO");
  }
}

}
}
