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
#include "MultiColvarBase.h"

namespace PLMD {
namespace multicolvar {

class CoordAngles : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit CoordAngles(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(CoordAngles,"COORD_ANGLES")

void CoordAngles::registerKeywords(Keywords& keys) {
  ActionShortcut::registerKeywords( keys );
  keys.add("atoms","CATOMS","all the angles between the bonds that radiate out from these central atom are computed");
  keys.add("atoms","GROUP","a list of angls between pairs of bonds connecting one of the atoms specified using the CATOM command and two of the atoms specified here are computed");
  keys.add("compulsory","SWITCH","the switching function specifies that only those bonds that have a length that is less than a certain threshold are considered");
  //keys.addFlag("MEAN",false,"calculate the mean of all the quantities.");
  keys.addOutputComponent("mean","MEAN","the mean of the colvars");
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
  for(unsigned i=0;i<catoms.size();++i) {
      for(unsigned j=0;j<group.size();++j) { std::string num; Tools::convert( k, num ); atlist += " ATOMS" + num + "=" + catoms[i] + "," + group[j]; k++; }
  }
  // Calculate the distances
  readInputLine( getShortcutLabel() + "_dd: DISTANCES" + atlist );
  // Transform with the switching function
  std::string switch_input; parse("SWITCH",switch_input);
  readInputLine( getShortcutLabel() + "_sw: LESS_THAN ARG1=" + getShortcutLabel() + "_dd SWITCH={" + switch_input +"}");
  // And compute some corrections
  readInputLine( getShortcutLabel() + "_sw2: CUSTOM ARG1=" + getShortcutLabel() + "_sw FUNC=x*x PERIODIC=NO");
  readInputLine( getShortcutLabel() + "_swdiag: COMBINE ARG=" + getShortcutLabel() + "_sw2 PERIODIC=NO");
  // Now get the normalised vectors
  readInputLine( getShortcutLabel() + "_comp: DISTANCES" + atlist + " COMPONENTS"); 
  readInputLine( getShortcutLabel() + "_norm: NORMALIZE ARG1=" + getShortcutLabel() + "_comp.x" + " ARG2=" + getShortcutLabel() + "_comp.y ARG3=" + getShortcutLabel() + "_comp.z");
  readInputLine( getShortcutLabel() + "_stack: VSTACK ARG1=" + getShortcutLabel() + "_norm.x" + " ARG2=" + getShortcutLabel() + "_norm.y ARG3=" + getShortcutLabel() + "_norm.z"); 
  readInputLine( getShortcutLabel() + "_stackT: TRANSPOSE ARG=" + getShortcutLabel() + "_stack");
  // Create the matrix of weights
  readInputLine( getShortcutLabel() + "_swd: DOT ARG1=" + getShortcutLabel() + "_sw ARG2=" + getShortcutLabel() + "_sw");
  // And the matrix of dot products and the angles
  readInputLine( getShortcutLabel() + "_dpmat: DOT ARG1=" + getShortcutLabel() + "_stack ARG2=" + getShortcutLabel() + "_stackT");
  readInputLine( getShortcutLabel() + "_ang: CUSTOM ARG1=" + getShortcutLabel() + "_dpmat FUNC=acos(x) PERIODIC=NO");
  // Now the weights matrix
  readInputLine( getShortcutLabel() + "_wang: CUSTOM ARG1=" + getShortcutLabel() + "_swd ARG2=" + getShortcutLabel() + "_ang FUNC=x*y PERIODIC=NO");
  // And the numerator for the average
  readInputLine( getShortcutLabel() + "_ncoo: COORDINATIONNUMBER WEIGHT=" + getShortcutLabel() + "_wang");
  readInputLine( getShortcutLabel() + "_numer: COMBINE ARG=" + getShortcutLabel() + "_ncoo PERIODIC=NO");
  // Get the denominator
  readInputLine( getShortcutLabel() + "_denom: COMBINE ARG=" + getShortcutLabel() + "_swd PERIODIC=NO");
  // And the mean
  readInputLine( getShortcutLabel() + "_mean: CUSTOM ARG1=" + getShortcutLabel() + "_numer ARG2=" + getShortcutLabel() + "_denom ARG3=" + getShortcutLabel() + "_swdiag FUNC=x/(y-z) PERIODIC=NO");
}

}
}
