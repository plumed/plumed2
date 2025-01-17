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

//+PLUMEDOC MCOLVAR BRIDGE
/*
Calculate a matrix with elements equal to one if there is a bridging atom between the two atoms

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace adjmat {

class Bridge : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit Bridge(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(Bridge,"BRIDGE")

void Bridge::registerKeywords(Keywords& keys) {
  ActionShortcut::registerKeywords(keys);
  keys.add("atoms","GROUP","the atoms for which you would like to calculate the adjacency matrix");
  keys.add("atoms","GROUPA","");
  keys.add("atoms","GROUPB","");
  keys.add("atoms","BRIDGING_ATOMS","The list of atoms that can form the bridge between the two interesting parts "
           "of the structure.");
  keys.add("optional","SWITCH","The parameters of the two switchingfunction in the above formula");
  keys.add("optional","SWITCHA","The switchingfunction on the distance between bridging atoms and the atoms in "
           "group A");
  keys.add("optional","SWITCHB","The switchingfunction on the distance between the bridging atoms and the atoms in "
           "group B");
  keys.needsAction("BRIDGE_MATRIX");
  keys.needsAction("SUM");
  keys.setValueDescription("the number of bridging atoms between the two groups");
}

Bridge::Bridge(const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  // Need to read in switch
  std::string s_inp, sfinput;
  parse("SWITCH",sfinput);
  if( sfinput.length()>0 ) {
    s_inp += "SWITCH={" + sfinput +"} ";
  }
  std::string sfinputa;
  parse("SWITCHA",sfinputa);
  if( sfinputa.length()>0 ) {
    s_inp += "SWITCHA={" + sfinputa +"} ";
  }
  std::string sfinputb;
  parse("SWITCHB",sfinputb);
  if( sfinputb.length()>0 ) {
    s_inp += "SWITCHB={" + sfinputb +"} ";
  }
  // Create the matrix object
  readInputLine( getShortcutLabel() + "_mat: BRIDGE_MATRIX " + s_inp + convertInputLineToString() );
  // Add all the elements of the matrix together
  readInputLine( getShortcutLabel() + ": SUM ARG=" + getShortcutLabel() + "_mat PERIODIC=NO");
}

}
}
