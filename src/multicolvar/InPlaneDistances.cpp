/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2023 The plumed team
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

#include <string>
#include <cmath>

//+PLUMEDOC MCOLVAR INPLANEDISTANCES
/*
Calculate the distance between a pair of atoms in the plane

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace multicolvar {

class InPlaneDistances : public ActionShortcut {
public:
  explicit InPlaneDistances(const ActionOptions&);
  static void registerKeywords( Keywords& keys );
};

PLUMED_REGISTER_ACTION(InPlaneDistances,"INPLANEDISTANCES")

void InPlaneDistances::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys );
  keys.add("atoms","GROUP","calculate distance for each distinct set of three atoms in the group");
  keys.add("atoms","VECTORSTART","The first atom position that is used to define the normal to the plane of interest");
  keys.add("atoms","VECTOREND","The second atom position that is used to defin the normal to the plane of interest");
  keys.setValueDescription("vector","the INPLANEDISTANCE between each of the atoms specified using the GROUP keyword and atom A in the plane perpendicular to the vector connecting the atoms specified using VECTORSTART and VECTOREND");
  MultiColvarShortcuts::shortcutKeywords( keys );
  keys.needsAction("DISTANCE");
  keys.needsAction("ANGLE");
}

InPlaneDistances::InPlaneDistances(const ActionOptions&ao):
  Action(ao),
  ActionShortcut(ao) {
  std::vector<std::string> str_atomsA;
  parseVector("VECTORSTART",str_atomsA);
  Tools::interpretRanges( str_atomsA );
  std::vector<std::string> str_atomsB;
  parseVector("VECTOREND",str_atomsB);
  Tools::interpretRanges( str_atomsB );
  std::vector<std::string> str_atomsC;
  parseVector("GROUP",str_atomsC);
  Tools::interpretRanges( str_atomsC );
  unsigned n=1;
  std::string dinput= getShortcutLabel() + "_dis: DISTANCE", ainput = getShortcutLabel() + "_ang: ANGLE";
  for(unsigned i=0; i<str_atomsA.size(); ++i ) {
    for(unsigned j=0; j<str_atomsB.size(); ++j ) {
      for(unsigned k=0; k<str_atomsC.size(); ++k) {
        std::string str_n;
        Tools::convert( n, str_n );
        n++;
        dinput += " ATOMS" + str_n + "=" + str_atomsA[j] + "," + str_atomsC[k];
        ainput += " ATOMS" + str_n + "=" + str_atomsB[j] + "," + str_atomsA[i] + "," + str_atomsC[k];
      }
    }
  }
  readInputLine( dinput );
  readInputLine( ainput );
  readInputLine( getShortcutLabel() + ": CUSTOM PERIODIC=NO FUNC=x*sin(y) ARG=" + getShortcutLabel() + "_dis," + getShortcutLabel() + "_ang");
  MultiColvarShortcuts::expandFunctions( getShortcutLabel(), getShortcutLabel(), "", this );
}

}
}
