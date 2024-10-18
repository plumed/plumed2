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

//+PLUMEDOC MCOLVAR UWALLS
/*
Add lower walls to a vector of quantities

Depracated action: use [UPPER_WALLS](UPPER_WALLS.md)

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace multicolvar {

class UWalls : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit UWalls(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(UWalls,"UWALLS")

void UWalls::registerKeywords(Keywords& keys) {
  ActionShortcut::registerKeywords( keys );
  keys.setDeprecated("UPPER_WALLS");
  keys.add("compulsory","DATA","the values you want to restrain");
  keys.add("compulsory","AT","the radius of the sphere");
  keys.add("compulsory","KAPPA","the force constant for the wall.  The k_i in the expression for a wall.");
  keys.add("compulsory","OFFSET","0.0","the offset for the start of the wall.  The o_i in the expression for a wall.");
  keys.add("compulsory","EXP","2.0","the powers for the walls.  The e_i in the expression for a wall.");
  keys.add("compulsory","EPS","1.0","the values for s_i in the expression for a wall");
  keys.add("atoms","CATOMS","all the angles between the bonds that radiate out from these central atom are computed");
  keys.add("atoms","GROUP","a list of angls between pairs of bonds connecting one of the atoms specified using the CATOM command and two of the atoms specified here are computed");
  keys.add("compulsory","SWITCH","the switching function specifies that only those bonds that have a length that is less than a certain threshold are considered");
  keys.addOutputComponent("bias","default","scalar","the instantaneous value of the bias potential");
  keys.addOutputComponent("force2","default","scalar","the instantaneous value of the squared force due to this bias potential");
  keys.needsAction("UPPER_WALLS");
}

UWalls::UWalls(const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  std::string data;
  parse("DATA",data);
  readInputLine( getShortcutLabel() + ": UPPER_WALLS ARG=" + data + " " + convertInputLineToString() );
}

}
}
