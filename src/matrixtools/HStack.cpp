/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2017 The plumed team
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

#include <complex>

using namespace std;

namespace PLMD {
namespace symfunc {

//+PLUMEDOC MCOLVAR HSTACK
/*

\par Examples


*/
//+ENDPLUMEDOC


class HStack : public ActionShortcut {
private: 
public:
  static void registerKeywords( Keywords& keys );
  explicit HStack(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(HStack,"HSTACK")

void HStack::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys );
}

HStack::HStack( const ActionOptions& ao):
Action(ao),
ActionShortcut(ao)
{
  // Create the vstack
  readInputLine( getShortcutLabel() + "T: VSTACK " + convertInputLineToString() );
  // And transpose that badboy
  readInputLine( getShortcutLabel() + ": TRANSPOSE ARG=" + getShortcutLabel() + "T");
}

}
}

