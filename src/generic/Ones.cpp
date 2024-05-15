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

namespace PLMD {
namespace generic {

//+PLUMEDOC COLVAR ONES
/*
Create a constant vector with all elements equal to one

\par Examples

*/
//+ENDPLUMEDOC

class Ones : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit Ones(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(Ones,"ONES")

void Ones::registerKeywords(Keywords& keys) {
  ActionShortcut::registerKeywords( keys );
  keys.add("compulsory","SIZE","the number of ones that you would like to create");
  keys.setValueDescription("a vector of ones with the required number of elements");
  keys.needsAction("CONSTANT");
}

Ones::Ones(const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao)
{
  unsigned size; parse("SIZE",size); if( size<1 ) error("size should be greater than 0");
  std::string ones="1"; for(unsigned i=1; i<size; ++i) ones +=",1";
  readInputLine( getShortcutLabel() + ": CONSTANT NOLOG VALUES=" + ones );
}

}
}
