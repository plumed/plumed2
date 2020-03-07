/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2020 The plumed team
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
#include "core/ActionRegister.h"
#include "core/ActionShortcut.h"

namespace PLMD {
namespace function {

class Product : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  explicit Product(const ActionOptions&ao);
};

PLUMED_REGISTER_ACTION(Product,"PRODUCT")

void Product::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords(keys);
  keys.add("compulsory","ARG","The point that we are calculating the distance from");
}

Product::Product( const ActionOptions& ao):
Action(ao),
ActionShortcut(ao)
{ 
  std::string arg; parse("ARG",arg); 
  readInputLine( getShortcutLabel() + "_vec: COMPOSE_VECTOR ARG=" + arg );
  //readInputLine("PRINT FILE=debug_1 ARG="+getShortcutLabel()+ "_vec STRIDE=1" ); //
  readInputLine( getShortcutLabel() + "_logs: MATHEVAL ARG1=" + getShortcutLabel() + "_vec FUNC=log(x) PERIODIC=NO");
  //readInputLine("PRINT FILE=debug_2 ARG="+getShortcutLabel()+ "_logs STRIDE=1" );//
  readInputLine( getShortcutLabel() + "_logsum: COMBINE ARG=" + getShortcutLabel() + "_logs PERIODIC=NO");
  //readInputLine("PRINT FILE=debug_3 ARG="+getShortcutLabel()+ "_logsum STRIDE=1" );//
  readInputLine( getShortcutLabel() + "_exp: MATHEVAL ARG1=" + getShortcutLabel() + "_logsum FUNC=exp(x) PERIODIC=NO");
  //readInputLine("PRINT FILE=debug_4 ARG="+getShortcutLabel()+ "_exp STRIDE=1" );//
}

}
}
