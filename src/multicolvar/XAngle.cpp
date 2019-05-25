/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2019 The plumed team
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

class XAngle : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit XAngle(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(XAngle,"XANGLES")
PLUMED_REGISTER_ACTION(XAngle,"YANGLES")
PLUMED_REGISTER_ACTION(XAngle,"ZANGLES")

void XAngle::registerKeywords(Keywords& keys) {
  ActionShortcut::registerKeywords( keys );
  keys.add("numbered","ATOMS","the pairs of atoms that you would like to calculate the angles for");
  keys.reset_style("ATOMS","atoms"); MultiColvarBase::shortcutKeywords( keys );
}

XAngle::XAngle(const ActionOptions& ao):
Action(ao),
ActionShortcut(ao)
{
  // Create distances
  std::string dline = getShortcutLabel() + ": DISTANCE COMPONENTS";
  for(unsigned i=1;;++i) {
      std::string atstring; parseNumbered("ATOMS",i,atstring);
      if( atstring.length()==0 ) break;
      std::string num; Tools::convert( i, num );
      dline += " ATOMS" + num + "=" + atstring;
  }
  readInputLine( dline );
  // Normalize the vectors
  readInputLine( getShortcutLabel() + "_norm: NORMALIZE ARG1=" + getShortcutLabel() + ".x ARG2=" + getShortcutLabel() + ".y ARG3=" + getShortcutLabel() + ".z"); 
  // Now compute the angles with matheval
  if( getName()=="XANGLES" ) readInputLine( getShortcutLabel() + "_ang: MATHEVAL FUNC=acos(x) PERIODIC=NO ARG1=" + getShortcutLabel() + "_norm.x");
  if( getName()=="YANGLES" ) readInputLine( getShortcutLabel() + "_ang: MATHEVAL FUNC=acos(x) PERIODIC=NO ARG1=" + getShortcutLabel() + "_norm.y");
  if( getName()=="ZANGLES" ) readInputLine( getShortcutLabel() + "_ang: MATHEVAL FUNC=acos(x) PERIODIC=NO ARG1=" + getShortcutLabel() + "_norm.z");   
  // Add shortcuts to label
  MultiColvarBase::expandFunctions( getShortcutLabel(), getShortcutLabel() + "_ang", "", this );
}

}
}
