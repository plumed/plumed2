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
#include "ActionShortcut.h"
#include "PlumedMain.h"
#include "ActionSet.h"

namespace PLMD {

void ActionShortcut::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
}

ActionShortcut::ActionShortcut(const ActionOptions&ao):
  Action(ao),
  shortcutlabel(label)
{
  std::string s; Tools::convert(plumed.getActionSet().size(),s);
  if( shortcutlabel==("@" + s) ) {
    std::string t; Tools::convert(plumed.getActionSet().size()+1,t);
    shortcutlabel="@" + t;
  } else label = ("@" + s);
}

const std::string & ActionShortcut::getShortcutLabel() const {
  return shortcutlabel;
}

std::string ActionShortcut::convertInputLineToString() {
  std::string output;
  for(auto p=line.begin(); p!=line.end(); ++p) output += " " + (*p); 
  line.resize(0); return output;
}

}
