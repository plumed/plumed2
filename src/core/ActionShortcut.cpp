/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2018-2020 The plumed team
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

void ActionShortcut::readInputLine( const std::string& input ) {
  std::string f_input = input;
  if( update_from!=std::numeric_limits<double>::max() ) {
    std::string ufrom; Tools::convert( update_from, ufrom ); f_input += " UPDATE_FROM=" + ufrom;
  }
  if( update_until!=std::numeric_limits<double>::max() ) {
    std::string util; Tools::convert( update_until, util ); f_input += " UPDATE_UNTIL=" + util;
  }
  if( keywords.exists("RESTART") ) {
    if( restart ) f_input += " RESTART=YES";
    if( !restart ) f_input += " RESTART=NO";
  }
  plumed.readInputLine( f_input );
}

const std::string & ActionShortcut::getShortcutLabel() const {
  return shortcutlabel;
}

}
