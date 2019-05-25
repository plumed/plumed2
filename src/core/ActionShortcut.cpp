/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2018,2019 The plumed team
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

void ActionShortcut::readShortcutKeywords( const Keywords& keys, std::map<std::string,std::string>& keymap ) {
  for(unsigned i=0; i<keys.size(); ++i) {
      std::string t, keyname = keys.get(i);
      if( keys.style( keyname, "optional") || keys.style( keyname, "compulsory") ) {
          parse(keyname,t); 
          if( t.length()>0 ) {
             keymap.insert(std::pair<std::string,std::string>(keyname,t));
          } else if( keys.numbered( keyname ) ) {
             for(unsigned i=1;; ++i) {
               std::string istr; Tools::convert( i, istr );
               if( !parseNumbered(keyname,i,t) ) break ;
               keymap.insert(std::pair<std::string,std::string>(keyname + istr,t));
             }
          }
      } else if( keys.style( keyname, "flag") ) {
          bool found=false; parseFlag(keyname,found);
          if( found ) keymap.insert(std::pair<std::string,std::string>(keyname,""));
      } else plumed_merror("shortcut keywords should be optional, compulsory or flags");
  }
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
  for(auto p=line.begin(); p!=line.end(); ++p) {
      if( (*p).find(" " )!=std::string::npos ) {
          std::size_t eq = (*p).find_first_of("=");
          output += " " + (*p).substr(0,eq) + "={" + (*p).substr(eq+1) + "}";  
      } else output += " " + (*p);
  } 
  line.resize(0); return output;
}

}
