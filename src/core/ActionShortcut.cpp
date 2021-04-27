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
#include "ActionRegister.h"
#include "ActionWithValue.h"
#include "PlumedMain.h"
#include "ActionSet.h"

namespace PLMD {

void ActionShortcut::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  keys.add("hidden","IS_SHORTCUT","hidden keyword to tell if actions are shortcuts so that example generator can provide expansions of shortcuts");
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

void ActionShortcut::readInputLine( const std::string& input, const bool never_update ) {
  std::string f_input = input; savedInputLines.push_back( input );
  if( !never_update ) {
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
  }
  plumed.readInputLine( f_input );
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

void ActionShortcut::interpretDataLabel( const std::string& mystr, Action* myuser, unsigned& nargs, std::vector<Value*>& arg ) const {
  std::size_t dot=mystr.find_first_of('.'); std::string a=mystr.substr(0,dot); std::string name=mystr.substr(dot+1);
  // Retrieve the keywords for the shortcut
  Keywords skeys; actionRegister().getKeywords( getName(), skeys );
  std::vector<std::string> out_comps( skeys.getAllOutputComponents() ); unsigned carg = nargs; 
  // Now get the output components
  if( name=="*" ) {
      for(unsigned k=0; k<out_comps.size(); ++k) {
          if( out_comps[k]=="" ) {
              ActionWithValue* action=plumed.getActionSet().selectWithLabel<ActionWithValue*>( a );
              if( action ) action->interpretDataLabel( a, myuser, nargs, arg );
          } else {
              ActionWithValue* action=plumed.getActionSet().selectWithLabel<ActionWithValue*>( a + "_" + out_comps[k] );
              if( action ) action->interpretDataLabel( a + "_" + out_comps[k], myuser, nargs, arg );
              else {
                 for(unsigned j=1;; ++j) {
                     std::string numstr; Tools::convert( j, numstr );
                     ActionWithValue* act=plumed.getActionSet().selectWithLabel<ActionWithValue*>( a + "_" + out_comps[k] + numstr ); 
                     if(!act) break;
                     act->interpretDataLabel( a + "_" + out_comps[k] + numstr, myuser, nargs, arg );
                 }
              }
          }
      }
  } else {
      for(unsigned k=0; k<out_comps.size(); ++k) {
          if( name.find(out_comps[k])!=std::string::npos ) {
              if( name==out_comps[k] ) {
                   ActionWithValue* action=plumed.getActionSet().selectWithLabel<ActionWithValue*>( a + "_" + name );
                   action->interpretDataLabel( a + "_" + name, myuser, nargs, arg );
              } else { 
                   for(unsigned j=1;; ++j) {
                       std::string numstr; Tools::convert( j, numstr );
                       if( name==out_comps[k] + numstr ) {
                           ActionWithValue* action=plumed.getActionSet().selectWithLabel<ActionWithValue*>( a + "_" + name + numstr );
                           action->interpretDataLabel( a + "_" + name + numstr, myuser, nargs, arg ); break;
                       }
                   }
              }
              break;
          }
      }
  }
}

std::vector<std::string> ActionShortcut::getSavedInputLines() const {
  return savedInputLines;
}

}
