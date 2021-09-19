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
#ifndef __PLUMED_core_ActionShortcut_h
#define __PLUMED_core_ActionShortcut_h

#include "ActionWithArguments.h"

namespace PLMD {

/**
\ingroup MULTIINHERIT
Action used to create a command that expands to multiple PLMD::Action commands when read in during input
*/
class ActionShortcut :
  public virtual Action {
friend class ActionSet;
private:
  bool wildcard;
  std::string shortcutlabel;
  std::vector<std::string> savedInputLines;
public:
  static void registerKeywords( Keywords& keys );
/// Read keywords
  void readShortcutKeywords( const Keywords& keys, std::map<std::string,std::string>& keymap );
/// Constructor
  explicit ActionShortcut(const ActionOptions&ao);
/// Read a line of input and create appropriate actions
  void readInputLine( const std::string& input, const bool never_update=false );
/// Get the label for the shortcut
  const std::string & getShortcutLabel() const ;
/// Take everything that was input to this action and convert it to a string
  std::string convertInputLineToString();
/// This sorts out the reading of arguments from shortcuts
  void interpretDataLabel( const std::string& mystr, Action* myuser, std::vector<Value*>& args ) const ;
/// It is a shortcut it should never need to be activated
  void activate() {}
/// Do nothing.
  void calculate() override {}
/// Do nothing.
  void apply() override {}
/// Get the lines of the shortcut that were read in
  std::vector<std::string> getSavedInputLines() const ;
/// Is the output from this shortcut a match for the * keyword
  bool matchWildcard() const ;
};

}

#endif
