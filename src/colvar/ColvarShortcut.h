/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2023 The plumed team
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
#ifndef __PLUMED_colvar_ColvarShortcut_h
#define __PLUMED_colvar_ColvarShortcut_h

#include "core/ActionShortcut.h"

namespace PLMD {
namespace colvar {

template <class T>
class ColvarShortcut : public ActionShortcut {
public:
  static void registerKeywords(Keywords&);
  explicit ColvarShortcut(const ActionOptions&);
};

template <class T>
void ColvarShortcut<T>::registerKeywords(Keywords& keys ) {
  T::registerKeywords( keys );
  keys.remove("NO_ACTION_LOG");
  unsigned nkeys = keys.size();
  for(unsigned i=0; i<nkeys; ++i) {
    if( keys.style( keys.get(i), "atoms" ) ) {
      keys.reset_style( keys.get(i), "numbered" );
    }
  }
  keys.addActionNameSuffix("_SCALAR");
  keys.addActionNameSuffix("_VECTOR");
}

template <class T>
ColvarShortcut<T>::ColvarShortcut(const ActionOptions&ao):
  Action(ao),
  ActionShortcut(ao) {
  bool scalar=true;
  unsigned nkeys = keywords.size();
  if( getName()=="MASS" || getName()=="CHARGE" || getName()=="POSITION" ) {
    std::string inpt;
    parse("ATOMS",inpt);
    if( inpt.length()>0 ) {
      readInputLine( getShortcutLabel() + ": " + getName() + "_VECTOR ATOMS=" + inpt + " " + convertInputLineToString() );
      scalar=false;
    }
  }
  for(unsigned i=0; i<nkeys; ++i) {
    if( keywords.style( keywords.get(i), "atoms" ) ) {
      std::string inpt;
      parseNumbered( keywords.get(i), 1, inpt );
      if( inpt.length()>0 ) {
        readInputLine( getShortcutLabel() + ": " + getName() + "_VECTOR " + keywords.get(i) + "1=" + inpt + " " + convertInputLineToString() );
        scalar=false;
        break;
      }
    }
  }
  if( scalar ) {
    readInputLine( getShortcutLabel() + ": " + getName() + "_SCALAR " + convertInputLineToString() );
  }
}

}
}
#endif
