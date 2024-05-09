/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2017 The plumed team
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
#include "core/ActionWithArguments.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"

//+PLUMEDOC PRINTANALYSIS SELECT_COMPONENTS
/*
Create a new value to hold a subset of the components that are in a vector or matrix

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace valtools {

class SelectComponents : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
/// Constructor
  explicit SelectComponents(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(SelectComponents,"SELECT_COMPONENTS")

void SelectComponents::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys );
  keys.add("compulsory","ARG","the argument we are using to build the shortcut");
  keys.add("compulsory","COMPONENTS","the components in the input value that you woul like to build a new vector from");
  keys.needsAction("FLATTEN"); keys.needsAction("CONSTANT"); keys.needsAction("SELECT_WITH_MASK");
  keys.setValueDescription("a vector containing the selected components");
}

SelectComponents::SelectComponents(const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao)
{
  std::vector<std::string> argn; parseVector("ARG",argn); std::vector<Value*> theargs;
  ActionWithArguments::interpretArgumentList( argn, plumed.getActionSet(), this, theargs );
  if( theargs.size()!=1 ) error("should only be one argument input to this action");
  // Create an array that will eventually hold the mask
  std::vector<double> mask( theargs[0]->getNumberOfValues(), 1 );
  std::vector<std::string> elements; parseVector("COMPONENTS",elements);
  if( theargs[0]->getRank()==1 ) {
    for(unsigned i=0; i<elements.size(); ++i) { unsigned sel; Tools::convert( elements[i], sel ); mask[sel-1]=0; }
  } else if( theargs[0]->getRank()==2 ) {
    for(unsigned i=0; i<elements.size(); ++i) {
      std::size_t dot = elements[i].find_first_of(".");
      if( dot==std::string::npos ) error("found no dot in specification of required matrix element");
      std::string istr=elements[i].substr(0,dot), jstr=elements[i].substr(dot+1);
      unsigned ival, jval; Tools::convert( istr, ival ); Tools::convert( jstr, jval );
      mask[(ival-1)*theargs[0]->getShape()[1] + jval - 1] = 0;
    }
    readInputLine( getShortcutLabel() + "_flat: FLATTEN ARG=" + theargs[0]->getName() );
  } else error("input to this argument should be a vector/matrix");
  // Now create the mask action
  std::string mask_str; Tools::convert( mask[0], mask_str ); unsigned check_mask=mask[0];
  for(unsigned i=1; i<mask.size(); ++i) { std::string num; Tools::convert( mask[i], num ); mask_str += "," + num; check_mask +=mask[i]; }
  if( mask.size()-check_mask!=elements.size() ) error("found repeated indexes in COMPONENTS");
  readInputLine( getShortcutLabel() + "_mask: CONSTANT VALUES=" + mask_str );
  // And finally create the selector
  if( theargs[0]->getRank()==1 ) readInputLine( getShortcutLabel() + ": SELECT_WITH_MASK ARG=" + theargs[0]->getName() + " MASK=" + getShortcutLabel() + "_mask");
  else if( theargs[0]->getRank()==2 ) readInputLine( getShortcutLabel() + ": SELECT_WITH_MASK ARG=" + getShortcutLabel() + "_flat MASK=" + getShortcutLabel() + "_mask");
}

}
}
