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

Output a scalar or vector that contains a subset of the elements in the input vector or matrix.
In the example below the value `s` is a scalar that contains the distance between atoms 3 and 4.

```plumed
d: DISTANCE ATOMS1=1,2 ATOMS2=3,4 ATOMS3=5,6 ATOMS4=7,8
s: SELECT_COMPONENTS ARG=d COMPONENTS=2
```

In this example the output z is a 2 dimensional vector containing the distances between atoms 3 and 4
and 7 and 8.

```plumed
d: DISTANCE ATOMS1=1,2 ATOMS2=3,4 ATOMS3=5,6 ATOMS4=7,8
s: SELECT_COMPONENTS ARG=d COMPONENTS=2,4
```

Lastly, in this example we calculate a matrix of distances.  The scalar `s` that is output from the
SELECT_COMPONENTS action contains the distance between atoms 1 and 3 and 3 and 4.

```plumed
d: DISTANCE_MATRIX GROUP=1-5
s: SELECT_COMPONENTS ARG=d COMPONENTS=1.3,3.4
```

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
  keys.addInputKeyword("compulsory","ARG","vector/matrix","the value from which we are selecting components");
  keys.add("compulsory","COMPONENTS","the components in the input value that you woul like to build a new vector from");
  keys.needsAction("FLATTEN");
  keys.needsAction("CONSTANT");
  keys.needsAction("SELECT_WITH_MASK");
  keys.setValueDescription("vector","a vector containing the selected components");
}

SelectComponents::SelectComponents(const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  std::vector<std::string> argn;
  parseVector("ARG",argn);
  std::vector<Value*> theargs;
  ActionWithArguments::interpretArgumentList( argn, plumed.getActionSet(), this, theargs );
  if( theargs.size()!=1 ) {
    error("should only be one argument input to this action");
  }
  // Create an array that will eventually hold the mask
  std::vector<double> mask( theargs[0]->getNumberOfValues(), 1 );
  std::vector<std::string> elements;
  parseVector("COMPONENTS",elements);
  if( theargs[0]->getRank()==1 ) {
    for(unsigned i=0; i<elements.size(); ++i) {
      unsigned sel;
      Tools::convert( elements[i], sel );
      mask[sel-1]=0;
    }
  } else if( theargs[0]->getRank()==2 ) {
    for(unsigned i=0; i<elements.size(); ++i) {
      std::size_t dot = elements[i].find_first_of(".");
      if( dot==std::string::npos ) {
        error("found no dot in specification of required matrix element");
      }
      std::string istr=elements[i].substr(0,dot), jstr=elements[i].substr(dot+1);
      unsigned ival, jval;
      Tools::convert( istr, ival );
      Tools::convert( jstr, jval );
      mask[(ival-1)*theargs[0]->getShape()[1] + jval - 1] = 0;
    }
    readInputLine( getShortcutLabel() + "_flat: FLATTEN ARG=" + theargs[0]->getName() );
  } else {
    error("input to this argument should be a vector/matrix");
  }
  // Now create the mask action
  std::string mask_str;
  Tools::convert( mask[0], mask_str );
  unsigned check_mask=mask[0];
  for(unsigned i=1; i<mask.size(); ++i) {
    std::string num;
    Tools::convert( mask[i], num );
    mask_str += "," + num;
    check_mask +=mask[i];
  }
  if( mask.size()-check_mask!=elements.size() ) {
    error("found repeated indexes in COMPONENTS");
  }
  readInputLine( getShortcutLabel() + "_mask: CONSTANT VALUES=" + mask_str );
  // And finally create the selector
  if( theargs[0]->getRank()==1 ) {
    readInputLine( getShortcutLabel() + ": SELECT_WITH_MASK ARG=" + theargs[0]->getName() + " MASK=" + getShortcutLabel() + "_mask");
  } else if( theargs[0]->getRank()==2 ) {
    readInputLine( getShortcutLabel() + ": SELECT_WITH_MASK ARG=" + getShortcutLabel() + "_flat MASK=" + getShortcutLabel() + "_mask");
  }
}

}
}
