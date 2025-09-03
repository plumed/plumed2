/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2023 The plumed team
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
#include "core/ActionAtomistic.h"

namespace PLMD {
namespace generic {

//+PLUMEDOC COLVAR ONES
/*
Create a constant vector with all elements equal to one

The following input creates and outputs a constant vector with 10 elements that are all equal to one

```plumed
ones: ONES SIZE=10
PRINT ARG=ones FILE=onesfile
```

Notice that the ONES action is a shortcut to [CONSTANT](CONSTANT.md).

This action is used extensively when calculating coordination numbers in inputs similar to this one:

```plumed
c1: CONTACT_MATRIX GROUP=1-7 SWITCH={RATIONAL R_0=2.6 NN=6 MM=12}
ones: ONES SIZE=7
cc: MATRIX_VECTOR_PRODUCT ARG=c1,ones
PRINT ARG=cc FILE=colvar
```

For more information on why this is useful see [CONTACT_MATRIX](CONTACT_MATRIX.md)

*/
//+ENDPLUMEDOC

class Ones : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit Ones(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(Ones,"ONES")

void Ones::registerKeywords(Keywords& keys) {
  ActionShortcut::registerKeywords( keys );
  keys.add("compulsory","SIZE","the number of ones that you would like to create");
  keys.setValueDescription("scalar/vector","a vector of ones with the required number of elements");
  keys.needsAction("CONSTANT");
}

Ones::Ones(const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  unsigned size;
  std::string matstr = "";
  std::vector<std::string> size_str;
  parseVector("SIZE",size_str);
  if( size_str.size()==2 ) {
    unsigned nr, nc;
    Tools::convert( size_str[0], nr );
    Tools::convert( size_str[1], nc );
    size = nr*nc;
    matstr = " NROWS=" + size_str[0] + " NCOLS=" + size_str[1];
  } else if( size_str[0]=="@natoms" ) {
    if( size_str.size()!=1 ) {
      error("should only be one @natoms string in input");
    }
    std::vector<Value*> xpos, ypos, zpos, masv, chargev;
    ActionAtomistic::getAtomValuesFromPlumedObject( plumed, xpos, ypos, zpos, masv, chargev );
    size = 0;
    for(unsigned i=0; i<xpos.size(); ++i ) {
      size += xpos[i]->getNumberOfValues();
    }
  } else {
    if( size_str.size()!=1 ) {
      error("size should be one or two values");
    }
    Tools::convert( size_str[0], size );
  }
  if( size<1 ) {
    error("size should be greater than 0");
  }
  std::string ones="1";
  for(unsigned i=1; i<size; ++i) {
    ones +=",1";
  }
  readInputLine( getShortcutLabel() + ": CONSTANT NOLOG VALUES=" + ones + matstr );
}

}
}
