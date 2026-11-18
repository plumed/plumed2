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
#include "core/ActionShortcut.h"
#include "core/ActionRegister.h"
#include "tools/IFile.h"

namespace PLMD {
namespace generic {

//+PLUMEDOC COLVAR READCOLUMN
/*
Create a vector valued constant by reading one column from a colvar file.

This action creates a [CONSTANT](CONSTANT.md) action that outputs a vector quantity. The elements in this constant
vector are equal to the numbers in the column of an input colvar file. Thus, for example, the following input:

```plumed
#SETTINGS INPUTFILES=regtest/multicolvar/rt21/colvar.reference
r: READCOLUMN ARG=regtest/multicolvar/rt21/colvar.reference VALUE=d1_mean
```

will read the final column of the following file:

````
! FIELDS time d1_morethan d1_altmin d1_min d1_max d1_highest d1_lowest d1_mean
 0.000000  1.53877  1.14858  1.18493  8.15639  1.26259  1.18711  1.22485
 0.050000  1.55097  1.14309  1.16211  8.17171  1.31759  1.16228  1.23993
 0.100000  1.50216  1.08228  1.08684  8.17276  1.39339  1.08684  1.24011
 0.150000  1.47690  1.04339  1.04472  8.19389  1.47548  1.04472  1.26010
 0.200000  1.53834  1.08469  1.08643  8.22211  1.49076  1.08643  1.28859
````

and create a [CONSTANT](CONSTANT.md) vector `r` with five elements that can be used elsewhere in the PLUMED input.

*/
//+ENDPLUMEDOC

class ReadColumnFromColvar : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit ReadColumnFromColvar(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(ReadColumnFromColvar,"READCOLUMN")

void ReadColumnFromColvar::registerKeywords(Keywords& keys) {
  ActionShortcut::registerKeywords( keys );
  keys.add("compulsory","FILE","input file that contains the masses and charges that should be used");
  keys.add("compulsory","VALUE","the title of the column that you would like to read in to the output vector");
  keys.setValueDescription("vector","a vector that contains the values in the column of your input file");
  keys.needsAction("CONSTANT");
}

ReadColumnFromColvar::ReadColumnFromColvar(const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  std::string fname,vname,allvals,svalue,pstring="";
  parse("FILE",fname);
  parse("VALUE",vname);
  IFile mfile;
  mfile.open(fname);
  mfile.allowIgnoredFields();
  mfile.scanField(vname,allvals);
  if( mfile.FieldExist("min_" + vname ) ) {
    std::string smin, smax;
    mfile.scanField("min_" + vname,smin);
    mfile.scanField("max_" + vname,smax);
    pstring = " PERIODIC=" + smin + "," + smax;
  }
  mfile.scanField();
  while( mfile.scanField(vname,svalue)) {
    allvals += "," + svalue;
    mfile.scanField();
  }
  mfile.close();
  readInputLine( getShortcutLabel() + ": CONSTANT NOLOG VALUES=" + allvals + pstring);
}

}
}
