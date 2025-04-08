/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2023 The plumed team
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
#include "core/ActionRegister.h"
#include "core/ActionShortcut.h"

//+PLUMEDOC FUNCTION PRODUCT
/*
Calculate the product of the input quantities

This shortcut can be used to calculate the product of a collection of scalar inputs as illustrated by the following
simple example.

```plumed
d1: DISTANCE ATOMS=1,2
d2: DISTANCE ATOMS=3,4
p: PRODUCT ARG=d1,d2
```

This action is currently only used in the [DETERMINANT](DETERMINANT.md) shortcut. In that
action the determinant of the input square matrix is found by finding the eigenvalues of the
input matrix using [DIAGONALIZE](DIAGONALIZE.md). The determinant is then found by evaluating
the product of these eigenvalues.

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace function {

class Product : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  explicit Product(const ActionOptions&ao);
};

PLUMED_REGISTER_ACTION(Product,"PRODUCT")

void Product::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords(keys);
  keys.add("compulsory","ARG","The set of scalars that you would like to multiply together");
  keys.setValueDescription("scalar","the product of all the elements in the input vector");
  keys.needsAction("CONCATENATE");
  keys.needsAction("CUSTOM");
  keys.needsAction("SUM");
}

Product::Product( const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  std::string arg;
  parse("ARG",arg);
  readInputLine( getShortcutLabel() + "_vec: CONCATENATE ARG=" + arg );
  readInputLine( getShortcutLabel() + "_logs: CUSTOM ARG=" + getShortcutLabel() + "_vec FUNC=log(x) PERIODIC=NO");
  readInputLine( getShortcutLabel() + "_logsum: SUM ARG=" + getShortcutLabel() + "_logs PERIODIC=NO");
  readInputLine( getShortcutLabel() + ": CUSTOM ARG=" + getShortcutLabel() + "_logsum FUNC=exp(x) PERIODIC=NO");
}

}
}
