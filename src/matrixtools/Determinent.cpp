/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2020 The plumed team
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

//+PLUMEDOC MCOLVAR DETERMINANT
/*
Calculate the determinant of a matrix

This shorctu allows you to calculate the [determinant](https://en.wikipedia.org/wiki/Determinant) of
a square matrix. The following example shows how the action is used:

```plumed
d1: DISTANCE_MATRIX ATOMS=1-5
det: DETERMINANT ARG=d1
PRINT ARG=det FILE=colvar
```

If you look at the expanded version of the input above you can see that PLUMED calculates the determinant
by first diagonalising the input matrix and then calculating the product of the eigenvalues.

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace matrixtools {

class Determinant : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  explicit Determinant(const ActionOptions&ao);
};

PLUMED_REGISTER_ACTION(Determinant,"DETERMINANT")

void Determinant::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords(keys);
  keys.addInputKeyword("compulsory","ARG","matrix","The matrix that we are calculating the determinant for");
  keys.setValueDescription("scalar","the determinant of the matrix");
  keys.needsAction("DIAGONALIZE");
  keys.needsAction("PRODUCT");
}

Determinant::Determinant( const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  std::string arg;
  parse("ARG",arg);
  // Compose a vector from the args
  readInputLine( getShortcutLabel() + "_diag: DIAGONALIZE ARG=" + arg + " VECTORS=all");
  // Not sure about the regexp here - check with matrix with more than 10 rows
  readInputLine( getShortcutLabel() + ": PRODUCT ARG=(" + getShortcutLabel() + R"=(_diag\.vals-[0-9]))=");
}

}
}
