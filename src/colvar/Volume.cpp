/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2019 The plumed team
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
#include "Colvar.h"
#include "ActionRegister.h"

#include <string>
#include <cmath>

using namespace std;

namespace PLMD {
namespace colvar {

//+PLUMEDOC COLVAR VOLUME
/*
Calculate the volume of the simulation box.

\par Examples

The following input tells plumed to print the volume of the system
\plumedfile
vol: VOLUME
PRINT ARG=vol
\endplumedfile

*/
//+ENDPLUMEDOC


class Volume : public Colvar {

public:
  explicit Volume(const ActionOptions&);
// active methods:
  void calculate() override;
/// Register all the keywords for this action
  static void registerKeywords( Keywords& keys );
};

PLUMED_REGISTER_ACTION(Volume,"VOLUME")

Volume::Volume(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao)
{
  std::vector<AtomNumber> atoms;
  checkRead();

  addValueWithDerivatives(); setNotPeriodic();
  requestAtoms(atoms);
}

void Volume::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys );
}


// calculator
void Volume::calculate() {

  double v=getBox().determinant();
  setBoxDerivatives(-v*Tensor::identity());
  setValue         (v);
}

}
}



