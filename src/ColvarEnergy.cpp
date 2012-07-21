/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

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
#include "PlumedMain.h"

#include <string>
#include <cmath>
#include <cassert>

namespace PLMD{

//+PLUMEDOC COLVAR ENERGY
/*
Calculate the total energy of the simulation box.

\par Examples
The following input instructs plumed to print the energy of the system
\verbatim
ENERGY LABEL=ene
PRINT ARG=ene
\endverbatim
(See also \ref PRINT).

\bug This command does not work with variable cell simulations. Should it also include long range corrections?

*/
//+ENDPLUMEDOC


class ColvarEnergy : public Colvar {
  bool components;

public:
  ColvarEnergy(const ActionOptions&);
// active methods:
  virtual void calculate();
  static void registerKeywords( Keywords& keys );
};


using namespace std;


PLUMED_REGISTER_ACTION(ColvarEnergy,"ENERGY")

ColvarEnergy::ColvarEnergy(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao),
components(false)
{
  assert(!checkNumericalDerivatives());
  isEnergy=true;
  addValueWithDerivatives(); setNotPeriodic();
  getPntrToValue()->resizeDerivatives(1);
}

void ColvarEnergy::registerKeywords( Keywords& keys ){
  Action::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  keys.remove("NUMERICAL_DERIVATIVES"); 
}


// calculator
void ColvarEnergy::calculate(){
  setValue( getEnergy() );
  getPntrToComponent(0)->addDerivative(0,1.0);
}

}



