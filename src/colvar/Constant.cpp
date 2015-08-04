/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2015 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

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
#include "core/PlumedMain.h"
#include "core/Atoms.h"

#include <string>
#include <cmath>

namespace PLMD{
namespace colvar{

//+PLUMEDOC COLVAR CONSTANT
/*
Return a constant quantity.

Useful in combination with functions.

\par Examples

The following input instructs plumed to compute the distance
between atoms 1 and 2. If this distance is between 1.0 and 2.0, it is
printed. If it is lower than 1.0 (larger than 2.0), 1.0 (2.0) is printed

\verbatim
one: CONSTANT VALUE=1.0
two: CONSTANT VALUE=2.0
dis: DISTANCE ATOMS=1,2
sss: SORT ARG=one,dis,two
PRINT ARG=sss.2
\endverbatim
(See also \ref DISTANCE, \ref SORT, and \ref PRINT).

*/
//+ENDPLUMEDOC


class Constant : public Colvar {
  double value;
public:
  explicit Constant(const ActionOptions&);
// active methods:
  virtual void calculate();
  static void registerKeywords( Keywords& keys );
};


using namespace std;


PLUMED_REGISTER_ACTION(Constant,"CONSTANT")

Constant::Constant(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao),
value(0.0)
{
  parse("VALUE",value);
  addValueWithDerivatives();
  setNotPeriodic();
// fake request to avoid errors:
  std::vector<AtomNumber> atoms;
  requestAtoms(atoms);
}

void Constant::registerKeywords( Keywords& keys ){
  Action::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  keys.remove("NUMERICAL_DERIVATIVES"); 
  keys.add("compulsory","VALUE","The value of the constant");
}

// calculator
void Constant::calculate(){
  setValue(value);
}

}
}



