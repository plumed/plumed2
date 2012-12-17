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
#include "ActionWithValue.h"
#include "ActionRegister.h"
#include <string>
#include <cmath>

using namespace std;

namespace PLMD{

//+PLUMEDOC GENERIC TIME
/*
retrieve the time of the simulation to be used elsewere

\par Examples

\verbatim
TIME            LABEL=t1
PRINT ARG=t1
\endverbatim
(See also \ref PRINT).

*/
//+ENDPLUMEDOC
    
class GenericTime : public ActionWithValue {
public:
  static void registerKeywords( Keywords& keys );
  GenericTime(const ActionOptions&);
// active methods:
  virtual void calculate();
  virtual void apply(){};
};

PLUMED_REGISTER_ACTION(GenericTime,"TIME")

void GenericTime::registerKeywords( Keywords& keys ){
  Action::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
}

GenericTime::GenericTime(const ActionOptions&ao):
Action(ao),ActionWithValue(ao)
{
  addValueWithDerivatives(); setNotPeriodic();
  // resize derivative by hand to a nonzero value
  getPntrToValue()->resizeDerivatives(1);
}

void GenericTime::calculate(){
    setValue           (getTime());
}

}





