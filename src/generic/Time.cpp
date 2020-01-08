/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2020 The plumed team
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
#include "core/ActionWithValue.h"
#include "core/ActionRegister.h"
#include <string>
#include <cmath>

using namespace std;

namespace PLMD {
namespace generic {

//+PLUMEDOC GENERIC TIME
/*
retrieve the time of the simulation to be used elsewhere

\par Examples

\plumedfile
TIME            LABEL=t1
PRINT ARG=t1
\endplumedfile

*/
//+ENDPLUMEDOC

class Time : public ActionWithValue {
public:
  static void registerKeywords( Keywords& keys );
  explicit Time(const ActionOptions&);
// active methods:
  void calculate() override;
  void apply() override {}
  unsigned getNumberOfDerivatives() override { return 0; }
};

PLUMED_REGISTER_ACTION(Time,"TIME")

void Time::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
}

Time::Time(const ActionOptions&ao):
  Action(ao),ActionWithValue(ao)
{
  addValueWithDerivatives(); setNotPeriodic();
  // resize derivative by hand to a nonzero value
  getPntrToValue()->resizeDerivatives(1);
}

void Time::calculate() {
  setValue           (getTime());
}

}





}
