/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2023 The plumed team
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

//+PLUMEDOC ANALYSIS LOGSUMEXP
/*
This action takes the exponential of a vector of logarithms and divides each element of the vector by the sum of the exponentials.

The log-exp-sum trick is used here

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace landmarks {

class LogSumExp : public ActionShortcut {
private:
  std::string fixArgumentName( const std::string& argin );
public:
  static void registerKeywords( Keywords& keys );
  explicit LogSumExp( const ActionOptions& ao );
};

PLUMED_REGISTER_ACTION(LogSumExp,"LOGSUMEXP")

void LogSumExp::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys );
  keys.add("compulsory","ARG","the vector of logweights that you would like to normalise using the logsumexp trick");
  keys.setValueDescription("the logarithms of the input weights logweights that are computed with the log-sum weights formula");
  keys.needsAction("HIGHEST"); keys.needsAction("CUSTOM"); keys.needsAction("SUM");
}


LogSumExp::LogSumExp( const ActionOptions& ao ):
  Action(ao),
  ActionShortcut(ao)
{
  // Find the argument name
  std::string argn; parse("ARG",argn);
  // Find the maximum weight
  readInputLine( getShortcutLabel() + "_maxlogweight: HIGHEST ARG=" + argn );
  readInputLine( getShortcutLabel() + "_maxweight: CUSTOM ARG=" + getShortcutLabel() + "_maxlogweight FUNC=exp(x) PERIODIC=NO");
  // Calculate the maximum
  readInputLine( getShortcutLabel() + "_shiftw: CUSTOM ARG=" + argn + "," + getShortcutLabel() + "_maxlogweight FUNC=exp(x-y) PERIODIC=NO");
  // compute the sum of all the exponentials
  readInputLine( getShortcutLabel() + "_sumw: SUM ARG=" + getShortcutLabel() + "_shiftw PERIODIC=NO");
  // and the logsum
  readInputLine( getShortcutLabel() + "_logsum: CUSTOM ARG=" + getShortcutLabel() + "_sumw," + getShortcutLabel() + "_maxlogweight FUNC=y+log(x) PERIODIC=NO");
  // And the final weights
  readInputLine( getShortcutLabel() + ": CUSTOM ARG=" + argn + "," +  getShortcutLabel() + "_logsum FUNC=exp(x-y) PERIODIC=NO");
}

}
}
