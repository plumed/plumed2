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
#include "ActionSetup.h"
#include "ActionRegister.h"
#include "PlumedMain.h"
#include "Atoms.h"
#include "PlumedException.h"

using namespace std;

namespace PLMD{

//+PLUMEDOC GENERIC RESTART
/*
Activate restart.

This is a Setup directive and, as such, should appear
at the beginning of the input file.

\par Examples

Using the following input:
\verbatim
d: DISTANCE ATOMS=1,2
PRINT ARG=d FILE=out
\endverbatim
a new 'out' file will be created. If an old one is on the way, it will be automatically backed up.
On the other hand, using the following input:
\verbatim
RESTART
d: DISTANCE ATOMS=1,2
PRINT ARG=d FILE=out
\endverbatim
the file 'out' will be appended.
(See also \ref DISTANCE and \ref PRINT).

\attention
This directive can have also other side effects, e.g. on \ref METAD

*/
//+ENDPLUMEDOC

class SetupRestart :
  public virtual ActionSetup
{
public:
  static void registerKeywords( Keywords& keys );
  SetupRestart(const ActionOptions&ao);
};

PLUMED_REGISTER_ACTION(SetupRestart,"RESTART")

void SetupRestart::registerKeywords( Keywords& keys ){
  ActionSetup::registerKeywords(keys);
}

SetupRestart::SetupRestart(const ActionOptions&ao):
Action(ao),
ActionSetup(ao)
{
  plumed.restart=true;
  log<<"Restarting simulation: files will be appended\n";
}

}
