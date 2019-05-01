/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2017-2019 The plumed team
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
#include "core/ActionSet.h"
#include "core/PlumedMain.h"

using namespace std;

namespace PLMD {
namespace generic {

//+PLUMEDOC GENERIC ENDPLUMED
/*
Terminate plumed input.

Can be used to effectively comment out the rest of the input file.
It can be useful to quickly ignore part of a long input file. However,
one should keep in mind that when opening the file it might be difficult to
find where the commented out part begins. Regular comments (with `#`) are
usually easier to read. Notice that \ref VimSyntax "VIM syntax" should be able
to detect this command and properly mark the rest of the file as a comment,
although since vim doesn't parse the whole file it might fail in doing so for long
input files.

\par Examples

\plumedfile
d: DISTANCE ATOMS=1,10
PRINT ARG=d FILE=COLVAR STRIDE=10
ENDPLUMED
commands here are ignored
PRINT ARG=d FILE=COLVAR STRIDE=1
\endplumedfile

*/
//+ENDPLUMEDOC
class EndPlumed:
  public Action
{
public:
  explicit EndPlumed(const ActionOptions&ao);
/// Register all the relevant keywords for the action
  static void registerKeywords( Keywords& keys );
  void calculate() {}
  void apply() {}
};

PLUMED_REGISTER_ACTION(EndPlumed,"ENDPLUMED")

void EndPlumed::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
}

EndPlumed::EndPlumed(const ActionOptions&ao):
  Action(ao)
{
  checkRead();
  plumed.setEndPlumed();
}

}
}

