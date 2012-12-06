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
#include "tools/PlumedException.h"

using namespace std;

namespace PLMD{

//+PLUMEDOC GENERIC LOAD
/*
Loads a library, possibly defining new actions.

It is available only
on systems allowing for dynamic loading. It can also be fed with a cpp file,
in which case the file is compiled first.

\par Examples

\verbatim
LOAD FILE=extensions.so
\endverbatim

*/
//+ENDPLUMEDOC

class SetupLoad :
  public virtual ActionSetup
{
public:
  static void registerKeywords( Keywords& keys );
  SetupLoad(const ActionOptions&ao);
};

PLUMED_REGISTER_ACTION(SetupLoad,"LOAD")

void SetupLoad::registerKeywords( Keywords& keys ){
  ActionSetup::registerKeywords(keys);
  keys.add("compulsory","FILE","-","file to be loaded");
}

SetupLoad::SetupLoad(const ActionOptions&ao):
Action(ao),
ActionSetup(ao)
{
  std::string f;
  parse("FILE",f);
  checkRead();
  plumed_assert(f!="-");
  plumed.load(f);
}

}

