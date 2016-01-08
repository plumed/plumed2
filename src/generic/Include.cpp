/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2016 The plumed team
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
#include "core/Action.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "tools/Exception.h"

using namespace std;

namespace PLMD{
namespace generic{

//+PLUMEDOC GENERIC INCLUDE
/*
Includes an external input file, similar to "#include" in C preprocessor.

Useful to split very large plumed.dat files.

\par Examples

This input
\verbatim
c1: COM ATOMS=1-100
c2: COM ATOMS=101-202
d: DISTANCE ARG=c1,c2
PRINT ARG=d
\endverbatim

can be replaced with
\verbatim
INCLUDE FILE=pippo.dat
d: DISTANCE ARG=c1,c2
PRINT ARG=d
\endverbatim

where the content of file pippo.dat is
\verbatim
c1: COM ATOMS=1-100
c2: COM ATOMS=101-202
\endverbatim

(see also \ref COM, \ref DISTANCE, and \ref PRINT).

*/
//+ENDPLUMEDOC

class Include :
  public Action
{
public:
  static void registerKeywords( Keywords& keys );
  explicit Include(const ActionOptions&ao);
  void calculate(){}
  void apply(){}
};

PLUMED_REGISTER_ACTION(Include,"INCLUDE")

void Include::registerKeywords( Keywords& keys ){
  Action::registerKeywords(keys);
  keys.add("compulsory","FILE","file to be included");
}

Include::Include(const ActionOptions&ao):
Action(ao)
{
  std::string f;
  parse("FILE",f);
  checkRead();
  plumed.readInputFile(f);
}

}
}

