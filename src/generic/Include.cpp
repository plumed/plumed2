/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2019 The plumed team
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
#include "core/ActionAnyorder.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "tools/Exception.h"

using namespace std;

namespace PLMD {
namespace generic {

//+PLUMEDOC GENERIC INCLUDE
/*
Includes an external input file, similar to "#include" in C preprocessor.

Useful to split very large plumed.dat files. Notice that in PLUMED 2.4 this action
cannot be used before the initial setup part of the file (e.g. in the part with \ref UNITS, \ref MOLINFO, etc).
As of PLUMED 2.5, \ref INCLUDE can be used in any position of the file.

\par Examples

This input:
\plumedfile
c1: COM ATOMS=1-100
c2: COM ATOMS=101-202
d: DISTANCE ATOMS=c1,c2
PRINT ARG=d
\endplumedfile
can be replaced with this input:
\plumedfile
INCLUDE FILE=pippo.dat
d: DISTANCE ATOMS=c1,c2
PRINT ARG=d
\endplumedfile
where the content of file pippo.dat is
\plumedfile
c1: COM ATOMS=1-100
c2: COM ATOMS=101-202
\endplumedfile

The files in this example are rather short, but imagine a case like this one:
\plumedfile
INCLUDE FILE=groups.dat
c: COORDINATION GROUPA=groupa GROUPB=groupb R_0=0.5
METAD ARG=c HEIGHT=0.2 PACE=100 SIGMA=0.2 BIASFACTOR=5
\endplumedfile
Here `groups.dat` could be huge file containing group definitions such as
\plumedfile
groupa: GROUP ...
  ATOMS={
    10
    50
    60
## imagine a long list here
    70
    80
    120
  }
...
groupb: GROUP ...
  ATOMS={
    11
    51
    61
## imagine a long list here
    71
    81
    121
  }
...
\endplumedfile
So, included files are the best place where one can store long definitions.

Another case where INCLUDE is very useful is when running multi-replica simulations.
Here different replicas might have different input files, but perhaps a large part of the
input is shared. This part can be put in a common included file. For instance you could have
`common.dat`:
\plumedfile
# this is common.dat
t: TORSION ATOMS=1,2,3,4
\endplumedfile
Then `plumed.0.dat`:
\plumedfile
# this is plumed.0.dat
INCLUDE FILE=common.dat
RESTRAINT ARG=t AT=1.0 KAPPA=10
\endplumedfile
And `plumed.1.dat`:
\plumedfile
# this is plumed.1.dat
INCLUDE FILE=common.dat
RESTRAINT ARG=t AT=1.2 KAPPA=10
\endplumedfile

\warning
Remember that when using multi replica simulations whenever plumed tried to open
a file for reading it looks for a file with the replica suffix first.
This is true also for files opened by INCLUDE!

As an example, the same result of the inputs above could have been obtained using
`plumed.dat`:
\plumedfile
# this is plumed.dat
t: TORSION ATOMS=1,2,3,4
INCLUDE FILE=other.dat
\endplumedfile
Then `other.0.dat`:
\plumedfile
# this is other.0.dat
RESTRAINT ARG=t AT=1.0 KAPPA=10
\endplumedfile
And `other.1.dat`:
\plumedfile
# this is other.1.dat
RESTRAINT ARG=t AT=1.2 KAPPA=10
\endplumedfile





*/
//+ENDPLUMEDOC

class Include :
  public ActionAnyorder
{
public:
  static void registerKeywords( Keywords& keys );
  explicit Include(const ActionOptions&ao);
  void calculate() {}
  void apply() {}
};

PLUMED_REGISTER_ACTION(Include,"INCLUDE")

void Include::registerKeywords( Keywords& keys ) {
  Action::registerKeywords(keys);
  keys.add("compulsory","FILE","file to be included");
}

Include::Include(const ActionOptions&ao):
  Action(ao),
  ActionAnyorder(ao)
{
  std::string f;
  parse("FILE",f);
  checkRead();
  plumed.readInputFile(f);
}

}
}

