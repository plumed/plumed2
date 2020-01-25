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
#include "core/ActionSetup.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "tools/Exception.h"

using namespace std;

namespace PLMD {
namespace setup {

//+PLUMEDOC GENERIC RESTART
/*
Activate restart.

This is a Setup directive and, as such, should appear
at the beginning of the input file. It influences the way
PLUMED treat files open for writing (see also \ref Files).

Notice that it is also possible to enable or disable restart on a per-action
basis using the RESTART keyword on a single action. In this case,
the keyword should be assigned a value. RESTART=AUTO means that global
settings are used, RESTART=YES or RESTART=NO respectively enable
and disable restart for that single action.

\attention
This directive can have also other side effects, e.g. on \ref METAD
and \ref PBMETAD and on some analysis action.

\par Examples

Using the following input:
\plumedfile
d: DISTANCE ATOMS=1,2
PRINT ARG=d FILE=out
\endplumedfile
a new 'out' file will be created. If an old one is on the way, it will be automatically backed up.

On the other hand, using the following input:
\plumedfile
RESTART
d: DISTANCE ATOMS=1,2
PRINT ARG=d FILE=out
\endplumedfile
the file 'out' will be appended.

In the following case, file out1 will be backed up and file out2 will be concatenated
\plumedfile
RESTART
d1: DISTANCE ATOMS=1,2
d2: DISTANCE ATOMS=1,2
PRINT ARG=d1 FILE=out1 RESTART=NO
PRINT ARG=d2 FILE=out2
\endplumedfile

In the following case, file out will backed up even if the MD code thinks that we
are restarting. Notice that not all the MD code send to PLUMED information about restarts.
If you are not sure, always put `RESTART` when you are restarting and nothing when you aren't
\plumedfile
RESTART NO
d1: DISTANCE ATOMS=1,2
PRINT ARG=d1 FILE=out1
\endplumedfile





*/
//+ENDPLUMEDOC

class Restart :
  public virtual ActionSetup
{
public:
  static void registerKeywords( Keywords& keys );
  explicit Restart(const ActionOptions&ao);
};

PLUMED_REGISTER_ACTION(Restart,"RESTART")

void Restart::registerKeywords( Keywords& keys ) {
  ActionSetup::registerKeywords(keys);
  keys.addFlag("NO",false,"switch off restart - can be used to override the behavior of the MD engine");
}

Restart::Restart(const ActionOptions&ao):
  Action(ao),
  ActionSetup(ao)
{
  bool no=false;
  parseFlag("NO",no);
  bool md=plumed.getRestart();
  log<<"  MD code "<<(md?"did":"didn't")<<" require restart\n";
  if(no) {
    if(md) log<<"  Switching off restart\n";
    plumed.setRestart(false);
    log<<"  Not restarting simulation: files will be backed up\n";
  } else {
    if(!md) log<<"  Switching on restart\n";
    plumed.setRestart(true);
    log<<"  Restarting simulation: files will be appended\n";
  }
}

}
}
