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
#include "ActionRegister.h"
#include "ActionPilot.h"
#include "PlumedMain.h"
#include "ActionSet.h"

namespace PLMD{

using namespace std;

//+PLUMEDOC GENERIC FLUSH
/*
This command instructs plumed to flush all the open files with a user specified frequency.  This
is useful for preventing data loss that would otherwise arrise as a consequence of the code
storing data for printing in the buffers

\par Examples
A command like this in the input will instruct plumed to flush all the output files every 100 steps
\verbatim
FLUSH STRIDE=100
\endverbatim
*/
//+ENDPLUMEDOC

class GenericFlush:
  public ActionPilot
{
public:
  GenericFlush(const ActionOptions&ao):
    Action(ao),
    ActionPilot(ao)
  {
    checkRead();
  }
  static void registerKeywords( Keywords& keys );
  void calculate(){};
  void apply(){
    plumed.fflush();
    log.flush();
    const ActionSet & actionSet(plumed.getActionSet());
    for(ActionSet::const_iterator p=actionSet.begin();p!=actionSet.end();++p)
    (*p)->fflush();
  }
};

PLUMED_REGISTER_ACTION(GenericFlush,"FLUSH")

void GenericFlush::registerKeywords( Keywords& keys ){
  Action::registerKeywords( keys );
  ActionPilot::registerKeywords( keys );
  keys.add("compulsory","STRIDE","the frequency with which all the open files should be flushed");
  keys.remove("LABEL");
}

}


