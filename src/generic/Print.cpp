/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2020 The plumed team
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
#include "core/ActionPilot.h"
#include "core/ActionWithArguments.h"
#include "core/ActionRegister.h"

namespace PLMD {
namespace generic {

//+PLUMEDOC PRINTANALYSIS PRINT
/*
Print quantities to a file.

This directive can be used multiple times
in the input so you can print files with different strides or print different quantities
to different files.  You can control the buffering of output using the \subpage FLUSH keyword.
Output file is either appended or backed up depending on the presence of the \ref RESTART action.
A per-action `RESTART` keyword can be used as well.

Notice that printing happens in the so-called "update" phase. This implies that printing
is affected by the presence of \ref UPDATE_IF actions. In addition, one might decide to start
and stop printing at preassigned values of time using the `UPDATE_FROM` and `UPDATE_UNTIL` keywords.
Keep into account that even on steps when the action is not updated (and thus the file is not printed)
the argument will be activated. In other words, if you use `UPDATE_FROM` to start printing at a given time,
the collective variables this PRINT statement depends on will be computed also before that time.

\par Examples

The following input instructs plumed to print the distance between atoms 3 and 5 on a file
called COLVAR every 10 steps, and the distance and total energy on a file called COLVAR_ALL
every 1000 steps.
\plumedfile
# compute distance:
distance: DISTANCE ATOMS=2,5
# compute total energy (potential)
energy: ENERGY
# print distance on a file
PRINT ARG=distance          STRIDE=10   FILE=COLVAR
# print both variables on another file
PRINT ARG=distance,energy   STRIDE=1000 FILE=COLVAR_ALL
\endplumedfile

Notice that \ref DISTANCE and \ref ENERGY are computed respectively every 10 and 1000 steps, that is
only when required.

*/
//+ENDPLUMEDOC

class Print :
  public ActionPilot,
  public ActionWithArguments
{
  std::string file;
  OFile ofile;
  std::string fmt;
// small internal utility
/////////////////////////////////////////
// these are crazy things just for debug:
// they allow to change regularly the
// printed argument
  int rotate;
  int rotateCountdown;
  int rotateLast;
  std::vector<Value*> rotateArguments;
/////////////////////////////////////////
public:
  void calculate() override {}
  void prepare() override;
  explicit Print(const ActionOptions&);
  static void registerKeywords(Keywords& keys);
  void apply() override {}
  void update() override;
  ~Print();
};

PLUMED_REGISTER_ACTION(Print,"PRINT")

void Print::registerKeywords(Keywords& keys) {
  Action::registerKeywords(keys);
  ActionPilot::registerKeywords(keys);
  ActionWithArguments::registerKeywords(keys);
  keys.use("ARG");
  keys.add("compulsory","STRIDE","1","the frequency with which the quantities of interest should be output");
  keys.add("optional","FILE","the name of the file on which to output these quantities");
  keys.add("optional","FMT","the format that should be used to output real numbers");
  keys.add("hidden","_ROTATE","some funky thing implemented by GBussi");
  keys.use("RESTART");
  keys.use("UPDATE_FROM");
  keys.use("UPDATE_UNTIL");
}

Print::Print(const ActionOptions&ao):
  Action(ao),
  ActionPilot(ao),
  ActionWithArguments(ao),
  fmt("%f"),
  rotate(0)
{
  ofile.link(*this);
  parse("FILE",file);
  if(file.length()>0) {
    ofile.open(file);
    log.printf("  on file %s\n",file.c_str());
  } else {
    log.printf("  on plumed log file\n");
    ofile.link(log);
  }
  parse("FMT",fmt);
  fmt=" "+fmt;
  log.printf("  with format %s\n",fmt.c_str());
  for(unsigned i=0; i<getNumberOfArguments(); ++i) ofile.setupPrintValue( getPntrToArgument(i) );
/////////////////////////////////////////
// these are crazy things just for debug:
// they allow to change regularly the
// printed argument
  parse("_ROTATE",rotate);
  if(rotate>0) {
    rotateCountdown=rotate;
    for(unsigned i=0; i<getNumberOfArguments(); ++i) rotateArguments.push_back( getPntrToArgument(i) );
    std::vector<Value*> a(1,rotateArguments[0]);
    requestArguments(std::vector<Value*>(1,rotateArguments[0]));
    rotateLast=0;
  }
/////////////////////////////////////////
  checkRead();
}

void Print::prepare() {
/////////////////////////////////////////
// these are crazy things just for debug:
// they allow to change regularly the
// printed argument
  if(rotate>0) {
    rotateCountdown--;
    if(rotateCountdown==0) {
      rotateCountdown=rotate;
      rotateLast++;
      rotateLast%=rotateArguments.size();
      requestArguments(std::vector<Value*>(1,rotateArguments[rotateLast]));
    }
  }
/////////////////////////////////////////
}

void Print::update() {
  ofile.fmtField(" %f");
  ofile.printField("time",getTime());
  for(unsigned i=0; i<getNumberOfArguments(); i++) {
    ofile.fmtField(fmt);
    ofile.printField( getPntrToArgument(i), getArgument(i) );
  }
  ofile.printField();
}

Print::~Print() {
}

}


}
