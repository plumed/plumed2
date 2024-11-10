/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2023 The plumed team
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

This command can be used to output quantities to a file. The following demonstrates a very 
simple example:

```plumed
d: DISTANCE ATOMS=1,2
PRINT ARG=d FILE=colvar
```

This command outputs the distance between atom 1 and 2 to a file called colvar every step of 
trajectory. If you want to output the distance less frequently you can use the STRIDE keyword
as illustrated below:

```plumed
d: DISTANCE ATOMS=1,2
PRINT ARG=d FILE=colvar STRIDE=10
```

With the input above the distance is only output on every 10th step. Notice this distance will only
be calculated on every step as this quantity is not used on every other step.

You can also use the PRINT command to output objects that have a rank that is greater than zero.  
For example, the following input calculates a vector of three distances. These three distances 
are then output to the output colvar file every step.

```plumed
d: DISTANCE ATOMS1=1,2 ATOMS2=3,4 ATOMS3=5,6
PRINT ARG=d FILE=colvar 
``` 

Alterantively, we can output a matrix by using an input like the one shown below:

```plumed
d: DISTANCE_MATRIX GROUP=1-4 
PRINT ARG=d FILE=colvar
```

The output `colvar` file here will contain the 16 elements of the matrix on each line output. You can 
also use PRINT to output functions on grids. However, you are probably better using [DUMPGRID](DUMPGRID.md) or
[DUMPCUBE](DUMPCUBE.md) to output functions on grids. 

The PRINT directive can be used multiple times in the input so you can print files with different strides or print different quantities
to different files as illustrated in the example input below:

```plumed
# compute distance:
distance: DISTANCE ATOMS=2,5
# compute total energy (potential)
energy: ENERGY
# print distance on a file
PRINT ARG=distance          STRIDE=10   FILE=COLVAR
# print both variables on another file
PRINT ARG=distance,energy   STRIDE=1000 FILE=COLVAR_ALL
```

The input above instructs plumed to print the distance between atoms 3 and 5 on a file
called COLVAR every 10 steps, and the distance and total energy on a file called COLVAR_ALL
every 1000 steps.

You can control the buffering of output using the [FLUSH](FLUSH.md) keyword. Furthermore, whether the output file is 
appended or backed up depends on whetehr or not the [RESTART](RESTART.md) action is present on the file. If you add `RESTART=yes` on the 
line containing the PRINT directive the file will be appended also.

Printing happens in the so-called "update" phase. This implies that printing
is affected by the presence of [UPDATE_IF](UPDATE_IF.md) actions. One can thus decide to start
and stop printing at preassigned values of time using the `UPDATE_FROM` and `UPDATE_UNTIL` keywords.
Keep into account that even on steps when the action is not updated (and thus the file is not printed)
the argument will be activated. In other words, if you use `UPDATE_FROM` to start printing at a given time,
the collective variables this PRINT statement depends on will be computed also before that time.

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
  std::string writeInGraph() const override;
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
  keys.addInputKeyword("compulsory","ARG","scalar/vector/matrix","the labels of the values that you would like to print to the file");
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
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    ofile.setupPrintValue( getPntrToArgument(i) );
    getPntrToArgument(i)->buildDataStore(true);
  }
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

std::string Print::writeInGraph() const {
  return getName() + "\n" + "FILE=" + file;
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
    ofile.fmtField(fmt); getPntrToArgument(i)->print( ofile );
  }
  ofile.printField();
}

Print::~Print() {
}

}


}
