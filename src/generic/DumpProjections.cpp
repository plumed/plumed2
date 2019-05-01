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
#include "core/ActionPilot.h"
#include "core/ActionWithValue.h"
#include "core/ActionWithArguments.h"
#include "core/ActionRegister.h"
#include "tools/File.h"

using namespace std;

namespace PLMD {
namespace generic {

//+PLUMEDOC PRINTANALYSIS DUMPPROJECTIONS
/*
Dump the derivatives with respect to the input parameters for one or more objects (generally CVs, functions or biases).

\par Examples

Compute the distance between two groups and write on a file the
derivatives of this distance with respect to all the atoms of the two groups

\plumedfile
x1: CENTER ATOMS=1-10
x2: CENTER ATOMS=11-20
d: DISTANCE ATOMS=x1,x2
DUMPPROJECTIONS ARG=d FILE=proj STRIDE=20
\endplumedfile

*/
//+ENDPLUMEDOC

class DumpProjections :
  public ActionPilot,
  public ActionWithArguments
{
  string file;
  string fmt;
  OFile of;
public:
  void calculate() {}
  explicit DumpProjections(const ActionOptions&);
  static void registerKeywords(Keywords& keys);
  void apply() {}
  void update();
  bool checkNeedsGradients()const {return true;}
  ~DumpProjections();
};

PLUMED_REGISTER_ACTION(DumpProjections,"DUMPPROJECTIONS")

void DumpProjections::registerKeywords(Keywords& keys) {
  Action::registerKeywords(keys);
  ActionPilot::registerKeywords(keys);
  ActionWithArguments::registerKeywords(keys);
  keys.use("ARG");
  keys.add("compulsory","STRIDE","1","the frequency with which the derivatives should be output");
  keys.add("compulsory","FILE","the name of the file on which to output the derivatives");
  keys.add("compulsory","FMT","%15.10f","the format with which the derivatives should be output");
  keys.use("RESTART");
  keys.use("UPDATE_FROM");
  keys.use("UPDATE_UNTIL");
}

DumpProjections::DumpProjections(const ActionOptions&ao):
  Action(ao),
  ActionPilot(ao),
  ActionWithArguments(ao),
  fmt("%15.10f")
{
  parse("FILE",file);
  if( file.length()==0 ) error("filename not specified");
  parse("FMT",fmt);
  fmt=" "+fmt;
  of.open(file);
  log.printf("  on file %s\n",file.c_str());
  log.printf("  with format %s\n",fmt.c_str());
  checkRead();

  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    (getPntrToArgument(i)->getPntrToAction())->turnOnDerivatives();
  }
}


void DumpProjections::update() {
  of.fmtField(" %f").printField("time",getTime());
  for(unsigned i=0; i<getNumberOfArguments(); i++) {
    for(unsigned j=0; j<getNumberOfArguments(); j++) {
      of.fmtField(fmt);
      of.printField(getPntrToArgument(i)->getName()+"-"+getPntrToArgument(j)->getName(),getProjection(i,j));
    }
  }
  of.printField();
}

DumpProjections::~DumpProjections() {
}

}


}
