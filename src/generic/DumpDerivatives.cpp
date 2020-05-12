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
#include "core/ActionWithValue.h"
#include "core/ActionWithArguments.h"
#include "core/ActionRegister.h"
#include "tools/File.h"

using namespace std;

namespace PLMD {
namespace generic {

//+PLUMEDOC PRINTANALYSIS DUMPDERIVATIVES
/*
Dump the derivatives with respect to the input parameters for one or more objects (generally CVs, functions or biases).

For a CV this line in input instructs plumed to print the derivative of the CV with respect to the atom positions
and the cell vectors (virial-like form).  In contrast, for a function or bias the derivative with respect to the input "CVs"
will be output.  This command is most often used to test whether or not analytic derivatives have been implemented correctly.  This
can be done by outputting the derivatives calculated analytically and numerically.  You can control the buffering of output using the \ref FLUSH keyword.

\par Examples

The following input instructs plumed to write a file called deriv that contains both the
analytical and numerical derivatives of the distance between atoms 1 and 2.
\plumedfile
DISTANCE ATOMS=1,2 LABEL=distance
DISTANCE ATOMS=1,2 LABEL=distanceN NUMERICAL_DERIVATIVES
DUMPDERIVATIVES ARG=distance,distanceN STRIDE=1 FILE=deriv
\endplumedfile

(See also \ref DISTANCE)

*/
//+ENDPLUMEDOC

class DumpDerivatives :
  public ActionPilot,
  public ActionWithArguments
{
  string file;
  string fmt;
  OFile of;
public:
  void calculate() override {}
  explicit DumpDerivatives(const ActionOptions&);
  static void registerKeywords(Keywords& keys);
  void apply() override {}
  void update() override;
  ~DumpDerivatives();
};

PLUMED_REGISTER_ACTION(DumpDerivatives,"DUMPDERIVATIVES")

void DumpDerivatives::registerKeywords(Keywords& keys) {
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

DumpDerivatives::DumpDerivatives(const ActionOptions&ao):
  Action(ao),
  ActionPilot(ao),
  ActionWithArguments(ao),
  fmt("%15.10f")
{
  parse("FILE",file);
  if( file.length()==0 ) error("name of output file was not specified");
  parse("FMT",fmt);
  fmt=" "+fmt;
  of.link(*this);
  of.open(file);
  log.printf("  on file %s\n",file.c_str());
  log.printf("  with format %s\n",fmt.c_str());
  unsigned nargs=getNumberOfArguments();
  if( nargs==0 ) error("no arguments specified");
  (getPntrToArgument(0)->getPntrToAction())->turnOnDerivatives();
  unsigned npar=getPntrToArgument(0)->getNumberOfDerivatives();
  if( npar==0 ) error("one or more arguments has no derivatives");
  for(unsigned i=1; i<nargs; i++) {
    (getPntrToArgument(i)->getPntrToAction())->turnOnDerivatives();
    if( npar!=getPntrToArgument(i)->getNumberOfDerivatives() ) error("the number of derivatives must be the same in all values being dumped");
  }
  checkRead();
}


void DumpDerivatives::update() {
  unsigned npar=getPntrToArgument(0)->getNumberOfDerivatives();
  for(unsigned ipar=0; ipar<npar; ipar++) {
    of.fmtField(" %f");
    of.printField("time",getTime());
    of.printField("parameter",(int)ipar);
    for(unsigned i=0; i<getNumberOfArguments(); i++) {
      of.fmtField(fmt);
      of.printField(getPntrToArgument(i)->getName(),getPntrToArgument(i)->getDerivative(ipar) );
    }
    of.printField();
  }
}

DumpDerivatives::~DumpDerivatives() {
}

}


}
