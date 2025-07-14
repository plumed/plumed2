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
#include "core/ActionWithValue.h"
#include "core/ActionWithArguments.h"
#include "core/ActionRegister.h"
#include "tools/File.h"

namespace PLMD {
namespace generic {

//+PLUMEDOC PRINTANALYSIS DUMPDERIVATIVES
/*
Dump the derivatives with respect to the input parameters for scalar values (generally CVs, functions or biases).

You can use this command to output the derivatives of a scalar with respect to whatever input parameters are used to
calculate that scalar.  For example, the following input instructs plumed to write a file called deriv that contains both the
analytical and numerical derivatives of the distance between atoms 1 and 2.

```plumed
distance: DISTANCE ATOMS=1,2
distanceN: DISTANCE ATOMS=1,2 NUMERICAL_DERIVATIVES
DUMPDERIVATIVES ARG=distance,distanceN STRIDE=1 FILE=deriv
```

This input outputs the derivative of the distance with respect to the atom positions and the cell vectors (virial-like form).
It will thus output 15 deriatives for the two quantities output - the derivatives with respect to the x,y and z components of the input atoms
and the 9 components of the virial.

By contast the following input will only output 1 derivative; namely, the derivative of the switching function, `f`, with respect to the input distance, `d`.

```plumed
d: DISTANCE ATOMS=1,2
f: LESS_THAN ARG=d SWITCH={RATIONAL R_0=0.2}
DUMPDERIVATIVES ARG=f STRIDE=1 FMT=%8.4f FILE=deriv
```

!!! warning ""

    You can only use this command to output the derivatives of rank 0 (scalar) values.
    You cannot use it to output derivatives of vector, matrix of function values.

## RESTART, UPDATE_FROM and UPDATE_UNTIL

Notice that the RESTART, UPDATE_FROM and UPDATE_UNTIL keywords keywords
can be used in this action in the same way as they are used for [PRINT](PRINT.md).
Consequently, if you would like to append derivatives to an existing file called `deriv` instead of backing that
file up at the start of the calculation and outputting the data from the calculation on a new file called `deriv`
you would use an input like the one shown below:

```plumed
distance: DISTANCE ATOMS=1,2
DUMPDERIVATIVES ARG=distance STRIDE=1 FILE=deriv RESTART=YES
```

Similarly, if you want to only output the derivatives during the 400 ps time interval after the first
100 ps of the simulation you would use an input like the one shown below:

```plumed
distance: DISTANCE ATOMS=1,2
DUMPDERIVATIVES ARG=distance STRIDE=1 FILE=deriv UPDATE_FROM=100 UPDATE_UNTIL=500
```

*/
//+ENDPLUMEDOC

class DumpDerivatives :
  public ActionPilot,
  public ActionWithArguments {
  std::string file;
  std::string fmt;
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
  keys.addInputKeyword("compulsory","ARG","scalar","the labels of the values whose derivatives should be output");
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
  fmt("%15.10f") {
  parse("FILE",file);
  if( file.length()==0 ) {
    error("name of output file was not specified");
  }
  parse("FMT",fmt);
  fmt=" "+fmt;
  of.link(*this);
  of.open(file);
  log.printf("  on file %s\n",file.c_str());
  log.printf("  with format %s\n",fmt.c_str());
  unsigned nargs=getNumberOfArguments();
  if( nargs==0 ) {
    error("no arguments specified");
  }
  (getPntrToArgument(0)->getPntrToAction())->turnOnDerivatives();
  if( getPntrToArgument(0)->getRank()>0 ) {
    error("cannot dump derivatives of non-scalar objects");
  }
  unsigned npar=getPntrToArgument(0)->getNumberOfDerivatives();
  if( npar==0 ) {
    error("one or more arguments has no derivatives");
  }
  for(unsigned i=1; i<nargs; i++) {
    (getPntrToArgument(i)->getPntrToAction())->turnOnDerivatives();
    if( getPntrToArgument(i)->getRank()>0 ) {
      error("cannot dump derivatives of non-scalar objects");
    }
    if( npar!=getPntrToArgument(i)->getNumberOfDerivatives() ) {
      error("the number of derivatives must be the same in all values being dumped");
    }
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
