/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2023 The plumed team
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

//+PLUMEDOC PRINTANALYSIS DUMPPROJECTIONS
/*
Dump the derivatives with respect to the input parameters for one or more objects (generally CVs, functions or biases).

Consider the following PLUMED input:

```plumed
x1: CENTER ATOMS=1-10
x2: CENTER ATOMS=11-20
r: DISTANCE ATOMS=x1,x2
```

The derivatives of the distance `r` with respect to the atomic positions, $X_i$, are given by:

$$
\frac{\partial r}{\partial X_i} = \frac{\partial r}{\partial x_1}\frac{\partial r}{\partial x_2}\frac{\partial x_1}{\partial X_i}\frac{\partial x_2}{\partial X_i}
$$

For the input above there are 69 of these derivatives in total. These include the derivatives of $r$ with respect to the $x$, $y$ and $z$ components of the 10 atomic
positions that are used to calculate $x_1$, the derivatives of $r$ with respect to the $x$, $y$ and $z$ components of the 10 atomic
positions that are used to calculate $x_2$ and the 9 derivatives of $r$ with respect to the cell vectors.

If we add a DUMPROJECTIONS command to the above input as shown below:

```plumed
x1: CENTER ATOMS=1-10
x2: CENTER ATOMS=11-20
d: DISTANCE ATOMS=x1,x2
DUMPPROJECTIONS ARG=d FILE=proj STRIDE=20
```

PLUMED will output these 69 derivatives to the proj file.

!!! caution "only works for scalars"

    This action cannot be used if non rank zero objects are passed between actions that are used to calculate the quantity (`d` in the above examples)
    whose projections are being output using the DUMPPROJECTIONS command

## RESTART, UPDATE_FROM and UPDATE_UNTIL

Notice that the RESTART, UPDATE_FROM and UPDATE_UNTIL keywords keywords
can be used in this action in the same way as they are used for [PRINT](PRINT.md).
Consequently, if you would like to append derivatives to an existing file called `proj` instead of backing that
file up at the start of the calculation and outputting the data from the calculation on a new file called `proj`
you would use an input like the one shown below:

```plumed
x1: CENTER ATOMS=1-10
x2: CENTER ATOMS=11-20
d: DISTANCE ATOMS=x1,x2
DUMPPROJECTIONS ARG=d STRIDE=1 FILE=proj FMT=%10.6f RESTART=YES
```

Similarly, if you want to only output the derivatives during the 400 ps time interval after the first
100 ps of the simulation you would use an input like the one shown below:

```plumed
x1: CENTER ATOMS=1-10
x2: CENTER ATOMS=11-20
d: DISTANCE ATOMS=x1,x2
DUMPPROJECTIONS ARG=d STRIDE=1 FILE=d UPDATE_FROM=100 UPDATE_UNTIL=500
```


*/
//+ENDPLUMEDOC

class DumpProjections :
  public ActionPilot,
  public ActionWithArguments {
  std::string file;
  std::string fmt;
  OFile of;
public:
  void calculate() override {}
  explicit DumpProjections(const ActionOptions&);
  static void registerKeywords(Keywords& keys);
  void apply() override {}
  void update() override;
  bool checkNeedsGradients()const override {
    return true;
  }
  ~DumpProjections();
};

PLUMED_REGISTER_ACTION(DumpProjections,"DUMPPROJECTIONS")

void DumpProjections::registerKeywords(Keywords& keys) {
  Action::registerKeywords(keys);
  ActionPilot::registerKeywords(keys);
  ActionWithArguments::registerKeywords(keys);
  keys.addInputKeyword("compulsory","ARG","scalar","the labels of the values whose gradients should be outpu");
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
  fmt("%15.10f") {
  parse("FILE",file);
  if( file.length()==0 ) {
    error("filename not specified");
  }
  parse("FMT",fmt);
  fmt=" "+fmt;
  of.open(file);
  log.printf("  on file %s\n",file.c_str());
  log.printf("  with format %s\n",fmt.c_str());
  checkRead();

  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    if( getPntrToArgument(i)->getRank()!=0 ) {
      error("cannot use DUMPPROJECTIONS with actions that are not scalars");
    }
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
