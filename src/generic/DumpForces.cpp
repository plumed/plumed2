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
#include "tools/File.h"

namespace PLMD {
namespace generic {

//+PLUMEDOC PRINTANALYSIS DUMPFORCES
/*
Dump the force acting on one of a values in a file.

Consider the following PLUMED input:

```plumed
r: DISTANCE ATOMS=1,2
V: RESTRAINT ARG=r AT=0.2 KAPPA=100
```

This input will ultimately apply forces on the $x$, $y$ and $z$ components of two
atoms as well as the 9 cell vectors.  The force on the $i$th of these 15 quantities
is given by:

$$
F_i = - \frac{\textrm{d}V}{\textrm{d}r}\frac{\partial r}{\partial x_i}
$$

where $x_i$ is the $x$, $y$ or $z$ component of the position of one of the two atoms or one of the cell
vectors.  If we modify the input above by adding the DUMPFORCES command as shown below:

```plumed
r: DISTANCE ATOMS=1,2
V: RESTRAINT ARG=r AT=0.2 KAPPA=100
DUMPFORCES ARG=r FILE=forces
```

we can then monitor the value of $-\frac{\textrm{d}V}{\textrm{d}r}\frac{\partial r}$ in the output file
`forces`.  As explained by the equation above to get the forces on the atom this 'input force' needs to
be multiplied by $\frac{\partial r}{\partial x_i}$.  To view the various components of the $\frac{\partial r}{\partial x_i}$
you would use the [DUMPDERIVATIVES](DUMPDERIVATIVES.md) command.  To control the buffering of output you use the
[FLUSH](FLUSH.md) command.

*/
//+ENDPLUMEDOC

class DumpForces :
  public ActionPilot,
  public ActionWithArguments {
  std::string file;
  std::string fmt;
  OFile of;
public:
  void calculate() override {}
  explicit DumpForces(const ActionOptions&);
  static void registerKeywords(Keywords& keys);
  void apply() override {}
  void update() override;
  ~DumpForces();
};

PLUMED_REGISTER_ACTION(DumpForces,"DUMPFORCES")

void DumpForces::registerKeywords(Keywords& keys) {
  Action::registerKeywords(keys);
  ActionPilot::registerKeywords(keys);
  ActionWithArguments::registerKeywords(keys);
  keys.addInputKeyword("compulsory","ARG",Keywords::argType::scalar,"the labels of the values whose forces should be output");
  keys.add("compulsory","STRIDE","1","the frequency with which the forces should be output");
  keys.add("compulsory","FILE","the name of the file on which to output the forces");
  keys.add("compulsory","FMT","%15.10f","the format with which the derivatives should be output");
  keys.use("RESTART");
  keys.use("UPDATE_FROM");
  keys.use("UPDATE_UNTIL");
}

DumpForces::DumpForces(const ActionOptions&ao):
  Action(ao),
  ActionPilot(ao),
  ActionWithArguments(ao),
  fmt("%15.10f") {
  parse("FILE",file);
  if( file.length()==0 ) {
    error("name of file was not specified");
  }
  parse("FMT",fmt);
  fmt=" "+fmt;
  of.link(*this);
  of.open(file);
  log.printf("  on file %s\n",file.c_str());
  log.printf("  with format %s\n",fmt.c_str());
  if( getNumberOfArguments()==0 ) {
    error("no arguments have been specified");
  }
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    if( getPntrToArgument(i)->getRank()>0 ) {
      error("can only use DUMPFORCES to output forces on rank zero objects");
    }
  }
  checkRead();
}


void DumpForces::update() {
  of.fmtField(" %f");
  of.printField("time",getTime());
  for(unsigned i=0; i<getNumberOfArguments(); i++) {
    of.fmtField(fmt);
    of.printField(getPntrToArgument(i)->getName(),getPntrToArgument(i)->getForce());
  }
  of.printField();
}

DumpForces::~DumpForces() {
}

}


}
