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
#include "core/ActionAtomistic.h"
#include "core/ActionRegister.h"
#include "tools/CheckInRange.h"

//+PLUMEDOC PRINTANALYSIS PRINT_NDX
/*
Print an ndx file

The following example shows how you can use this command to print out the indices of all the atoms
that have a coordination number that is greater than or equal to 4.

```plumed
# These three lines calculate the coordination numbers of 100 atoms
c1: CONTACT_MATRIX GROUP=1-100 SWITCH={RATIONAL R_0=0.1 NN=6 MM=12}
ones: ONES SIZE=100
cc: MATRIX_VECTOR_PRODUCT ARG=c1,ones

# This command then prints the indices of the atoms that have a coordination number that is greater than 4
# every step
PRINT_NDX ATOMS=1-100 ARG=cc FILE=index.ndx GREATER_THAN_OR_EQUAL=4
```

Obviously, if you want to print the indices of the atoms that have a coordination number that is equal to 4
you use a command like the one shown below:

```plumed
# These three lines calculate the coordination numbers of 100 atoms
c1: CONTACT_MATRIX GROUP=1-100 SWITCH={RATIONAL R_0=0.1 NN=6 MM=12}
ones: ONES SIZE=100
cc: MATRIX_VECTOR_PRODUCT ARG=c1,ones

# This command then prints the indices of the atoms that have a coordination number that is equal to 4
# on every 10th step of the simulation
PRINT_NDX ...
   ATOMS=1-100 ARG=cc FILE=index.ndx
   LESS_THAN_OR_EQUAL=4
   GREATER_THAN_OR_EQUAL=4
   STRIDE=10
...
```

A command like the one above is used in the [OUTPUT_CLUSTER](OUTPUT_CLUSTER.md) command that you can use to output the indices
of the atoms that are in a particular cluster.

## RESTART, UPDATE_FROM and UPDATE_UNTIL

Notice that the RESTART, UPDATE_FROM and UPDATE_UNTIL keywords keywords
can be used in this action in the same way as they are used for [PRINT](PRINT.md).
Consequently, if you would like to append to an existing file called `index.ndx` instead of backing that
file up at the start of the calculation and outputting the data from the calculation on a new file called `index.ndx`
you would use an input like the one shown below:

```plumed
c1: CONTACT_MATRIX GROUP=1-100 SWITCH={RATIONAL R_0=0.1 NN=6 MM=12}
ones: ONES SIZE=100
cc: MATRIX_VECTOR_PRODUCT ARG=c1,ones

PRINT_NDX ...
  ATOMS=1-100 ARG=cc FILE=index.ndx
  GREATER_THAN_OR_EQUAL=4 RESTART=YES
...
```

Similarly, if you want to only output the `index.ndx` file during the 400 ps time interval after the first
100 ps of the simulation you would use an input like the one shown below:

```plumed
c1: CONTACT_MATRIX GROUP=1-100 SWITCH={RATIONAL R_0=0.1 NN=6 MM=12}
ones: ONES SIZE=100
cc: MATRIX_VECTOR_PRODUCT ARG=c1,ones

PRINT_NDX ...
  ATOMS=1-100 ARG=cc FILE=index.ndx
  GREATER_THAN_OR_EQUAL=4
  UPDATE_FROM=100 UPDATE_UNTIL=500
...
```

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace generic {

class PrintNDX :
  public ActionPilot,
  public ActionAtomistic,
  public ActionWithArguments {
  std::string file;
  OFile ofile;
  CheckInRange bounds;
public:
  void calculate() override {}
  std::string writeInGraph() const override;
  explicit PrintNDX(const ActionOptions&);
  static void registerKeywords(Keywords& keys);
  bool actionHasForces() override {
    return false;
  }
  void calculateNumericalDerivatives( ActionWithValue* a=NULL ) override {
    plumed_error();
  }
  void lockRequests() override;
  void unlockRequests() override;
  void apply() override {}
  void update() override;
  ~PrintNDX() {}
};

PLUMED_REGISTER_ACTION(PrintNDX,"PRINT_NDX")

void PrintNDX::registerKeywords(Keywords& keys) {
  Action::registerKeywords(keys);
  ActionPilot::registerKeywords(keys);
  ActionAtomistic::registerKeywords( keys );
  ActionWithArguments::registerKeywords(keys);
  keys.addInputKeyword("optional","ARG","vector","the labels of vectors that should be used when printind the NDX file");
  keys.add("atoms","ATOMS","the list of atoms that have the corresponding arguments");
  keys.add("compulsory","STRIDE","1","the frequency with which the quantities of interest should be output");
  keys.add("optional","FILE","the name of the file on which to output these quantities");
  keys.add("optional","LESS_THAN_OR_EQUAL","when printing with arguments that are vectors only print components of vectors have a value less than or equal to this value");
  keys.add("optional","GREATER_THAN_OR_EQUAL","when printing with arguments that are vectors only print components of vectors have a value greater than or equal to this value");
  keys.use("RESTART");
  keys.use("UPDATE_FROM");
  keys.use("UPDATE_UNTIL");
}

PrintNDX::PrintNDX(const ActionOptions&ao):
  Action(ao),
  ActionPilot(ao),
  ActionAtomistic(ao),
  ActionWithArguments(ao) {
  ofile.link(*this);
  parse("FILE",file);
  if(file.length()>0) {
    ofile.open(file);
    log.printf("  on file %s\n",file.c_str());
  } else {
    log.printf("  on plumed log file\n");
    ofile.link(log);
  }
  std::vector<AtomNumber> all_atoms;
  parseAtomList("ATOMS",all_atoms);
  std::vector<std::string> argnames( getNumberOfArguments() );
  requestAtoms( all_atoms, false );
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    if( getPntrToArgument(i)->getRank()!=1 || getPntrToArgument(i)->hasDerivatives() ) {
      error("arguments for print ndx should be vector");
    }
    if( getPntrToArgument(i)->getShape()[0]!=all_atoms.size() ) {
      error("mismatch between number of arguments and number of input atoms");
    }
    argnames[i] = getPntrToArgument(i)->getName();
  }
  log.printf("  printing ndx file containing indices of atoms that have arguments in ranges prescribed below \n");
  log.printf("  full set of atom indices investigated are : ");
  for(unsigned int i=0; i<all_atoms.size(); ++i) {
    if ( (i+1) % 25 == 0 ) {
      log.printf("  \n");
    }
    log.printf("  %d", all_atoms[i].serial());
  }
  log.printf("\n");
  std::vector<std::string> str_upper, str_lower;
  std::string errors;
  parseVector("LESS_THAN_OR_EQUAL",str_upper);
  parseVector("GREATER_THAN_OR_EQUAL",str_lower);
  if( !bounds.setBounds( getNumberOfArguments(), str_lower, str_upper, errors ) ) {
    error( errors );
  }
  if( bounds.wereSet() ) {
    log.printf("  %s \n", bounds.report( argnames ).c_str() );
  }
  checkRead();
}

std::string PrintNDX::writeInGraph() const {
  return getName() + "\n" + "FILE=" + file;
}

void PrintNDX::lockRequests() {
  ActionWithArguments::lockRequests();
  ActionAtomistic::lockRequests();
}

void PrintNDX::unlockRequests() {
  ActionWithArguments::unlockRequests();
  ActionAtomistic::unlockRequests();
}

void PrintNDX::update() {
  unsigned n=0;
  std::vector<double> argvals( getNumberOfArguments() );
  ofile.printf("[ %s step %ld ] \n", getLabel().c_str(), getStep() );
  for(unsigned i=0; i<getNumberOfAtoms(); ++i) {
    for(unsigned j=0; j<getNumberOfArguments(); ++j) {
      argvals[j] = getPntrToArgument(j)->get(i);
    }
    if( bounds.check( argvals ) ) {
      ofile.printf("%6d", getAbsoluteIndexes()[i].serial() );
      n++;
      if( n%15==0 ) {
        ofile.printf("\n");
      }
    }
  }
  ofile.printf("\n");
}

}


}
