/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2023 The plumed team
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
#include "core/ActionWithArguments.h"
#include "core/ActionWithValue.h"
#include "core/ActionPilot.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "tools/Communicator.h"
#include "tools/OFile.h"

namespace PLMD {
namespace generic {

//+PLUMEDOC PRINTANALYSIS DUMPVECTOR
/*
Print a vector to a file

In the following input the four distances calculated by the [DISTANCE](DISTANCE.md) command
are output to files using this command and using [PRINT](PRINT.md)

```plumed
d: DISTANCE ATOMS1=1,2 ATOMS2=3,4 ATOMS3=5,6 ATOMS4=7,8
PRINT ARG=d FILE=colvar1 STRIDE=1
DUMPVECTOR ARG=d FILE=colvar2 STRIDE=1
```

The PRINT command outputs the instantaneous values of the four distances on a single row.  If we have a trajectory
with four frames the `colvar1` file will look like the one shown below:

````
! FIELDS time d.1 d.2 d.3 d.4
 0.000000 1.000000 1.000000 1.414214 1.000000
 1.000000 1.000000 1.000000 1.414214 1.000000
 2.000000 1.000000 1.000000 1.414214 1.000000
 3.000000 1.000000 1.000000 1.414214 1.000000
````

By contrast the DUMPVECTOR command will produce four output files - for each step of the simulation. Each of these
output files looks like this:

````
! FIELDS time parameter d
 3.000000 0 1.000000
 3.000000 1 1.000000
 3.000000 2 1.414214
 3.000000 3 1.000000
````

In other words, the four elements of the vector are printed on different rows of the output.
The latest file output will be called `colvar2` and all earlier files will be called `analysis.n.colvar` where
n is an integer that indicates the order in which files were produced. By adding the flag `PRINT_ONE_FILE` as shown below:

```plumed
d: DISTANCE ATOMS1=1,2 ATOMS2=3,4 ATOMS3=5,6 ATOMS4=7,8
DUMPVECTOR ARG=d FILE=colvar2 STRIDE=1 PRINT_ONE_FILE
```

We can ensure that all the data from all four frames is concatenated to a single file that looks as follows:

````
#! FIELDS time parameter d
 3.000000 0 1.000000
 3.000000 1 1.000000
 3.000000 2 1.414214
 3.000000 3 1.000000
#! FIELDS time parameter d
 0.000000 0 1.000000
 0.000000 1 1.000000
 0.000000 2 1.414214
 0.000000 3 1.000000
#! FIELDS time parameter d
 1.000000 0 1.000000
 1.000000 1 1.000000
 1.000000 2 1.414214
 1.000000 3 1.000000
#! FIELDS time parameter d
 2.000000 0 1.000000
 2.000000 1 1.000000
 2.000000 2 1.414214
 2.000000 3 1.000000
````

This command is useful for printing out time series that have been stored using the [COLLECT](COLLECT.md)
action or for printing out projections that have been generating using the tools in the [dimred](module_dimred.md)
module.  The following input shows how you can use it to calculate and print the time series of values for the
distances between atoms 1 and 2 and atoms 3 and 4.

```plumed
d1: DISTANCE ATOMS=1,2
d2: DISTANCE ATOMS=3,4
c1: COLLECT ARG=d1 STRIDE=1
c2: COLLECT ARG=d2 STRIDE=1
DUMPVECTOR ARG=c1,c2 FMT=%10.5f FILE=timeseries
```

In the example input above the time series is output at the end of the calculation.  Further notice how we have used
the FMT keyword here to control the appearance of the numbers in the output using a [C format specifier](https://www.simplilearn.com/tutorials/c-tutorial/format-specifiers-in-c)

## Outputing matrices

You can also use this command to output matrices as the following input demonstrates:

```plumed
d: DISTANCE_MATRIX GROUPA=1,2 GROUPB=3-6
DUMPVECTOR ARG=d FILE=matrix STRIDE=1
```

The files `matrix` and `analysis.n.matrix` that are output on each step here looks as follows:

````
! FIELDS time parameter d.1 d.2 d.3 d.4
 4.000000 0 2.000000 2.000000 1.000000 1.000000
 4.000000 1 1.000000 2.000000 2.000000 1.414214
````

In other words, the rows and columns of the file are used to display the rows and columns of the input matrix.

Further note that if your input matrix was constructed using the [VSTACK](VSTACK.md) and [COLLECT](COLLECT.md) commands
as shown below:

```plumed
d1: DISTANCE ATOMS=1,2
d2: DISTANCE ATOMS=3,4
d3: DISTANCE ATOMS=5,6
c1: COLLECT ARG=d1 STRIDE=1
c2: COLLECT ARG=d2 STRIDE=1
c3: COLLECT ARG=d3 STRIDE=1
v: VSTACK ARG=c1,c2,c3
DUMPVECTOR ARG=v FILE=matrix
```

The output file looks as follows:

````
#! FIELDS time parameter d1 d2 d3
 3.000000 0 1.000000 1.000000 1.414214
 3.000000 1 1.000000 1.000000 1.414214
 3.000000 2 1.000000 1.000000 1.414214
````

In other words, the names of the columns in the output reflect the names of the underlying variables that were
collected and stacked together to form the input matrix.

*/
//+ENDPLUMEDOC

class DumpVector :
  public ActionWithArguments,
  public ActionPilot {
private:
  bool onefile;
  std::vector<std::string> argnames;
  std::string fmt, filename;
  void buildArgnames();
public:
  static void registerKeywords( Keywords& keys );
  explicit DumpVector(const ActionOptions&ao);
  ~DumpVector() {}
  void calculate() override {}
  void apply() override {}
  void update() override ;
};

PLUMED_REGISTER_ACTION(DumpVector,"DUMPVECTOR")

void DumpVector::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionPilot::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys );
  keys.addInputKeyword("compulsory","ARG","vector/matrix","the labels of vectors/matrices that should be output in the file");
  keys.add("compulsory","STRIDE","0","the frequency with which the grid should be output to the file.");
  keys.add("compulsory","FILE","density","the file on which to write the vetors");
  keys.add("optional","FMT","the format that should be used to output real numbers");
  keys.addFlag("PRINT_ONE_FILE",false,"output vectors one after the other in a single file");
}

DumpVector::DumpVector(const ActionOptions&ao):
  Action(ao),
  ActionWithArguments(ao),
  ActionPilot(ao),
  fmt("%f") {
  if( getNumberOfArguments()==0 ) {
    error("found no arguments");
  }
  buildArgnames();
  parse("FILE",filename);
  parseFlag("PRINT_ONE_FILE", onefile);
  if(filename.length()==0) {
    error("name out output file was not specified");
  }

  log.printf("  outputting data with label %s to file named %s",getPntrToArgument(0)->getName().c_str(), filename.c_str() );
  parse("FMT",fmt);
  log.printf(" with format %s \n", fmt.c_str() );
  fmt = " " + fmt;
  if( onefile ) {
    log.printf("  printing all grids on a single file \n");
  } else {
    log.printf("  printing all grids on separate files \n");
  }
}

void DumpVector::buildArgnames() {
  argnames.resize(0);
  unsigned nvals = getPntrToArgument(0)->getShape()[0];
  if( getPntrToArgument(0)->getRank()==2 ) {
    nvals = getPntrToArgument(0)->getShape()[0];
  }
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    if( getPntrToArgument(i)->getShape()[0]!=nvals ) {
      error("all arguments should have same number of values");
    }
    if( getPntrToArgument(i)->getRank()==1 ) {
      argnames.push_back( getPntrToArgument(i)->getName() );
    } else if( getPntrToArgument(i)->getRank()==2 ) {
      (getPntrToArgument(i)->getPntrToAction())->getMatrixColumnTitles( argnames );
    }
  }
}

void DumpVector::update() {
  OFile ofile;
  ofile.link(*this);
  if( onefile ) {
    ofile.enforceRestart();
  } else {
    ofile.setBackupString("analysis");
  }
  ofile.open( filename );

  unsigned totargs = 0;
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    if( getPntrToArgument(i)->getRank()==1 ) {
      totargs += 1;
    } else if( getPntrToArgument(i)->getRank()==2 ) {
      totargs += getPntrToArgument(i)->getShape()[1];
    }
  }
  if( totargs!=argnames.size() ) {
    buildArgnames();
  }

  unsigned nvals = getPntrToArgument(0)->getShape()[0];
  for(unsigned i=0; i<nvals; ++i) {
    unsigned n=0;
    ofile.fmtField(" %f");
    ofile.printField("time",getTime());
    ofile.printField("parameter",int(i));
    for(unsigned j=0; j<getNumberOfArguments(); j++) {
      if( getPntrToArgument(j)->getRank()==1 ) {
        ofile.fmtField(fmt);
        ofile.printField(argnames[n],getPntrToArgument(j)->get(i) );
        n++;
      } else if( getPntrToArgument(j)->getRank()==2 ) {
        unsigned ncols = getPntrToArgument(j)->getShape()[1];
        for(unsigned k=0; k<ncols; ++k) {
          ofile.fmtField(fmt);
          ofile.printField(argnames[n],getPntrToArgument(j)->get(i*ncols+k));
          n++;
        }
      }
    }
    ofile.printField();
  }
}

}
}
