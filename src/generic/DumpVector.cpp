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

\par Examples

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
  keys.use("ARG");
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
    getPntrToArgument(i)->buildDataStore();
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
