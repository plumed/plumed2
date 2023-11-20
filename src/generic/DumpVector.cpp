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
#include "core/ActionPilot.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "tools/Communicator.h"
#include "tools/OFile.h"

namespace PLMD {
namespace generic {

//+PLUMEDOC PRINTANALYSIS DUMPVECTOR
/*

\par Examples

*/
//+ENDPLUMEDOC

class DumpVector : 
public ActionWithArguments,
public ActionPilot {
private:
  std::string fmt, filename;
  std::vector<unsigned> preps;
  bool output_for_all_replicas, onefile, xyzfile;
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
  ActionWithArguments::registerKeywords( keys ); keys.use("ARG");
  keys.add("compulsory","STRIDE","1","the frequency with which the grid should be output to the file.");
  keys.add("compulsory","FILE","density","the file on which to write the vetors");
  keys.add("compulsory","REPLICA","all","the replica for which you would like to output this information");
  keys.add("optional","FMT","the format that should be used to output real numbers");
}

DumpVector::DumpVector(const ActionOptions&ao):
  Action(ao),
  ActionWithArguments(ao),
  ActionPilot(ao),
  fmt("%f"),
  output_for_all_replicas(false)
{
  if( getNumberOfArguments()==0 ) error("found no arguments");
  unsigned nvals = getPntrToArgument(0)->getNumberOfValues();
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
      if( getPntrToArgument(i)->getNumberOfValues()!=nvals ) error("all arguments should have same number of values");
      getPntrToArgument(i)->buildDataStore();
  }
  parse("FILE",filename);
  if(filename.length()==0) error("name out output file was not specified");

  log.printf("  outputting data with label %s to file named %s",getPntrToArgument(0)->getName().c_str(), filename.c_str() );
  parse("FMT",fmt); log.printf(" with format %s \n", fmt.c_str() ); fmt = " " + fmt; 

  std::string rep_data; parse("REPLICA",rep_data);
  if( rep_data=="all" ) output_for_all_replicas=true;
  else { preps.resize(1); Tools::convert( rep_data, preps[0] ); }
  if( output_for_all_replicas ) log.printf("  outputting files for all replicas \n");
  else {
    log.printf("  outputting data for replicas ");
    for(unsigned i=0; i<preps.size(); ++i) log.printf("%d ", preps[i] );
  }

  checkRead();
}

void DumpVector::update() {
  if( !output_for_all_replicas ) {
      bool found=false; unsigned myrep=plumed.multi_sim_comm.Get_rank();
      for(unsigned i=0; i<preps.size(); ++i) {
        if( myrep==preps[i] ) { found=true; break; }
      }
      if( !found ) return;
  }
  OFile ofile; ofile.link(*this);
  ofile.enforceRestart();
  ofile.open( filename );

  unsigned nvals = getPntrToArgument(0)->getNumberOfValues();
  for(unsigned i=0; i<nvals; ++i) {
      ofile.fmtField(" %f"); 
      ofile.printField("time",getTime()); 
      ofile.printField("parameter",int(i));
      for(unsigned j=0; j<getNumberOfArguments(); j++) {
        ofile.fmtField(fmt);
        ofile.printField(getPntrToArgument(j)->getName(),getPntrToArgument(j)->get(i) );
      }
      ofile.printField();
  } 
}

}
}