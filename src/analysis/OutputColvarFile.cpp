/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2019 The plumed team
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
#include "AnalysisBase.h"
#include "reference/ReferenceAtoms.h"
#include "reference/ReferenceArguments.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"

namespace PLMD {
namespace analysis {

//+PLUMEDOC ANALYSIS OUTPUT_ANALYSIS_DATA_TO_COLVAR
/*
This can be used to output the data that has been stored in an Analysis object.

The most useful application of this method is to output all projections of all the
points that were stored in an analysis object that performs some form of dimensionality
reduction.  If you use the USE_DIMRED_DATA_FROM option below projections of all the
stored points will be output to a file.  The positions of these projections will be calculated
using that dimensionality reduction algorithms out-of-sample extension algorithm.

\par Examples

*/
//+ENDPLUMEDOC

class OutputColvarFile : public AnalysisBase {
private:
  std::string fmt;
  std::string filename;
  bool output_for_all_replicas;
  std::vector<unsigned> preps;
  std::vector<std::string> req_vals;
public:
  static void registerKeywords( Keywords& keys );
  OutputColvarFile( const ActionOptions& );
  void performTask( const unsigned&, const unsigned&, MultiValue& ) const { plumed_error(); }
  void performAnalysis();
};

PLUMED_REGISTER_ACTION(OutputColvarFile,"OUTPUT_ANALYSIS_DATA_TO_COLVAR")

void OutputColvarFile::registerKeywords( Keywords& keys ) {
  AnalysisBase::registerKeywords( keys ); keys.use("ARG");
  keys.add("compulsory","FILE","the name of the file to output to");
  keys.add("compulsory","REPLICA","0","the replicas for which you would like to output this information");
  keys.add("compulsory","STRIDE","0","the frequency with which to perform the required analysis and to output the data.  The default value of 0 tells plumed to use all the data");
  keys.add("optional","FMT","the format to output the data using");
}

OutputColvarFile::OutputColvarFile( const ActionOptions& ao ):
  Action(ao),
  AnalysisBase(ao),
  fmt("%f"),
  output_for_all_replicas(false)
{
  parse("FILE",filename); parse("FMT",fmt);
  if( !getRestart() ) { OFile ofile; ofile.link(*this); ofile.setBackupString("analysis"); ofile.backupAllFiles(filename); }
  log.printf("  printing data to file named %s \n",filename.c_str() );
  if( getArguments().size()==0 ) {
    std::vector<std::string> tmp_vals( my_input_data->getArgumentNames() );
    req_vals.resize( tmp_vals.size() ); for(unsigned i=0; i<tmp_vals.size(); ++i) req_vals[i]=tmp_vals[i];
  } else {
    req_vals.resize( getArguments().size() ); for(unsigned i=0; i<req_vals.size(); ++i) req_vals[i]=getPntrToArgument(i)->getName();
  }
  if( req_vals.size()==0 ) {
    log.printf("  outputting weights from input action \n");
  } else {
    log.printf("  outputting %s", req_vals[0].c_str() );
    for(unsigned i=1; i<req_vals.size(); ++i) log.printf(",", req_vals[i].c_str() );
    log.printf("\n");
  }
  std::vector<std::string> rep_data; parseVector("REPLICA",rep_data);
  if( rep_data.size()==1 ) {
    if( rep_data[0]=="all" ) output_for_all_replicas=true;
    else {
      preps.resize(1); Tools::convert( rep_data[0], preps[0] );
    }
  } else {
    preps.resize( rep_data.size() );
    for(unsigned i=0; i<rep_data.size(); ++i) Tools::convert( rep_data[i], preps[i] );
  }
}

void OutputColvarFile::performAnalysis() {
  if( !output_for_all_replicas ) {
    bool found=false; unsigned myrep=plumed.multi_sim_comm.Get_rank();
    for(unsigned i=0; i<preps.size(); ++i) {
      if( myrep==preps[i] ) { found=true; break; }
    }
    if( !found ) return;
  }
  // Output the embedding as long lists of data
  OFile gfile; gfile.link(*this);
  gfile.setBackupString("analysis");
  gfile.fmtField(fmt+" ");
  gfile.open( filename.c_str() );

  // Print embedding coordinates
  for(unsigned i=0; i<getNumberOfDataPoints(); ++i) {
    DataCollectionObject& mydata=getStoredData(i, false);
    for(unsigned j=0; j<req_vals.size(); ++j) gfile.printField( req_vals[j], mydata.getArgumentValue(req_vals[j]) );
    gfile.printField( "weight", getWeight(i) ); gfile.printField();
  }
  gfile.close();
}

}
}
