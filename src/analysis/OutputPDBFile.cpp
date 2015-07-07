/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

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
#include "AnalysisWithDataCollection.h"
#include "reference/ReferenceAtoms.h"
#include "reference/ReferenceArguments.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"

namespace PLMD {
namespace analysis {

//+PLUMEDOC ANALYSIS OUTPUT_ANALYSIS_DATA_TO_PDB
/*
This can be used to output the data that has been stored in an Analysis object.

\par Examples

*/
//+ENDPLUMEDOC

class OutputPDBFile : public AnalysisWithDataCollection {
private:
  std::string fmt;
  std::string filename;
public:
  static void registerKeywords( Keywords& keys );
  OutputPDBFile( const ActionOptions& );
  void performTask( const unsigned& , const unsigned& , MultiValue& ) const { plumed_error(); }
  void performAnalysis();
};

PLUMED_REGISTER_ACTION(OutputPDBFile,"OUTPUT_ANALYSIS_DATA_TO_PDB")

void OutputPDBFile::registerKeywords( Keywords& keys ){
  AnalysisWithDataCollection::registerKeywords( keys );
  keys.add("compulsory","FILE","the name of the file to output to");
  keys.add("optional","FMT","the format to use in the output file");
  keys.reset_style("ATOMS","hidden"); keys.reset_style("STRIDE","hidden");
  keys.reset_style("RUN","hidden"); keys.reset_style("USE_ALL_DATA","hidden");
  keys.reset_style("REWEIGHT_BIAS","hidden"); keys.reset_style("REWEIGHT_TEMP","hidden");
  keys.reset_style("TEMP","hidden"); keys.reset_style("REUSE_INPUT_DATA_FROM","hidden");
  keys.reset_style("WRITE_CHECKPOINT","hidden"); keys.reset_style("NOMEMORY","hidden");
  keys.reset_style("RESTART","hidden"); keys.reset_style("UPDATE_FROM","hidden");
  keys.reset_style("UPDATE_UNTIL","hidden"); keys.reset_style("ARG","hidden");
}

OutputPDBFile::OutputPDBFile( const ActionOptions& ao ):
Action(ao),
AnalysisWithDataCollection(ao),
fmt("%f")
{
  parse("FILE",filename); parse("FMT",fmt);
  if( !getRestart() ){ OFile ofile; ofile.link(*this); ofile.setBackupString("analysis"); ofile.backupAllFiles(filename); }
  log.printf("  printing data to file named %s \n",filename.c_str() );
}

void OutputPDBFile::performAnalysis(){
  // Output the embedding in plumed pdb format
  OFile afile; afile.link(*this); afile.setBackupString("analysis"); std::size_t psign=fmt.find("%");
  afile.open( filename.c_str() ); std::string descr="REMARK WEIGHT=%-" + fmt.substr(psign+1) + " TYPE=" + getMetricName() + "\n";
  for(unsigned j=0;j<getNumberOfDataPoints();++j){
      afile.printf("DESCRIPTION: analysis data from calculation done at time %f \n",getLabel().c_str(),getTime() );
      afile.printf(descr.c_str(),getWeight(j) ); 
      if( plumed.getAtoms().usingNaturalUnits() ) getReferenceConfiguration(j,false)->print( 1.0, afile, fmt );
      else getReferenceConfiguration(j,false)->print( plumed.getAtoms().getUnits().getLength()/0.1, afile, fmt );
  }
  afile.close();
}

}
}
