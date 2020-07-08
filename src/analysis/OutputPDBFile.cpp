/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2020 The plumed team
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
#include "core/ActionSet.h"
#include "core/Atoms.h"
#include "core/GenericMolInfo.h"
#include "tools/PDB.h"

namespace PLMD {
namespace analysis {

//+PLUMEDOC ANALYSIS OUTPUT_ANALYSIS_DATA_TO_PDB
/*
This can be used to output the data that has been stored in an Analysis object.

\par Examples

*/
//+ENDPLUMEDOC

class OutputPDBFile : public AnalysisBase {
private:
  PDB mypdb;
  std::string fmt;
  std::string filename;
public:
  static void registerKeywords( Keywords& keys );
  explicit OutputPDBFile( const ActionOptions& );
  void performTask( const unsigned&, const unsigned&, MultiValue& ) const override { plumed_error(); }
  void performAnalysis() override;
};

PLUMED_REGISTER_ACTION(OutputPDBFile,"OUTPUT_ANALYSIS_DATA_TO_PDB")

void OutputPDBFile::registerKeywords( Keywords& keys ) {
  AnalysisBase::registerKeywords( keys );
  keys.add("compulsory","FILE","the name of the file to output to");
  keys.add("optional","FMT","the format to use in the output file");
  keys.add("compulsory","STRIDE","0","the frequency with which to perform the required analysis and to output the data.  The default value of 0 tells plumed to use all the data");
}

OutputPDBFile::OutputPDBFile( const ActionOptions& ao ):
  Action(ao),
  AnalysisBase(ao),
  fmt("%f")
{
  // Get setup the pdb
  mypdb.setAtomNumbers( my_input_data->getAtomIndexes() );
  mypdb.setArgumentNames( my_input_data->getArgumentNames() );

  // Find a moldata object
  auto* moldat=plumed.getActionSet().selectLatest<GenericMolInfo*>(this);
  if( ! moldat ) warning("PDB output files do not have atom types unless you use MOLDATA");

  parse("FILE",filename); parse("FMT",fmt);
  if( !getRestart() ) { OFile ofile; ofile.link(*this); ofile.setBackupString("analysis"); ofile.backupAllFiles(filename); }
  log.printf("  printing data to file named %s \n",filename.c_str() );
}

void OutputPDBFile::performAnalysis() {
  // Find a moldata object
  auto* mymoldat=plumed.getActionSet().selectLatest<GenericMolInfo*>(this);
  // Output the embedding in plumed pdb format
  OFile afile; afile.link(*this); afile.setBackupString("analysis"); std::size_t psign=fmt.find("%");
  afile.open( filename.c_str() ); std::string descr="REMARK WEIGHT=%-" + fmt.substr(psign+1) + "\n";
  for(unsigned j=0; j<getNumberOfDataPoints(); ++j) {
    afile.printf("DESCRIPTION: analysis data from calculation done by %s at time %f \n",getLabel().c_str(),getTime() );
    if( dissimilaritiesWereSet() ) afile.printf("REMARK %s \n", getDissimilarityInstruction().c_str() );
    afile.printf(descr.c_str(),getWeight(j) ); getStoredData(j,false).transferDataToPDB( mypdb );
    if( plumed.getAtoms().usingNaturalUnits() ) mypdb.print( 1.0, mymoldat, afile, fmt );
    else mypdb.print( plumed.getAtoms().getUnits().getLength()/0.1, mymoldat, afile, fmt );
  }
  afile.close();
}

}
}
