/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2020 The plumed team
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
#include "analysis/AnalysisBase.h"
#include "reference/ReferenceAtoms.h"
#include "reference/ReferenceArguments.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/Atoms.h"
#include "core/SetupMolInfo.h"
#include "tools/PDB.h"
#include "PCA.h"

namespace PLMD {
namespace dimred {

//+PLUMEDOC DIMRED OUTPUT_PCA_PROJECTION
/*
This is used to output the projection calculated by principle component analysis

\par Examples

*/
//+ENDPLUMEDOC

class OutputPCAProjection : public analysis::AnalysisBase {
private:
  PDB mypdb;
  PCA* mypca;
  std::string fmt;
  std::string filename;
public:
  static void registerKeywords( Keywords& keys );
  explicit OutputPCAProjection( const ActionOptions& );
  void performTask( const unsigned&, const unsigned&, MultiValue& ) const { plumed_error(); }
  void performAnalysis();
};

PLUMED_REGISTER_ACTION(OutputPCAProjection,"OUTPUT_PCA_PROJECTION")

void OutputPCAProjection::registerKeywords( Keywords& keys ) {
  analysis::AnalysisBase::registerKeywords( keys );
  keys.add("compulsory","FILE","the name of the file to output to");
  keys.add("optional","FMT","the format to use in the output file");
  keys.add("compulsory","STRIDE","0","the frequency with which to perform the required analysis and to output the data.  The default value of 0 tells plumed to use all the data");
}

OutputPCAProjection::OutputPCAProjection( const ActionOptions& ao ):
  Action(ao),
  analysis::AnalysisBase(ao),
  fmt("%f")
{
  // Setup the PCA object
  mypca = dynamic_cast<PCA*>( my_input_data );
  if( !mypca ) error("input must be a PCA object");

  // Get setup the pdb
  mypdb.setAtomNumbers( my_input_data->getAtomIndexes() );
  mypdb.setArgumentNames( (mypca->my_input_data)->getArgumentNames() );

  // Find a moldata object
  std::vector<SetupMolInfo*> moldat=plumed.getActionSet().select<SetupMolInfo*>();
  if( moldat.empty() ) warning("PDB output files do not have atom types unless you use MOLDATA");

  parse("FILE",filename); parse("FMT",fmt);
  if( !getRestart() ) { OFile ofile; ofile.link(*this); ofile.setBackupString("analysis"); ofile.backupAllFiles(filename); }
  log.printf("  printing data to file named %s \n",filename.c_str() );
}

void OutputPCAProjection::performAnalysis() {
  // Find a moldata object
  std::vector<SetupMolInfo*> moldat=plumed.getActionSet().select<SetupMolInfo*>();
  if( moldat.size()>1 ) error("you should only have one MOLINFO action in your input file");
  SetupMolInfo* mymoldat=NULL; if( moldat.size()==1 ) mymoldat=moldat[0];
  // Output the embedding in plumed pdb format
  OFile afile; afile.link(*this); afile.setBackupString("analysis");
  mypdb.setAtomPositions( (mypca->myref)->getReferencePositions() );
  for(unsigned j=0; j<mypca->getArguments().size(); ++j) mypdb.setArgumentValue( (mypca->getArguments()[j])->getName(), (mypca->myref)->getReferenceArgument(j) );
  // And output the first frame
  afile.open( filename.c_str() ); afile.printf("REMARK TYPE=%s \n", mypca->mtype.c_str() );
  if( plumed.getAtoms().usingNaturalUnits() ) mypdb.print( 1.0, mymoldat, afile, fmt );
  else mypdb.print( atoms.getUnits().getLength()/0.1, mymoldat, afile, fmt );
  // And now ouput the eigenvectors
  for(unsigned dim=0; dim<mypca->nlow; ++dim) {
    afile.printf("REMARK TYPE=DIRECTION \n");
    mypdb.setAtomPositions( mypca->directions[dim].getReferencePositions() );
    for(unsigned j=0; j<mypca->getArguments().size(); ++j) mypdb.setArgumentValue( (mypca->getArguments()[j])->getName(), mypca->directions[dim].getReferenceArgument(j) );
    if( plumed.getAtoms().usingNaturalUnits() ) mypdb.print( 1.0, mymoldat, afile, fmt );
    else mypdb.print( atoms.getUnits().getLength()/0.1, mymoldat, afile, fmt );
  }
  afile.close();
}

}
}
