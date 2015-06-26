/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013,2014 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

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
#include "tools/OFile.h"
#include "core/ActionRegister.h"
#include "DissimilarityMatrixBase.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"

namespace PLMD {
namespace analysis {

class PrintDissimilarityMatrix : public ActionPilot {
private:
  std::string fmt;
  std::string fname;
  DissimilarityMatrixBase* mydissim;
  void printToFile( OFile& ofile );
public:
  static void registerKeywords( Keywords& keys );
  PrintDissimilarityMatrix( const ActionOptions& ao );
  void calculate(){};
  void apply(){};
  void update();
  void runFinalJobs();
};

PLUMED_REGISTER_ACTION(PrintDissimilarityMatrix,"PRINT_DISSIMILARITY_MATRIX")

void PrintDissimilarityMatrix::registerKeywords( Keywords& keys ){
  Action::registerKeywords( keys );
  ActionPilot::registerKeywords( keys );
  keys.add("compulsory","MATRIX","the action that calculates the matrix of dissimiliarities");
  keys.add("compulsory","FILE","name of file on which to output the data");
  keys.add("optional","FMT","the format to use for numbers in the output file");
}

PrintDissimilarityMatrix::PrintDissimilarityMatrix( const ActionOptions& ao ):
Action(ao),
ActionPilot(ao),
fmt("%f")
{
  std::string aname; parse("MATRIX",aname);
  mydissim = plumed.getActionSet().selectWithLabel<DissimilarityMatrixBase*>(aname);
  if(!mydissim) error("action labelled " +  aname + " does not exist or is not of Analysis type");
  addDependency( mydissim );

  if( mydissim->runFinalAnalysisOnly() ) setStride( mydissim->getStride() );

  log.printf("  printing dissimilarity matrix calculated by action with label %s \n",aname.c_str() );

  parse("FILE",fname); parse("FMT",fmt); fmt=" "+fmt;
  log.printf("  printing to file named %s with formt%s \n",fname.c_str(), fmt.c_str() );
}

void PrintDissimilarityMatrix::printToFile( OFile& ofile ){
  ofile.open(fname); 
  for(unsigned i=0;i<mydissim->getNumberOfDataPoints();++i){
      for(unsigned j=0;j<mydissim->getNumberOfDataPoints();++j) ofile.printf(fmt.c_str(), mydissim->getDissimilarity( i,j ) );
      ofile.printf("\n");
  }   
  ofile.close();
}

void PrintDissimilarityMatrix::update(){
  if( mydissim->runFinalAnalysisOnly() ) return ;
  OFile ofile; ofile.setBackupString("analysis"); printToFile( ofile );
}

void PrintDissimilarityMatrix::runFinalJobs(){
  OFile ofile; printToFile( ofile );
}

}
}
