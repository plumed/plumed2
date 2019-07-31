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
#include "tools/OFile.h"
#include "AnalysisBase.h"
#include "core/ActionRegister.h"

//+PLUMEDOC ANALYSIS PRINT_DISSIMILARITY_MATRIX
/*
Print the matrix of dissimilarities between a trajectory of atomic configurations.

\par Examples

*/
//+ENDPLUMEDOC


namespace PLMD {
namespace analysis {

class PrintDissimilarityMatrix : public AnalysisBase {
private:
  std::string fmt;
  std::string fname;
public:
  static void registerKeywords( Keywords& keys );
  explicit PrintDissimilarityMatrix( const ActionOptions& ao );
  void performAnalysis() override;
  void performTask( const unsigned&, const unsigned&, MultiValue& ) const override { plumed_error(); }
};

PLUMED_REGISTER_ACTION(PrintDissimilarityMatrix,"PRINT_DISSIMILARITY_MATRIX")

void PrintDissimilarityMatrix::registerKeywords( Keywords& keys ) {
  AnalysisBase::registerKeywords( keys );
  keys.add("compulsory","FILE","name of file on which to output the data");
  keys.add("optional","FMT","the format to use for the output of numbers");
  keys.add("compulsory","STRIDE","0","the frequency with which to perform the required analysis and to output the data.  The default value of 0 tells plumed to use all the data");
}

PrintDissimilarityMatrix::PrintDissimilarityMatrix( const ActionOptions& ao ):
  Action(ao),
  AnalysisBase(ao),
  fmt("%f")
{
  if( !dissimilaritiesWereSet() ) error("dissimilarities have not been set in base classes");

  parse("FILE",fname); parse("FMT",fmt);
  if( !getRestart() ) { OFile ofile; ofile.link(*this); ofile.setBackupString("analysis"); ofile.backupAllFiles(fname); }
  log.printf("  printing to file named %s with formt %s \n",fname.c_str(), fmt.c_str() );
}

void PrintDissimilarityMatrix::performAnalysis() {
  std::string ofmt=" "+fmt;
  OFile ofile; ofile.setBackupString("analysis"); ofile.open(fname);
  for(unsigned i=0; i<getNumberOfDataPoints(); ++i) {
    for(unsigned j=0; j<getNumberOfDataPoints(); ++j) ofile.printf(ofmt.c_str(), sqrt( my_input_data->getDissimilarity( i,j ) ) );
    ofile.printf("\n");
  }
  ofile.close();
}

}
}
