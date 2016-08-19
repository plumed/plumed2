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

class OutputColvarFile : public AnalysisWithDataCollection {
private:
  std::string fmt;
  std::string filename;
public:
  static void registerKeywords( Keywords& keys );
  OutputColvarFile( const ActionOptions& );
  void performTask( const unsigned& , const unsigned& , MultiValue& ) const { plumed_error(); }
  void performAnalysis();
};

PLUMED_REGISTER_ACTION(OutputColvarFile,"OUTPUT_ANALYSIS_DATA_TO_COLVAR")

void OutputColvarFile::registerKeywords( Keywords& keys ){
  AnalysisWithDataCollection::registerKeywords( keys );
  keys.add("compulsory","FILE","the name of the file to output to");
  keys.add("optional","FMT","the format to output the data using");
  keys.reset_style("ATOMS","hidden"); keys.reset_style("STRIDE","hidden"); 
  keys.reset_style("RUN","hidden"); keys.reset_style("USE_ALL_DATA","hidden");
  keys.reset_style("REWEIGHT_BIAS","hidden"); keys.reset_style("REWEIGHT_TEMP","hidden");
  keys.reset_style("TEMP","hidden"); keys.reset_style("REUSE_INPUT_DATA_FROM","hidden");
  keys.reset_style("WRITE_CHECKPOINT","hidden"); keys.reset_style("NOMEMORY","hidden");
  keys.reset_style("RESTART","hidden"); keys.reset_style("UPDATE_FROM","hidden");
  keys.reset_style("UPDATE_UNTIL","hidden"); keys.reset_style("ARG","hidden");
} 

OutputColvarFile::OutputColvarFile( const ActionOptions& ao ):
Action(ao),
AnalysisWithDataCollection(ao),
fmt("%f")
{
  parse("FILE",filename); parse("FMT",fmt);
  if( !getRestart() ){ OFile ofile; ofile.link(*this); ofile.setBackupString("analysis"); ofile.backupAllFiles(filename); }
  log.printf("  printing data to file named %s \n",filename.c_str() );
}

void OutputColvarFile::performAnalysis(){

  // Output the embedding as long lists of data
  OFile gfile; gfile.link(*this); 
  gfile.setBackupString("analysis");
  gfile.fmtField(fmt+" ");
  gfile.open( filename.c_str() );

  ReferenceConfiguration* myp = getReferenceConfiguration(0,false);
  if( myp->getNumberOfProperties()==0 ){
      plumed_assert( !dynamic_cast<ReferenceAtoms*>( myp ) );
      for(unsigned j=0;j<myp->getReferenceArguments().size();++j) gfile.setupPrintValue( getArguments()[j] );
  }

  // Print embedding coordinates
  for(unsigned i=0;i<getNumberOfDataPoints();++i){
      ReferenceConfiguration* mypoint=getReferenceConfiguration(i,false);
      for(unsigned j=0;j<mypoint->getNumberOfProperties();++j){
          gfile.printField( mypoint->getPropertyName(j), mypoint->getPropertyValue(j) );
      }
      if( mypoint->getNumberOfProperties()==0 ){
          ReferenceArguments* myref=dynamic_cast<ReferenceArguments*>( mypoint );
          plumed_assert( myref );
          for(unsigned j=0;j<myref->getReferenceArguments().size();++j){
              gfile.printField( getArguments()[j], myref->getReferenceArgument(j) );
          }
      }
      gfile.printField( "weight", getWeight(i) );
      gfile.printField();
  }  
  gfile.close();
}

}
}
