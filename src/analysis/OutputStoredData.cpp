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
#include "Analysis.h"
#include "reference/ReferenceAtoms.h"
#include "reference/ReferenceArguments.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"

namespace PLMD {
namespace analysis {

//+PLUMEDOC ANALYSIS OUTPUT_ANALYSIS_DATA
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

class OutputStoredData : public Analysis {
private:
  std::string ofilename;
  std::string efilename;
public:
  static void registerKeywords( Keywords& keys );
  OutputStoredData( const ActionOptions& );
  void performTask(){ plumed_error(); }
  void performAnalysis();
};

PLUMED_REGISTER_ACTION(OutputStoredData,"OUTPUT_ANALYSIS_DATA")

void OutputStoredData::registerKeywords( Keywords& keys ){
  Analysis::registerKeywords( keys );
  keys.add("compulsory","OUTPUT_FILE","dont output","this file will contain the selected landmarks in pdb format");
  keys.add("compulsory","LIST_FILE","dont output","this file contains a list of landmark points.  This can only be used when the distances are calculated from "
                                                  "the values of arguments.  That is to say when the metric does not involve atom positions");
}

OutputStoredData::OutputStoredData( const ActionOptions& ao ):
Action(ao),
Analysis(ao)
{
  parseOutputFile("OUTPUT_FILE",efilename);
  parseOutputFile("LIST_FILE",ofilename);
}

void OutputStoredData::performAnalysis(){

  if( ofilename!="dont output"){
      // Output the embedding as long lists of data
      OFile gfile; gfile.link(*this); 
      gfile.setBackupString("analysis");
      gfile.fmtField(getOutputFormat()+" ");
      gfile.open( ofilename.c_str() );

      // Can't print out all landmark data if we have reference atom positions
      ReferenceAtoms* myat=dynamic_cast<ReferenceAtoms*>( getReferenceConfiguration(0) );
      plumed_assert( !myat );
      
      // Print embedding coordinates
      for(unsigned i=0;i<getNumberOfDataPoints();++i){
          ReferenceArguments* myref=dynamic_cast<ReferenceArguments*>( getReferenceConfiguration(i) );
          plumed_assert( myref );
          for(unsigned j=0;j<myref->getReferenceArguments().size();++j){
              gfile.printField( myref->getArgumentNames()[j], myref->getReferenceArgument(j) );
          }
          gfile.printField();
      }  
      gfile.close();
  }

  // Output the embedding in plumed format
  if( efilename!="dont output"){
     OFile afile; afile.link(*this); afile.setBackupString("analysis"); std::size_t psign=getOutputFormat().find("%");
     afile.open( efilename.c_str() ); std::string descr="REMARK WEIGHT=%-" + getOutputFormat().substr(psign+1);
     for(unsigned j=0;j<getNumberOfDataPoints();++j){
         afile.printf("DESCRIPTION: landmark configuration from %s performed at time %f",getLabel().c_str(),getTime() );
         afile.printf(descr.c_str(),getReferenceConfiguration(j)->getWeight() ); 
         if( plumed.getAtoms().usingNaturalUnits() ) getReferenceConfiguration(j)->print( 1.0, afile, getOutputFormat() );
         else getReferenceConfiguration(j)->print( plumed.getAtoms().getUnits().getLength()/0.1, afile, getOutputFormat() );
     }
     afile.close();
  }
}

}
}
