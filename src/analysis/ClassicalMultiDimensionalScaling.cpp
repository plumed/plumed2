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
#include "AnalysisWithLandmarks.h"
#include "ClassicalScaling.h"
#include "reference/PointWiseMapping.h"
#include "core/ActionRegister.h"

namespace PLMD {
namespace analysis {

class ClassicalMultiDimensionalScaling : public AnalysisWithLandmarks {
private:
  unsigned nlow;
  std::string ofilename;
  std::string efilename;
  PointWiseMapping* myembedding;
public:
  static void registerKeywords( Keywords& keys );
  ClassicalMultiDimensionalScaling( const ActionOptions& ao );
  ~ClassicalMultiDimensionalScaling();
  void analyzeLandmarks();
};

PLUMED_REGISTER_ACTION(ClassicalMultiDimensionalScaling,"CLASSICAL_MDS")

void ClassicalMultiDimensionalScaling::registerKeywords( Keywords& keys ){
  AnalysisWithLandmarks::registerKeywords( keys );
  keys.add("compulsory","NLOW_DIM","number of low-dimensional coordinates required");
  keys.add("compulsory","OUTPUT_FILE","file on which to output the final embedding coordinates");
  keys.add("compulsory","EMBEDDING_OFILE","dont output","file on which to output the embedding in plumed input format");
}

ClassicalMultiDimensionalScaling::ClassicalMultiDimensionalScaling( const ActionOptions& ao ):
Action(ao),
AnalysisWithLandmarks(ao)
{
  myembedding = new PointWiseMapping( getMetricName(), false );
  setDataToAnalyze( dynamic_cast<MultiReferenceBase*>(myembedding) );

  parse("NLOW_DIM",nlow);
  if( nlow<1 ) error("dimensionality of low dimensional space must be at least one");
  std::vector<std::string> propnames( nlow ); std::string num;
  for(unsigned i=0;i<propnames.size();++i){
     Tools::convert(i+1,num); std::string lab=getLabel();
     if(lab.find("@")!=std::string::npos) propnames[i]=getName() + "." + num;
     else propnames[i]=getLabel() + "." + num;
  }
  myembedding->setPropertyNames( propnames, false );

  parseOutputFile("EMBEDDING_OFILE",efilename);
  parseOutputFile("OUTPUT_FILE",ofilename);
}

ClassicalMultiDimensionalScaling::~ClassicalMultiDimensionalScaling(){
  delete myembedding;
}

void ClassicalMultiDimensionalScaling::analyzeLandmarks(){
  // Calculate all pairwise diatances
  myembedding->calculateAllDistances( getPbc(), getArguments(), comm, myembedding->modifyDmat(), true );

  // Run multidimensional scaling
  ClassicalScaling::run( myembedding );

  // Output the embedding as long lists of data
//  std::string gfname=saveResultsFromPreviousAnalyses( ofilename );
  OFile gfile; gfile.link(*this); 
  gfile.setBackupString("analysis");
  gfile.fmtField(getOutputFormat()+" ");
  gfile.open( ofilename.c_str() );
  
  // Print embedding coordinates
  for(unsigned i=0;i<myembedding->getNumberOfReferenceFrames();++i){
      for(unsigned j=0;j<nlow;++j){
          std::string num; Tools::convert(j+1,num);
          gfile.printField( getLabel() + "." + num , myembedding->getProjectionCoordinate(i,j) );
      }
      gfile.printField();
  }  
  gfile.close();

  // Output the embedding in plumed format
  if( efilename!="dont output"){
     OFile afile; afile.link(*this); afile.setBackupString("analysis");
     afile.open( efilename.c_str() );
     myembedding->print( "classical mds", getTime(), afile, getOutputFormat() );
     afile.close();
  }
}

}
}
