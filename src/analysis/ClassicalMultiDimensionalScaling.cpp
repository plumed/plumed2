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
     Tools::convert(i+1,num);
     propnames[i]=getLabel() + "." + num;
  }
  myembedding->setPropertyNames( propnames, false );

  parse("EMBEDDING_OFILE",efilename);
  parse("OUTPUT_FILE",ofilename);
}

ClassicalMultiDimensionalScaling::~ClassicalMultiDimensionalScaling(){
  delete myembedding;
}

void ClassicalMultiDimensionalScaling::analyzeLandmarks(){
  printf("RUNNING ANALYSIS AT STEP %f\n",getStep() );

  // Output the embedding in plumed format
  if( efilename!="dont output"){
     //std::string ifname=saveResultsFromPreviousAnalyses( efilename );
     // OFile afile; afile.link(*this);
     // afile.open( efilename.c_str(), "w" );
     // myembedding->print( "classical mds", getTime(), afile );
     // afile.close();
  }
}

}
}
