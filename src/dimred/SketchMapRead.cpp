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
#include "core/ActionRegister.h"
#include "SketchMapBase.h"
#include "reference/ReferenceConfiguration.h"
#include "reference/MetricRegister.h"
#include "tools/ConjugateGradient.h"
#include "tools/GridSearch.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "tools/PDB.h"

//+PLUMEDOC DIMRED SKETCHMAP_READ
/*
Read in a sketch-map projection from an input file

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace dimred {

class SketchMapRead : public SketchMapBase {
private:
  std::vector<std::string> property;
  std::vector<ReferenceConfiguration*> myframes;
public:
  static void registerKeywords( Keywords& keys ); 
  SketchMapRead( const ActionOptions& ao );
  void minimise( Matrix<double>& );
  void getDataPoint( const unsigned& idata, std::vector<double>& point, double& weight ) const ;
  ReferenceConfiguration* getReferenceConfiguration( const unsigned& idata, const bool& calcdist );
  unsigned getNumberOfDataPoints() const ;
  unsigned getDataPointIndexInBase( const unsigned& idata ) const ;
  double getDissimilarity( const unsigned& i, const unsigned& j );
  double getWeight( const unsigned& idata ) const ;
};

PLUMED_REGISTER_ACTION(SketchMapRead,"SKETCHMAP_READ")

void SketchMapRead::registerKeywords( Keywords& keys ){
  SketchMapBase::registerKeywords( keys ); keys.remove("USE_OUTPUT_DATA_FROM");
  keys.add("compulsory","TYPE","OPTIMAL-FAST","the manner in which distances are calculated. More information on the different "
                                              "metrics that are available in PLUMED can be found in the section of the manual on "
                                              "\\ref dists");
  keys.add("compulsory","REFERENCE","the file containing the sketch-map projection");
  keys.add("compulsory","PROPERTY","the property to be used in the index. This should be in the REMARK of the reference");
}

SketchMapRead::SketchMapRead( const ActionOptions& ao ):
Action(ao),
SketchMapBase(ao)
{
  std::string mtype; parse("TYPE",mtype); 
  parseVector("PROPERTY",property); nlow=property.size();
  std::string ifilename; parse("REFERENCE",ifilename);
  FILE* fp=fopen(ifilename.c_str(),"r"); 
  if(!fp) error("could not open reference file " + ifilename ); 

  // Read in the embedding
  bool do_read=true; std::vector<double> weights; 
  unsigned nfram=0, wnorm=0., ww;
  while (do_read){
     PDB mypdb;
     // Read the pdb file
     do_read=mypdb.readFromFilepointer(fp,plumed.getAtoms().usingNaturalUnits(),0.1/atoms.getUnits().getLength());
     // Break if we are done
     if( !do_read ) break ;
     // Check for required properties
     if( !mypdb.hasRequiredProperties( property ) ) error("pdb input does not have contain required properties");
     // And read the frame
     myframes.push_back( metricRegister().create<ReferenceConfiguration>( mtype, mypdb ) );
     myframes[nfram]->checkRead(); ww=myframes[nfram]->getWeight();
     weights.push_back( ww );
     wnorm+=ww; nfram++;
  }
  fclose(fp); 

  if(nfram==0 ) error("no reference configurations were specified");
  log.printf("  found %u configurations in file %s\n",nfram,ifilename.c_str() );
  for(unsigned i=0;i<weights.size();++i) myframes[i]->setWeight( weights[i]/wnorm );
}

void SketchMapRead::minimise( Matrix<double>& projections ){
  for(unsigned i=0;i<myframes.size();++i){
      for(unsigned j=0;j<property.size();++j) projections(i,j) = myframes[i]->getPropertyValue( property[j] );
  }
}

void SketchMapRead::getDataPoint( const unsigned& idata, std::vector<double>& point, double& weight ) const {
  for(unsigned j=0;j<property.size();++j) point[j] = myframes[idata]->getPropertyValue( property[j] );
  weight = myframes[idata]->getWeight();
}

ReferenceConfiguration* SketchMapRead::getReferenceConfiguration( const unsigned& idata, const bool& calcdist ){
  return myframes[idata];
}

unsigned SketchMapRead::getNumberOfDataPoints() const {
  return myframes.size();
}

unsigned SketchMapRead::getDataPointIndexInBase( const unsigned& idata ) const {
  return idata;
}

// Highly unsatisfactory solution to problem GAT
double SketchMapRead::getDissimilarity( const unsigned& i, const unsigned& j ){
  plumed_merror("Cannot read in dissimilarities");
  return 0.0;
}

double SketchMapRead::getWeight( const unsigned& idata ) const {
  return myframes[idata]->getWeight();
}

}
}
