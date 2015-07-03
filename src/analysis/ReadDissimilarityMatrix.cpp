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
#include "AnalysisBase.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/ActionRegister.h"
#include "tools/IFile.h"

//+PLUMEDOC ANALYSIS READ_DISSIMILARITY_MATRIX
/* 
Read a matrix of dissimilarities between a trajectory of atomic configurations from a file.

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace analysis {

class ReadDissimilarityMatrix : public AnalysisBase {
private:
  unsigned nnodes;
  std::string fname, wfile;
  Matrix<double> dissimilarities;
  std::vector<double> weights;
public:
  static void registerKeywords( Keywords& keys );
  ReadDissimilarityMatrix( const ActionOptions& ao );
  unsigned getNumberOfDataPoints() const { return nnodes; }
/// This gives an error as if we read in the matrix we dont have the coordinates
  ReferenceConfiguration* getReferenceConfiguration( const unsigned& idata );
  ReferenceConfiguration* getInputReferenceConfiguration( const unsigned& idata );
/// This gives an error as if we read in the matrix we dont have the coordinates
  void getDataPoint( const unsigned& idata, std::vector<double>& point, double& weight ) const ;
/// Tell everyone we have dissimilarities
  bool dissimilaritiesWereSet() const { return true; }
/// Get the dissimilarity between two data points 
  double getDissimilarity( const unsigned& , const unsigned& );
/// Get the weight from the input file
  double getWeight( const unsigned& idata ) const ;
/// Just tell plumed to stop
  void update();
/// Read in the dissimilarity matrix
  void performAnalysis();
/// Overwrite virtual function in base class
  void performTask( const unsigned& , const unsigned& , MultiValue& ) const { plumed_error(); }
};

PLUMED_REGISTER_ACTION(ReadDissimilarityMatrix,"READ_DISSIMILARITY_MATRIX")

void ReadDissimilarityMatrix::registerKeywords( Keywords& keys ){
  AnalysisBase::registerKeywords( keys ); keys.remove("USE_OUTPUT_DATA_FROM");
  keys.add("compulsory","FILE","an input file containing the matrix of dissimilarities");
  keys.add("optional","WFILE","input file containing weights of points");
}

ReadDissimilarityMatrix::ReadDissimilarityMatrix( const ActionOptions& ao ):
Action(ao),
AnalysisBase(ao),
nnodes(1)
{
  if( plumed.getActionSet().size()!=0 ) error("read dissimilarity matrix command must be at top of input file");

  parse("FILE",fname);
  log.printf("  reading dissimilarity matrix from file %s \n",fname.c_str() );
  parse("WFILE",wfile);
  if( wfile.length()>0 ) log.printf("  reading weights of nodes from file named %s \n",wfile.c_str() );
  else log.printf("  setting weights of all nodes equal to one\n");

  // We have to set the information on how often we are reading data from the trajectory and how 
  // often we are running here so that the code behaves
  use_all_data=true; freq=1; setStride(1);
}

void ReadDissimilarityMatrix::update(){ plumed.stop(); }

void ReadDissimilarityMatrix::performAnalysis(){
  IFile mfile; mfile.open(fname); 
  // Read in first line
  std::vector<std::string> words; Tools::getParsedLine( mfile, words );
  nnodes=words.size(); dissimilarities.resize( nnodes, nnodes ); 
  for(unsigned j=0;j<nnodes;++j) Tools::convert( words[j], dissimilarities(0,j) );

  for(unsigned i=1;i<nnodes;++i){
       Tools::getParsedLine( mfile, words );
       if( words.size()!=nnodes ) error("bad formatting in matrix file");
       for(unsigned j=0;j<nnodes;++j) Tools::convert( words[j], dissimilarities(i,j) ); 
  }
  mfile.close();

  weights.resize( nnodes );
  if( wfile.length()>0 ){
     IFile wfilef; wfilef.open(wfile); 
     for(unsigned i=0;i<nnodes;++i){ 
       Tools::getParsedLine( wfilef, words ); Tools::convert( words[0], weights[i] );
     }
     wfilef.close();
  } else {
     weights.assign(weights.size(),1.0);
  }
}

double ReadDissimilarityMatrix::getDissimilarity( const unsigned& iframe, const unsigned& jframe ){
  return dissimilarities( iframe, jframe );
}

ReferenceConfiguration* ReadDissimilarityMatrix::getReferenceConfiguration( const unsigned& idata ){
  plumed_merror("cannot get reference configurations from read in dissimilarity matrix");
  return NULL;
}

ReferenceConfiguration* ReadDissimilarityMatrix::getInputReferenceConfiguration( const unsigned& idata ){
  plumed_merror("cannot get reference configurations from read in dissimilarity matrix");
  return NULL;
}

void ReadDissimilarityMatrix::getDataPoint( const unsigned& idata, std::vector<double>& point, double& weight ) const {
  plumed_merror("cannot get data points from read in dissmimilarity matrix");
}

double ReadDissimilarityMatrix::getWeight( const unsigned& idata ) const {
  plumed_assert( idata<nnodes ); return weights[idata]; 
}

}
}
