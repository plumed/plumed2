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
#include "DissimilarityMatrixBase.h"
#include "core/PlumedMain.h"
#include "core/ActionRegister.h"
#include "tools/IFile.h"

namespace PLMD {
namespace analysis {

class ReadDissimilarityMatrix : public DissimilarityMatrixBase {
private:
  unsigned nnodes;
  std::string fname, wfile;
  std::vector<double> weights;
public:
  static void registerKeywords( Keywords& keys );
  ReadDissimilarityMatrix( const ActionOptions& ao );
  unsigned getNumberOfDataPoints() const { return nnodes; }
  ReferenceConfiguration* getReferenceConfiguration( const unsigned& idata );
  void getDataPoint( const unsigned& idata, std::vector<double>& point, double& weight ) const ;
  void calcDissimilarity( const unsigned& , const unsigned& ){ plumed_error(); }
  void update();
  void performAnalysis();
  double getOutputWeight( const unsigned& idata ) const ;
};

PLUMED_REGISTER_ACTION(ReadDissimilarityMatrix,"READ_DISSIMILARITY_MATRIX")

void ReadDissimilarityMatrix::registerKeywords( Keywords& keys ){
  DissimilarityMatrixBase::registerKeywords( keys );
  keys.add("compulsory","FILE","an input file containing the matrix of dissimilarities");
  keys.add("optional","WFILE","input file containing weights of points");
  keys.remove("ATOMS"); keys.remove("METRIC"); keys.remove("STRIDE"); keys.remove("RUN");
  keys.remove("REUSE_INPUT_DATA_FROM"); keys.remove("USE_OUTPUT_DATA_FROM"); keys.remove("USE_ALL_DATA");
  keys.remove("FMT"); keys.remove("REWEIGHT_BIAS"); keys.remove("TEMP"); keys.remove("REWEIGHT_TEMP");
  keys.remove("WRITE_CHECKPOINT");
}

ReadDissimilarityMatrix::ReadDissimilarityMatrix( const ActionOptions& ao ):
Action(ao),
DissimilarityMatrixBase(ao),
nnodes(1)
{
  parse("FILE",fname);
  log.printf("  reading dissimilarity matrix from file %s \n",fname.c_str() );
  parse("WFILE",wfile);
  if( wfile.length()>0 ) log.printf("  reading weights of nodes from file named %s \n",wfile.c_str() );
  else log.printf("  setting weights of all nodes equal to one\n");
}

void ReadDissimilarityMatrix::update(){ plumed.stop(); }

void ReadDissimilarityMatrix::performAnalysis(){
  IFile mfile; mfile.open(fname); 
  // Read in first line
  std::vector<std::string> words; Tools::getParsedLine( mfile, words );
  nnodes=words.size(); mydissimilarities.resize( nnodes, nnodes ); 
  for(unsigned j=0;j<nnodes;++j) Tools::convert( words[j], mydissimilarities(0,j) );

  for(unsigned i=1;i<nnodes;++i){
       Tools::getParsedLine( mfile, words );
       if( words.size()!=nnodes ) error("bad formatting in matrix file");
       for(unsigned j=0;j<nnodes;++j) Tools::convert( words[j], mydissimilarities(i,j) ); 
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

ReferenceConfiguration* ReadDissimilarityMatrix::getReferenceConfiguration( const unsigned& idata ){
  plumed_merror("cannot get reference configurations from read in dissimilarity matrix");
  return NULL;
}

void ReadDissimilarityMatrix::getDataPoint( const unsigned& idata, std::vector<double>& point, double& weight ) const {
  plumed_merror("cannot get data points from read in dissmimilarity matrix");
}

double ReadDissimilarityMatrix::getOutputWeight( const unsigned& idata ) const {
  plumed_dbg_assert( idata<nnodes ); return weights[idata]; 
}


}
}
