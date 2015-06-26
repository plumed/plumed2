/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2015 The plumed team
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
#include "LandmarkSelectionBase.h"
#include "DissimilarityMatrixBase.h"

namespace PLMD {
namespace analysis {

void LandmarkSelectionBase::registerKeywords( Keywords& keys ){
  AnalysisWithAnalysableOutput::registerKeywords( keys );
  keys.add("compulsory","NLANDMARKS","the number of landmarks that you would like to select");
  keys.addFlag("NOVORONOI",false,"do not do a Voronoi analysis of the data to determine weights of final points");
  keys.addFlag("IGNORE_WEIGHTS",false,"ignore the weights in the underlying analysis object");
}

LandmarkSelectionBase::LandmarkSelectionBase( const ActionOptions& ao ):
Action(ao),
AnalysisWithAnalysableOutput(ao)
{
  parse("NLANDMARKS",nlandmarks); 
  log.printf("  selecting %d landmark points \n");
  setNumberOfOutputPoints(nlandmarks);

  parseFlag("NOVORONOI",novoronoi); 
  if( !novoronoi && !dissimilaritiesWereSet() ) error("cannot calculate voronoi weights without dissimilarity mesaure use NOVORONOI or DISSIMILARITIES keyword");

  if( !novoronoi ) log.printf("  ascribing weights to landmarks using voronoi analysis\n");
  else log.printf("  ascribing weights of original points to landmark\n");  
}

void LandmarkSelectionBase::selectFrame( const unsigned& iframe ){
  landmark_indices.push_back( iframe );
}

ReferenceConfiguration* LandmarkSelectionBase::getOutputConfiguration( const unsigned& idata ){
  return getReferenceConfiguration( landmark_indices[idata] );
}

double LandmarkSelectionBase::getOutputDissimilarity( const unsigned& idata, const unsigned& jdata ){
  return getDissimilarity( landmark_indices[idata], landmark_indices[jdata] );
}

void LandmarkSelectionBase::performAnalysis(){
  landmark_indices.resize(0); selectLandmarks(); 
  plumed_dbg_assert( nlandmarks==getNumberOfOutputPoints() );

  if( !novoronoi ){
      std::vector<double> lweights( nlandmarks, 0.0 );
      unsigned rank=comm.Get_rank(), size=comm.Get_size();
      for(unsigned i=rank;i<getNumberOfDataPoints();i+=size){
          unsigned closest=0; 
          double mindist=getDissimilarity( i, landmark_indices[0] );
          for(unsigned j=1;j<nlandmarks;++j){
              double dist=getDissimilarity( i, landmark_indices[j] );
              if( dist<mindist ){ mindist=dist; closest=j; }
          } 
          lweights[closest] += getWeight(i);
      }
      comm.Sum( &lweights[0], lweights.size() ); setOutputWeights( lweights );
  } else {
      std::vector<double> lweights( nlandmarks );
      for(unsigned i=0;i<nlandmarks;++i) lweights[i]=getWeight( landmark_indices[i] );
      setOutputWeights( lweights );
  }
}

}
}
