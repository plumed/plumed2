/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2019 The plumed team
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
#include "LandmarkSelectionBase.h"

namespace PLMD {
namespace analysis {

void LandmarkSelectionBase::registerKeywords( Keywords& keys ) {
  AnalysisBase::registerKeywords( keys );
  keys.add("compulsory","NLANDMARKS","the number of landmarks that you would like to select");
  keys.addFlag("NOVORONOI",false,"do not do a Voronoi analysis of the data to determine weights of final points");
  keys.addFlag("IGNORE_WEIGHTS",false,"ignore the weights in the underlying analysis object");
}

LandmarkSelectionBase::LandmarkSelectionBase( const ActionOptions& ao ):
  Action(ao),
  AnalysisBase(ao),
  nlandmarks(0)
{
  if( keywords.exists("NLANDMARKS") ) parse("NLANDMARKS",nlandmarks);
  log.printf("  selecting %u landmark points \n",nlandmarks);
  lweights.resize( nlandmarks );

  parseFlag("NOVORONOI",novoronoi);
  if( !novoronoi && !dissimilaritiesWereSet() ) error("cannot calculate voronoi weights without dissimilarity mesaure");

  if( !novoronoi ) log.printf("  ascribing weights to landmarks using voronoi analysis\n");
  else log.printf("  ascribing weights of original points to landmark\n");
}

void LandmarkSelectionBase::selectFrame( const unsigned& iframe ) {
  landmark_indices.push_back( iframe );
}

void LandmarkSelectionBase::performAnalysis() {
  landmark_indices.resize(0); selectLandmarks();
  plumed_dbg_assert( nlandmarks==getNumberOfDataPoints() );
  if( lweights.size()!=nlandmarks ) lweights.resize( nlandmarks );

  if( !novoronoi ) {
    lweights.assign(lweights.size(),0.0);
    std::vector<unsigned> tmpass( my_input_data->getNumberOfDataPoints() );
    voronoiAnalysis( landmark_indices, lweights, tmpass );
  } else {
    for(unsigned i=0; i<nlandmarks; ++i) lweights[i]=my_input_data->getWeight( landmark_indices[i] );
  }
}

void LandmarkSelectionBase::voronoiAnalysis( const std::vector<unsigned>& myindices, std::vector<double>& lweights, std::vector<unsigned>& assignments ) const {
  plumed_dbg_assert( myindices.size()==lweights.size() && assignments.size()==my_input_data->getNumberOfDataPoints() );
  lweights.assign( lweights.size(), 0 );
  unsigned rank=comm.Get_rank(), size=comm.Get_size();
  for(unsigned i=rank; i<my_input_data->getNumberOfDataPoints(); i+=size) {
    assignments[i]=0;
    double mindist=my_input_data->getDissimilarity( i, myindices[0] );
    for(unsigned j=1; j<nlandmarks; ++j) {
      double dist=my_input_data->getDissimilarity( i, myindices[j] );
      if( dist<mindist ) { mindist=dist; assignments[i]=j; }
    }
    lweights[ assignments[i] ] += my_input_data->getWeight(i);
  }
  comm.Sum( &lweights[0], lweights.size() );
  comm.Sum( &assignments[0], assignments.size() );
}

}
}
