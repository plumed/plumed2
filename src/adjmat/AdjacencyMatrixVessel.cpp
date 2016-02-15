/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2016 The plumed team
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
#include "vesselbase/ActionWithVessel.h"
#include "AdjacencyMatrixVessel.h"
#include "AdjacencyMatrixBase.h" 

namespace PLMD {
namespace adjmat {

void AdjacencyMatrixVessel::registerKeywords( Keywords& keys ){
  StoreDataVessel::registerKeywords(keys);
  keys.addFlag("SYMMETRIC",false,"is the matrix symmetric");
  keys.addFlag("HBONDS",false,"can we think of the matrix as a undirected graph");
}

AdjacencyMatrixVessel::AdjacencyMatrixVessel( const vesselbase::VesselOptions& da ):
StoreDataVessel(da)
{
  function=dynamic_cast<AdjacencyMatrixBase*>( getAction() );
  plumed_assert( function );
  parseFlag("SYMMETRIC",symmetric); parseFlag("HBONDS",hbonds);
  if( symmetric && hbonds ) error("matrix should be either symmetric or hbonds");
  if( symmetric && function->ablocks[0].size()!=function->ablocks[1].size() ) error("matrix is supposed to be symmetric but nrows!=ncols");
  if( hbonds &&  function->ablocks[0].size()!=function->ablocks[1].size() ) error("matrix is supposed to be hbonds but nrows!=ncols");
}

bool AdjacencyMatrixVessel::isSymmetric() const {
  return symmetric;
}

bool AdjacencyMatrixVessel::undirectedGraph() const {
  return ( symmetric || hbonds );
}

unsigned AdjacencyMatrixVessel::getNumberOfRows() const {
  return function->ablocks[0].size();
}

unsigned AdjacencyMatrixVessel::getNumberOfColumns() const {
  return function->ablocks[1].size();
}

bool AdjacencyMatrixVessel::matrixElementIsActive( const unsigned& ielem, const unsigned& jelem ) const {
  return StoreDataVessel::storedValueIsActive( getStoreIndexFromMatrixIndices( ielem, jelem ) );
}

unsigned AdjacencyMatrixVessel::getStoreIndexFromMatrixIndices( const unsigned& ielem, const unsigned& jelem ) const {
  if( !symmetric ) return (function->ablocks[1].size())*ielem + jelem;
  if( ielem>jelem ) return 0.5*ielem*(ielem-1)+jelem;
  return 0.5*jelem*(jelem-1) + ielem; 
}

AdjacencyMatrixBase* AdjacencyMatrixVessel::getMatrixAction() {
  return function;
}

void AdjacencyMatrixVessel::getMatrixIndices( const unsigned& code, unsigned& i, unsigned& j ) const {
  std::vector<unsigned> myatoms; function->decodeIndexToAtoms( function->getTaskCode(code), myatoms ); 
  i=myatoms[0]; j=myatoms[1];   
  if( !undirectedGraph() ) j -= function->ablocks[0].size(); // Have to remove number of columns as returns number in ablocks[1]  
}

void AdjacencyMatrixVessel::retrieveMatrix( DynamicList<unsigned>& myactive_elements, Matrix<double>& mymatrix ){
  unsigned vin; double df;
  myactive_elements.deactivateAll(); std::vector<double> vals( getNumberOfComponents() ); 
  for(unsigned i=0;i<getNumberOfStoredValues();++i){
      retrieveSequentialValue( i, false, vals );
      if( vals[0]<=wtol ) continue ;

      myactive_elements.activate(i);
      unsigned j, k; getMatrixIndices( function->getPositionInFullTaskList(i), k, j );

      if( symmetric ) mymatrix(k,j)=mymatrix(j,k)=function->transformStoredValues( vals, vin, df );      
      else mymatrix(k,j)=function->transformStoredValues( vals, vin, df );                                 
  }
  myactive_elements.updateActiveMembers();  
}

void AdjacencyMatrixVessel::retrieveAdjacencyLists( std::vector<unsigned>& nneigh, Matrix<unsigned>& adj_list ){
  plumed_dbg_assert( undirectedGraph() );
  // Currently everything has zero neighbors
  for(unsigned i=0;i<nneigh.size();++i) nneigh[i]=0;

  // And set up the adjacency list
  std::vector<double> myvals( getNumberOfComponents() );
  for(unsigned i=0;i<getNumberOfStoredValues();++i){
      // Check if atoms are connected 
      retrieveSequentialValue( i, false, myvals );
      if( myvals[0]<=wtol || !function->checkForConnection( myvals ) ) continue ; 

      unsigned j, k; getMatrixIndices( function->getPositionInFullTaskList(i), k, j ); 

      if( nneigh[j]>=adj_list.ncols() || nneigh[k]>=adj_list.ncols() ) error("adjacency lists are not large enough, increase maxconnections"); 
      // Store if atoms are connected
      // unsigned j, k; getMatrixIndices( i, k, j );
      adj_list(k,nneigh[k])=j; nneigh[k]++;
      adj_list(j,nneigh[j])=k; nneigh[j]++;
  }
}

void AdjacencyMatrixVessel::retrieveEdgeList( unsigned& nedge, std::vector<std::pair<unsigned,unsigned> >& edge_list ){
  plumed_dbg_assert( undirectedGraph() ); nedge=0;
  std::vector<double> myvals( getNumberOfComponents() );
  if( getNumberOfStoredValues()>edge_list.size() ) error("adjacency lists are not large enough, increase maxconnections");

  for(unsigned i=0;i<getNumberOfStoredValues();++i){
      // Check if atoms are connected 
      retrieveSequentialValue( i, false, myvals );
      if( myvals[0]<=wtol || !function->checkForConnection( myvals ) ) continue ;

      getMatrixIndices( function->getPositionInFullTaskList(i), edge_list[nedge].first, edge_list[nedge].second );
      nedge++;
  }
}

void AdjacencyMatrixVessel::retrieveDerivatives( const unsigned& myelem, const bool& normed, MultiValue& myvals ){
  StoreDataVessel::retrieveDerivatives( myelem, normed, myvals );
  if( !function->weightHasDerivatives ) return ;

  unsigned vi; std::vector<double> vals( getNumberOfComponents() ); retrieveValueWithIndex( myelem, normed, vals ); 
  double df, max=function->transformStoredValues( vals, vi, df );

  double pref = max/(vals[0]*vals[0]);
  for(unsigned i=0;i<myvals.getNumberActive();++i){
      unsigned jder=myvals.getActiveIndex(i);
      myvals.setDerivative( 1, jder, df*myvals.getDerivative(vi, jder)/vals[0] - pref*myvals.getDerivative(0, jder) );
  }
}

void AdjacencyMatrixVessel::recalculateStoredQuantity( const unsigned& myelem, MultiValue& myvals ){
  function->recalculateMatrixElement( myelem, myvals );
}

}
}

