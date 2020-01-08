/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2020 The plumed team
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
#include "AdjacencyMatrixVessel.h"
#include "AdjacencyMatrixBase.h"
#include "vesselbase/ActionWithVessel.h"

namespace PLMD {
namespace adjmat {

void AdjacencyMatrixVessel::registerKeywords( Keywords& keys ) {
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
  if( !symmetric && !hbonds ) return (function->ablocks[1].size())*ielem + jelem;
  if( !symmetric ) {
    plumed_dbg_assert( ielem!=jelem );
    if( jelem<ielem ) return (function->ablocks[1].size()-1)*ielem + jelem;
    return (function->ablocks[1].size()-1)*ielem + jelem - 1;
  }
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

void AdjacencyMatrixVessel::retrieveMatrix( DynamicList<unsigned>& myactive_elements, Matrix<double>& mymatrix ) {
  myactive_elements.deactivateAll(); std::vector<double> vals( getNumberOfComponents() );
  for(unsigned i=0; i<getNumberOfStoredValues(); ++i) {
    retrieveSequentialValue( i, false, vals );
    if( vals[0]<epsilon ) continue ;

    myactive_elements.activate(i);
    unsigned j, k; getMatrixIndices( function->getPositionInFullTaskList(i), k, j );

    if( symmetric ) mymatrix(k,j)=mymatrix(j,k)=vals[0]*vals[1];
    else mymatrix(k,j)=vals[0]*vals[1];
  }
  myactive_elements.updateActiveMembers();
}

void AdjacencyMatrixVessel::retrieveAdjacencyLists( std::vector<unsigned>& nneigh, Matrix<unsigned>& adj_list ) {
  plumed_dbg_assert( undirectedGraph() );
  // Currently everything has zero neighbors
  for(unsigned i=0; i<nneigh.size(); ++i) nneigh[i]=0;

  // And set up the adjacency list
  std::vector<double> myvals( getNumberOfComponents() );
  for(unsigned i=0; i<getNumberOfStoredValues(); ++i) {
    // Check if atoms are connected
    retrieveSequentialValue( i, false, myvals );
    if( myvals[0]<epsilon || myvals[1]<epsilon ) continue ;

    unsigned j, k; getMatrixIndices( function->getPositionInFullTaskList(i), k, j );

    if( nneigh[j]>=adj_list.ncols() || nneigh[k]>=adj_list.ncols() ) error("adjacency lists are not large enough, increase maxconnections");
    // Store if atoms are connected
    // unsigned j, k; getMatrixIndices( i, k, j );
    adj_list(k,nneigh[k])=j; nneigh[k]++;
    adj_list(j,nneigh[j])=k; nneigh[j]++;
  }
}

void AdjacencyMatrixVessel::retrieveEdgeList( unsigned& nedge, std::vector<std::pair<unsigned,unsigned> >& edge_list ) {
  plumed_dbg_assert( undirectedGraph() ); nedge=0;
  std::vector<double> myvals( getNumberOfComponents() );
  if( getNumberOfStoredValues()>edge_list.size() ) error("adjacency lists are not large enough, increase maxconnections");

  for(unsigned i=0; i<getNumberOfStoredValues(); ++i) {
    // Check if atoms are connected
    retrieveSequentialValue( i, false, myvals );
    if( myvals[0]<epsilon || myvals[1]<epsilon ) continue ;

    getMatrixIndices( function->getPositionInFullTaskList(i), edge_list[nedge].first, edge_list[nedge].second );
    nedge++;
  }
}

bool AdjacencyMatrixVessel::nodesAreConnected( const unsigned& iatom, const unsigned& jatom ) const {
  if( !matrixElementIsActive( iatom, jatom ) ) return false;
  unsigned ind=getStoreIndexFromMatrixIndices( iatom, jatom );

  std::vector<double> myvals( getNumberOfComponents() );
  retrieveValueWithIndex( ind, false, myvals );
  return ( myvals[0]>epsilon && myvals[1]>epsilon );
}

double AdjacencyMatrixVessel::getCutoffForConnection() const {
  return function->getLinkCellCutoff();
}

}
}

