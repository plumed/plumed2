/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2015 The plumed team
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
#include "vesselbase/ActionWithVessel.h"
#include "AdjacencyMatrixVessel.h"
#include "AdjacencyMatrixBase.h" 

namespace PLMD {
namespace adjmat {

void AdjacencyMatrixVessel::registerKeywords( Keywords& keys ){
  StoreDataVessel::registerKeywords(keys);
}

AdjacencyMatrixVessel::AdjacencyMatrixVessel( const vesselbase::VesselOptions& da ):
StoreDataVessel(da)
{
  function=dynamic_cast<AdjacencyMatrixBase*>( getAction() );
  plumed_assert( function );
}

void AdjacencyMatrixVessel::prepare(){
  finished=false; 
  StoreDataVessel::prepare();
}

void AdjacencyMatrixVessel::setFinishedTrue(){
  finished=true;
}

bool AdjacencyMatrixVessel::calculate( const unsigned& current, MultiValue& myvals, std::vector<double>& buffer, std::vector<unsigned>& der_list ) const {
  if( !finished ) return StoreDataVessel::calculate( current, myvals, buffer, der_list );
  return false;
}

void AdjacencyMatrixVessel::finish( const std::vector<double>& buffer ){
  if( !finished ){
     finished=true;
     StoreDataVessel::finish( buffer );
     function->dertime=true;
  }
}

AdjacencyMatrixBase* AdjacencyMatrixVessel::getMatrixAction() {
  return function;
}

void AdjacencyMatrixVessel::getMatrixIndices( const unsigned& code, unsigned& i, unsigned& j ){
  std::vector<unsigned> myatoms(2); function->decodeIndexToAtoms( function->getTaskCode(code), myatoms ); i=myatoms[0]; j=myatoms[1]; 
}

void AdjacencyMatrixVessel::retrieveMatrix( DynamicList<unsigned>& myactive_elements, Matrix<double>& mymatrix ){
  myactive_elements.deactivateAll();
  std::vector<unsigned> myatoms(2); std::vector<double> vals(2);
  for(unsigned i=0;i<getNumberOfStoredValues();++i){
      // Ignore any non active members
      if( !storedValueIsActive(i) ) continue ;
      myactive_elements.activate(i);
      function->decodeIndexToAtoms( function->getTaskCode(i), myatoms ); 
      unsigned j = myatoms[1], k = myatoms[0];
      retrieveValue( i, false, vals );
      mymatrix(k,j)=mymatrix(j,k)=vals[1];
  }
  myactive_elements.updateActiveMembers();  
}

void AdjacencyMatrixVessel::retrieveAdjacencyLists( std::vector<unsigned>& nneigh, Matrix<unsigned>& adj_list ){
  // Currently everything has zero neighbors
  for(unsigned i=0;i<nneigh.size();++i) nneigh[i]=0;

  // And set up the adjacency list
  std::vector<unsigned> myatoms(2);
  for(unsigned i=0;i<getNumberOfStoredValues();++i){
      // Ignore any non active members
      if( !storedValueIsActive(i) ) continue ;
      function->decodeIndexToAtoms( function->getTaskCode(i), myatoms );
      unsigned j = myatoms[1], k = myatoms[0];
      adj_list(k,nneigh[k])=j; nneigh[k]++;
      adj_list(j,nneigh[j])=k; nneigh[j]++;
  }
}

}
}

