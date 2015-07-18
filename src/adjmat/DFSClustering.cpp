/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014,2015 The plumed team
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
#include "DFSClustering.h"
#include "AdjacencyMatrixBase.h"
#include "AdjacencyMatrixVessel.h"

namespace PLMD {
namespace adjmat {

void DFSClustering::registerKeywords( Keywords& keys ){
  ActionWithInputMatrix::registerKeywords( keys );
}

DFSClustering::DFSClustering(const ActionOptions&ao):
Action(ao),
ActionWithInputMatrix(ao),
number_of_cluster(-1),
nneigh(getNumberOfNodes()),
adj_list(getNumberOfNodes(),getNumberOfNodes()),
cluster_sizes(getNumberOfNodes()),
which_cluster(getNumberOfNodes()),
color(getNumberOfNodes())
{
   if( getAdjacencyVessel()->getNumberOfBaseColvars()!=1 ) error("should only be running DFS Clustering with one base multicolvar in function");
}

void DFSClustering::turnOnDerivatives(){
   // Check base multicolvar isn't density probably other things shouldn't be allowed here as well
   if( (getAdjacencyVessel()->getBaseMultiColvar(0))->isDensity() ) error("DFS clustering cannot be differentiated if base multicolvar is DENSITY");

   // Check for dubious vessels
   for(unsigned i=0;i<getNumberOfVessels();++i){
      if( getPntrToVessel(i)->getName()=="MEAN" ) error("MEAN of cluster is not differentiable");
      if( getPntrToVessel(i)->getName()=="VMEAN" ) error("VMEAN of cluster is not differentiable");  
   }
   
   // Ensure that derivatives are turned on in base classes
   ActionWithValue::turnOnDerivatives();
   needsDerivatives();
}

unsigned DFSClustering::getNumberOfQuantities(){
  return (getAdjacencyVessel()->getBaseMultiColvar(0))->getNumberOfQuantities();
} 

void DFSClustering::calculate(){
   // Get the adjacency matrix
   getAdjacencyVessel()->retrieveAdjacencyLists( nneigh, adj_list ); 

   // All the clusters have zero size initially
   for(unsigned i=0;i<cluster_sizes.size();++i){ cluster_sizes[i].first=0; cluster_sizes[i].second=i;}

   // Perform clustering
   number_of_cluster=-1; color.assign(color.size(),0);
   for(unsigned i=0;i<getNumberOfNodes();++i){
      if( color[i]==0 ){ number_of_cluster++; color[i]=explore(i); } 
   }

   // Order the clusters in the system by size (this returns ascending order )
   std::sort( cluster_sizes.begin(), cluster_sizes.end() );
   // Finish the calculation (what we do depends on what we are inheriting into)
   doCalculationOnCluster();

}

int DFSClustering::explore( const unsigned& index ){

   color[index]=1;
   for(unsigned i=0;i<nneigh[index];++i){
       unsigned j=adj_list(index,i);
       if( color[j]==0 ) color[j]=explore(j);
   }

   // Count the size of the cluster
   cluster_sizes[number_of_cluster].first++;
   which_cluster[index] = number_of_cluster;
   return color[index];
}

void DFSClustering::retrieveAtomsInCluster( const unsigned& clust, std::vector<unsigned>& myatoms ) const {
   unsigned n=0; myatoms.resize( cluster_sizes[cluster_sizes.size() - clust].first );
   for(unsigned i=0;i<getNumberOfNodes();++i){
      if( which_cluster[i]==cluster_sizes[cluster_sizes.size() - clust].second ){ myatoms[n]=i; n++; }
   }
}

bool DFSClustering::isCurrentlyActive( const unsigned& ind ) const {
   return (getAdjacencyVessel()->getBaseData(0))->storedValueIsActive( ind );
}

void DFSClustering::getVectorForTask( const unsigned& ind, const bool& normed, std::vector<double>& orient0 ) const {
   plumed_dbg_assert( (getAdjacencyVessel()->getBaseData(0))->storedValueIsActive( ind ) ); 
   (getAdjacencyVessel()->getBaseData(0))->retrieveValue( ind, normed, orient0 );
}

void DFSClustering::getVectorDerivatives( const unsigned& ind, const bool& normed, MultiValue& myder0 ) const {
   plumed_dbg_assert( (getAdjacencyVessel()->getBaseData(0))->storedValueIsActive( ind ) ); 
   (getAdjacencyVessel()->getBaseData(0))->retrieveDerivatives(  ind, normed, myder0 );
}

}
}
