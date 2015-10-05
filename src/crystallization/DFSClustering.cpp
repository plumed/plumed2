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


namespace PLMD {
namespace crystallization {

void DFSClustering::registerKeywords( Keywords& keys ){
  multicolvar::AdjacencyMatrixAction::registerKeywords( keys );
}

DFSClustering::DFSClustering(const ActionOptions&ao):
Action(ao),
AdjacencyMatrixAction(ao),
number_of_cluster(-1),
nneigh(getFullNumberOfBaseTasks()),
adj_list(getFullNumberOfBaseTasks(),getFullNumberOfBaseTasks()),
cluster_sizes(getFullNumberOfBaseTasks()),
which_cluster(getFullNumberOfBaseTasks()),
color(getFullNumberOfBaseTasks())
{
   if( getNumberOfBaseMultiColvars()!=1 ) error("should only be running DFS Clustering with one base multicolvar");
}

void DFSClustering::turnOnDerivatives(){
   // Check base multicolvar isn't density probably other things shouldn't be allowed here as well
   if( getBaseMultiColvar(0)->isDensity() ) error("DFS clustering cannot be differentiated if base multicolvar is DENSITY");

   // Check for dubious vessels
   for(unsigned i=0;i<getNumberOfVessels();++i){
      if( getPntrToVessel(i)->getName()=="MEAN" ) error("MEAN of cluster is not differentiable");
      if( getPntrToVessel(i)->getName()=="VMEAN" ) error("VMEAN of cluster is not differentiable");  
   }

   MultiColvarBase::turnOnDerivatives();
}

unsigned DFSClustering::getNumberOfQuantities(){
  return getBaseMultiColvar(0)->getNumberOfQuantities();
} 

void DFSClustering::completeCalculation(){
   // Get the adjacency matrix
   retrieveAdjacencyLists( nneigh, adj_list ); 
//   for(unsigned i=0;i<nneigh.size();++i){
//       printf("ADJACENCY LIST FOR %d HAS %d MEMBERS : ",i,nneigh[i]);
//       for(unsigned j=0;j<nneigh[i];++j) printf("%d ",adj_list(i,j) );
//       printf("\n"); 
//   }

   // All the clusters have zero size initially
   for(unsigned i=0;i<cluster_sizes.size();++i){ cluster_sizes[i].first=0; cluster_sizes[i].second=i;}

   // Perform clustering
   number_of_cluster=-1; color.assign(color.size(),0);
   for(unsigned i=0;i<getFullNumberOfBaseTasks();++i){
      if( color[i]==0 ){ number_of_cluster++; color[i]=explore(i); } 
   }

   // Order the clusters in the system by size (this returns ascending order )
   std::sort( cluster_sizes.begin(), cluster_sizes.end() );

   // Finish the calculation (what we do depends on what we are inheriting into)
   doCalculationOnCluster();
   

//   // Work out which atoms are in the largest cluster
//   unsigned n=0; std::vector<unsigned> myatoms( cluster_sizes[cluster_sizes.size() - clustr].first );
//   for(unsigned i=0;i<getFullNumberOfBaseTasks();++i){
//      if( which_cluster[i]==cluster_sizes[cluster_sizes.size() - clustr].second ){ myatoms[n]=i; n++; }
//   }
 //  unsigned size=comm.Get_size(), rank=comm.Get_rank();

//    // Now calculate properties of the largest cluster 
//    ActionWithVessel::doJobsRequiredBeforeTaskList();  // Note we loose adjacency data by doing this
//    // Get size for buffer
//    unsigned bsize=0; std::vector<double> buffer( getSizeOfBuffer( bsize ), 0.0 );
//    std::vector<double> vals( getNumberOfQuantities() ); std::vector<unsigned> der_index;
//    MultiValue myvals( getNumberOfQuantities(), getNumberOfDerivatives() );
//    MultiValue bvals( getNumberOfQuantities(), getNumberOfDerivatives() );
//    // Get rid of bogus derivatives
//    clearDerivatives(); getAdjacencyVessel()->setFinishedTrue(); 
//    for(unsigned j=rank;j<myatoms.size();j+=size){
//        // Note loop above over array containing atoms so this is load balanced
//        unsigned i=myatoms[j];
//        // Need to copy values from base action
//        getVectorForTask( i, false, vals );
//        for(unsigned i=0;i<vals.size();++i) myvals.setValue( i, vals[i] );
//        if( !doNotCalculateDerivatives() ) getVectorDerivatives( i, false, myvals );
//        // Run calculate all vessels
//        calculateAllVessels( i, myvals, bvals, buffer, der_index );
//        myvals.clearAll();
//    }
//    // MPI Gather everything
//    if( buffer.size()>0 ) comm.Sum( buffer );
//    finishComputations( buffer );
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
   for(unsigned i=0;i<getFullNumberOfBaseTasks();++i){
      if( which_cluster[i]==cluster_sizes[cluster_sizes.size() - clust].second ){ myatoms[n]=i; n++; }
   }
}

}
}
