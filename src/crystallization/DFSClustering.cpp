/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014 The plumed team
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
#include "multicolvar/AdjacencyMatrixAction.h"
#include "core/ActionRegister.h"

//+PLUMEDOC MCOLVARF DFSCLUSTERING 
/*
Find average properites of atoms in a cluster

\par Examples

*/
//+ENDPLUMEDOC


namespace PLMD {
namespace crystallization {

class DFSClustering : public multicolvar::AdjacencyMatrixAction {
private:
/// The cluster we are looking for
  unsigned clustr;
/// Used to identify the cluster we are working on 
  int number_of_cluster;
/// The values from the underlying colvar
  std::vector<double> myvals;
/// The number of neighbors each atom has
  std::vector<unsigned> nneigh;
/// The adjacency list
  Matrix<unsigned> adj_list;
/// Vector that stores the sizes of the current set of clusters
  std::vector< std::pair<unsigned,unsigned> > cluster_sizes;
/// Vector that identifies the cluster each atom belongs to
  std::vector<unsigned> which_cluster;
/// The color that tells us whether a node has been visited
  std::vector<unsigned> color;
/// The recursive function at the heart of this method
  int explore( const unsigned& index );
public:
/// Create manual
  static void registerKeywords( Keywords& keys );
/// Constructor
  DFSClustering(const ActionOptions&);
/// Required as we have to be able to deal with vectors
  unsigned getNumberOfQuantities();
/// This checks whether derivatives can be computed given the base multicolvar
  void turnOnDerivatives();
/// Do the matrix calculation
  void completeCalculation();
/// Derivatives of elements of adjacency matrix are unimportant.  We thus
/// overwrite this routine as this makes the code faster
  void updateActiveAtoms(){}
};

PLUMED_REGISTER_ACTION(DFSClustering,"DFSCLUSTERING")

void DFSClustering::registerKeywords( Keywords& keys ){
  multicolvar::AdjacencyMatrixAction::registerKeywords( keys );
  keys.add("compulsory","CLUSTER","1","which cluster would you like to look at 1 is the largest cluster, 2 is the second largest, 3 is the the third largest and so on.");
  keys.use("WTOL"); keys.use("USE_ORIENTATION");
  keys.use("MEAN"); keys.use("MORE_THAN"); keys.use("LESS_THAN"); 
  if( keys.reserved("VMEAN") ) keys.use("VMEAN");
  if( keys.reserved("VSUM") ) keys.use("VSUM");
  keys.use("BETWEEN"); keys.use("HISTOGRAM"); keys.use("MOMENTS"); 
  keys.use("MIN"); keys.use("MAX"); keys.use("SUM"); keys.remove("LOWMEM"); keys.use("HIGHMEM");
}

DFSClustering::DFSClustering(const ActionOptions&ao):
Action(ao),
AdjacencyMatrixAction(ao),
nneigh(getFullNumberOfBaseTasks()),
adj_list(getFullNumberOfBaseTasks(),getFullNumberOfBaseTasks()),
cluster_sizes(getFullNumberOfBaseTasks()),
which_cluster(getFullNumberOfBaseTasks()),
color(getFullNumberOfBaseTasks())
{
   if( getNumberOfBaseMultiColvars()!=1 ) error("should only be running DFS Clustering with one base multicolvar");
   // Find out which cluster we want
   parse("CLUSTER",clustr); 

   if( clustr<1 ) error("cannot look for a cluster larger than the largest cluster");
   if( clustr>getFullNumberOfBaseTasks() ) error("cluster selected is invalid - too few atoms in system");

   // Set myvals size equal to number of components in vector + 1 (for norm)
   myvals.resize( getBaseMultiColvar(0)->getNumberOfQuantities() - 4 );
   // Setup the various things this will calculate
   readVesselKeywords(); 
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
   // Work out which atoms are in the largest cluster
   unsigned n=0; std::vector<unsigned> myatoms( cluster_sizes[cluster_sizes.size() - clustr].first );
   for(unsigned i=0;i<getFullNumberOfBaseTasks();++i){
      if( which_cluster[i]==cluster_sizes[cluster_sizes.size() - clustr].second ){ myatoms[n]=i; n++; }
   }
   unsigned size=comm.Get_size(), rank=comm.Get_rank();

//   for(unsigned i=0;i<cluster_sizes.size();++i) printf("HELLO CLUSTER %d %d \n",i,cluster_sizes[i].first );
   // Now calculate properties of the largest cluster 
   ActionWithVessel::doJobsRequiredBeforeTaskList();  // Note we loose adjacency data by doing this
   // Get rid of bogus derivatives
   clearDerivatives(); getAdjacencyVessel()->setFinishedTrue(); 
   for(unsigned j=rank;j<myatoms.size();j+=size){
       // Note loop above over array containing atoms so this is load balanced
       unsigned i=myatoms[j];
       // Need to copy values from base action
       extractValueForBaseTask( i, myvals );
       //  getValueForBaseTask( i, myvals );
       setElementValue(0, myvals[0] ); setElementValue(1, 1.0 );
       for(unsigned i=0;i<myvals.size()-1;++i) setElementValue(2+i, myvals[1+i] );
       // Prepare dynamic lists
       atoms_with_derivatives.deactivateAll();
       // Copy derivatives from base action
       extractWeightedAverageAndDerivatives( i, 1.0 ); 
       // Update all dynamic lists
       atoms_with_derivatives.updateActiveMembers();
       // Run calculate all vessels
       calculateAllVessels();
       // Must clear element values and derivatives
       clearAfterTask();
   }
   finishComputations();
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

}
}
