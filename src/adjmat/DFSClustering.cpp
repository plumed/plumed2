/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2019 The plumed team
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
#include "ClusteringBase.h"
#include "AdjacencyMatrixVessel.h"
#include "core/ActionRegister.h"

#ifdef __PLUMED_HAS_BOOST_GRAPH
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/graph_utility.hpp>
#endif

//+PLUMEDOC MATRIXF DFSCLUSTERING
/*
Find the connected components of the matrix using the depth first search clustering algorithm.

As discussed in the section of the manual on \ref contactmatrix a useful tool for developing complex collective variables is the notion of the
so called adjacency matrix.  An adjacency matrix is an \f$N \times N\f$ matrix in which the \f$i\f$th, \f$j\f$th element tells you whether
or not the \f$i\f$th and \f$j\f$th atoms/molecules from a set of \f$N\f$ atoms/molecules are adjacent or not.  As detailed in \cite tribello-clustering
these matrices provide a representation of a graph and can thus can be analyzed using tools from graph theory.  This particular action performs
a depth first search clustering to find the connected components of this graph.  You can read more about depth first search here:

https://en.wikipedia.org/wiki/Depth-first_search

This action is useful if you are looking at a phenomenon such as nucleation where the aim is to detect the sizes of the crystalline nuclei that have formed
in your simulation cell.

\par Examples

The input below calculates the coordination numbers of atoms 1-100 and then computes the an adjacency
matrix whose elements measures whether atoms \f$i\f$ and \f$j\f$ are within 0.55 nm of each other.  The action
labelled dfs then treats the elements of this matrix as zero or ones and thus thinks of the matrix as defining
a graph.  This dfs action then finds the largest connected component in this graph.  The sum of the coordination
numbers for the atoms in this largest connected component are then computed and this quantity is output to a colvar
file.  The way this input can be used is described in detail in \cite tribello-clustering.

\plumedfile
lq: COORDINATIONNUMBER SPECIES=1-100 SWITCH={CUBIC D_0=0.45  D_MAX=0.55} LOWMEM
cm: CONTACT_MATRIX ATOMS=lq  SWITCH={CUBIC D_0=0.45  D_MAX=0.55}
dfs: DFSCLUSTERING MATRIX=cm
clust1: CLUSTER_PROPERTIES CLUSTERS=dfs CLUSTER=1 SUM
PRINT ARG=clust1.* FILE=colvar
\endplumedfile

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace adjmat {

class DFSClustering : public ClusteringBase {
private:
#ifdef __PLUMED_HAS_BOOST_GRAPH
/// The list of edges in the graph
  std::vector<std::pair<unsigned,unsigned> > edge_list;
#else
/// The number of neighbors each atom has
  std::vector<unsigned> nneigh;
/// The adjacency list
  Matrix<unsigned> adj_list;
/// The color that tells us whether a node has been visited
  std::vector<unsigned> color;
/// The recursive function at the heart of this method
  int explore( const unsigned& index );
#endif
public:
/// Create manual
  static void registerKeywords( Keywords& keys );
/// Constructor
  explicit DFSClustering(const ActionOptions&);
/// Do the clustering
  void performClustering() override;
};

PLUMED_REGISTER_ACTION(DFSClustering,"DFSCLUSTERING")

void DFSClustering::registerKeywords( Keywords& keys ) {
  ClusteringBase::registerKeywords( keys );
  keys.add("compulsory","MAXCONNECT","0","maximum number of connections that can be formed by any given node in the graph. "
           "By default this is set equal to zero and the number of connections is set equal to the number "
           "of nodes.  You only really need to set this if you are working with a very large system and "
           "memory is at a premium");
}

DFSClustering::DFSClustering(const ActionOptions&ao):
  Action(ao),
  ClusteringBase(ao)
{
  unsigned maxconnections; parse("MAXCONNECT",maxconnections);
#ifdef __PLUMED_HAS_BOOST_GRAPH
  if( maxconnections>0 ) edge_list.resize( getNumberOfNodes()*maxconnections );
  else edge_list.resize(0.5*getNumberOfNodes()*(getNumberOfNodes()-1));
#else
  nneigh.resize( getNumberOfNodes() ); color.resize(getNumberOfNodes());
  if( maxconnections>0 ) adj_list.resize(getNumberOfNodes(),maxconnections);
  else adj_list.resize(getNumberOfNodes(),getNumberOfNodes());
#endif
}

void DFSClustering::performClustering() {
#ifdef __PLUMED_HAS_BOOST_GRAPH
  // Get the list of edges
  unsigned nedges=0; getAdjacencyVessel()->retrieveEdgeList( nedges, edge_list );

  // Build the graph using boost
  boost::adjacency_list<boost::vecS,boost::vecS,boost::undirectedS> sg(&edge_list[0],&edge_list[nedges],getNumberOfNodes());

  // Find the connected components using boost (-1 here for compatibility with non-boost version)
  number_of_cluster=boost::connected_components(sg,&which_cluster[0]) - 1;

  // And work out the size of each cluster
  for(unsigned i=0; i<which_cluster.size(); ++i) cluster_sizes[which_cluster[i]].first++;
#else
  // Get the adjacency matrix
  getAdjacencyVessel()->retrieveAdjacencyLists( nneigh, adj_list );

  // Perform clustering
  number_of_cluster=-1; color.assign(color.size(),0);
  for(unsigned i=0; i<getNumberOfNodes(); ++i) {
    if( color[i]==0 ) { number_of_cluster++; color[i]=explore(i); }
  }
#endif
}

#ifndef __PLUMED_HAS_BOOST_GRAPH
int DFSClustering::explore( const unsigned& index ) {

  color[index]=1;
  for(unsigned i=0; i<nneigh[index]; ++i) {
    unsigned j=adj_list(index,i);
    if( color[j]==0 ) color[j]=explore(j);
  }

  // Count the size of the cluster
  cluster_sizes[number_of_cluster].first++;
  which_cluster[index] = number_of_cluster;
  return color[index];
}
#endif

}
}
