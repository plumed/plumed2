/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2020 The plumed team
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
#include "core/ActionRegister.h"

#ifdef __PLUMED_HAS_BOOST_GRAPH
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/graph_utility.hpp>
#endif

//+PLUMEDOC MATRIXF DFSCLUSTERING
/*
Find the connected components of the matrix using the depth first search clustering algorithm.

This action is useful if you are looking at a phenomenon such as nucleation where the aim is to detect the sizes of the crystalline nuclei that have formed
in your simulation cell and is used in conjustion with tools from the adjmat module.  The adjmat modlule contains a number of different methods that can be
used to calculate adjacency matrices. An adjacency matrix is an $N \times N$ matrix
in which the $i$th, $j$th element tells you whether or not the $i$th and $j$th atoms/molecules from a set of $N$ atoms/molecules are adjacent or not.  The
simplest of these adjcency matrix methods is the one for calculating a [CONTACT_MATRIX](CONTACT_MATRIX.md).  All adjacency matrices provide a representation of a graph and can thus
can be analyzed using tools from graph theory. This particular action thus performs
[a depth first search clustering](https://en.wikipedia.org/wiki/Depth-first_search) to find the connected components of this graph.

The input below calculates a CONTACT_MATRIX for 100 atoms. In the graph that is created from this CONTACT_MATRIX atoms are connected if they are within
a distance of D_MAX of each other and are disconnected otherwise.  The [DFSCLUSTERING](DFSCLUSTERING.md) method is used to find the connected components in this graph.  This
method outputs a 100-dimensional vector.  The 1st element in this vector tells you which cluster the first atom is within, the second element tells you which
cluster the second atom is in and so on.  If an atom is in the largest cluster its corresponding element in the vector `dfs` will be one. We can thus print
the positions of the atoms in the largest cluster by using a [DUMPATOMS](DUMPATOMS.md) command like the one shown below:

```plumed
cm: CONTACT_MATRIX GROUP=1-100 SWITCH={CUBIC D_0=0.45  D_MAX=0.55}
dfs: DFSCLUSTERING ARG=cm
DUMPATOMS FILE=cluster.xyz ATOMS=1-100 ARG=dfs LESS_THAN_OR_EQUAL=1.5 GREATER_THAN_OR_EQUAL=0.5
```

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace clusters {

class DFSClustering : public ClusteringBase {
private:
#ifndef __PLUMED_HAS_BOOST_GRAPH
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
  keys.addDeprecatedFlag("LOWMEM","");
}

DFSClustering::DFSClustering(const ActionOptions&ao):
  Action(ao),
  ClusteringBase(ao) {
#ifndef __PLUMED_HAS_BOOST_GRAPH
  nneigh.resize( getNumberOfNodes() );
  color.resize(getNumberOfNodes());
#endif
  bool lowmem;
  parseFlag("LOWMEM",lowmem);
  if( lowmem ) {
    warning("LOWMEM flag is deprecated and is no longer required for this action");
  }
}

void DFSClustering::performClustering() {
#ifdef __PLUMED_HAS_BOOST_GRAPH
  // Get the list of edges
  unsigned nedges=0;
  retrieveEdgeList( 0, nedges );

  // Build the graph using boost
  boost::adjacency_list<boost::vecS,boost::vecS,boost::undirectedS> sg(&pairs[0],&pairs[nedges],getNumberOfNodes());

  // Find the connected components using boost (-1 here for compatibility with non-boost version)
  number_of_cluster=boost::connected_components(sg,&which_cluster[0]) - 1;

  // And work out the size of each cluster
  for(unsigned i=0; i<which_cluster.size(); ++i) {
    cluster_sizes[which_cluster[i]].first++;
  }
#else
  // Get the adjacency matrix
  retrieveAdjacencyLists( nneigh, adj_list );

  // Perform clustering
  number_of_cluster=-1;
  color.assign(color.size(),0);
  for(unsigned i=0; i<getNumberOfNodes(); ++i) {
    if( color[i]==0 ) {
      number_of_cluster++;
      color[i]=explore(i);
    }
  }
#endif
}

#ifndef __PLUMED_HAS_BOOST_GRAPH
int DFSClustering::explore( const unsigned& index ) {

  color[index]=1;
  for(unsigned i=0; i<nneigh[index]; ++i) {
    unsigned j=adj_list(index,i);
    if( color[j]==0 ) {
      color[j]=explore(j);
    }
  }

  // Count the size of the cluster
  cluster_sizes[number_of_cluster].first++;
  which_cluster[index] = number_of_cluster;
  return color[index];
}
#endif

}
}
