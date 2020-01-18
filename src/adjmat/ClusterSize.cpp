/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2020 The plumed team
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
#include "ClusterAnalysisBase.h"
#include "core/ActionRegister.h"

//+PLUMEDOC CONCOMP CLUSTER_NATOMS
/*
Gives the number of atoms in the connected component

As discussed in the section of the manual on \ref contactmatrix a useful tool for developing complex collective variables is the notion of the
so called adjacency matrix.  An adjacency matrix is an \f$N \times N\f$ matrix in which the \f$i\f$th, \f$j\f$th element tells you whether
or not the \f$i\f$th and \f$j\f$th atoms/molecules from a set of \f$N\f$ atoms/molecules are adjacent or not.  When analyzing these matrix
we can treat them as a graph and find connected components using some clustering algorithm.  This action is used in tandem with this form of analysis
to output the number of atoms that are connected together in a particular connected component.  It is important to note that the quantity that is
output by this action cannot be differentiated.  As such it cannot be used as a collective variable in a biased simulation.

\par Examples

The following input uses PLUMED to calculate a adjacency matrix that connects a pair of atoms if they both have a coordination number that is greater
than 2.0 and if they are within 6.0 nm of each other.  Depth first search clustering is used to find the connected components in this matrix and then
the number of atoms in the largest cluster is found.  This quantity is then output to a file called colvar

\plumedfile
# Calculate coordination numbers
c1: COORDINATIONNUMBER SPECIES=1-512 SWITCH={EXP D_0=4.0 R_0=0.5 D_MAX=6.0}
# Select coordination numbers that are more than 2.0
cf: MFILTER_MORE DATA=c1 SWITCH={RATIONAL D_0=2.0 R_0=0.1} LOWMEM
# Build a contact matrix
mat: CONTACT_MATRIX ATOMS=cf SWITCH={EXP D_0=4.0 R_0=0.5 D_MAX=6.0}
# Find largest cluster
dfs: DFSCLUSTERING MATRIX=mat
clust1: CLUSTER_PROPERTIES CLUSTERS=dfs CLUSTER=1
nat: CLUSTER_NATOMS CLUSTERS=dfs CLUSTER=1
PRINT ARG=nat FILE=COLVAR
\endplumedfile

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace adjmat {

class ClusterSize : public ClusterAnalysisBase {
private:
/// The cluster we are looking for
  unsigned clustr;
public:
/// Create manual
  static void registerKeywords( Keywords& keys );
/// Constructor
  explicit ClusterSize(const ActionOptions&);
///
  void calculate() override;
///
  void performTask( const unsigned& task_index, const unsigned& current, MultiValue& myvals ) const override { plumed_error(); }
///
  void turnOnDerivatives() override;
};

PLUMED_REGISTER_ACTION(ClusterSize,"CLUSTER_NATOMS")

void ClusterSize::registerKeywords( Keywords& keys ) {
  ClusterAnalysisBase::registerKeywords( keys );
  keys.add("compulsory","CLUSTER","1","which cluster would you like to look at 1 is the largest cluster, 2 is the second largest, 3 is the the third largest and so on.");
}

ClusterSize::ClusterSize(const ActionOptions&ao):
  Action(ao),
  ClusterAnalysisBase(ao)
{
  // Find out which cluster we want
  parse("CLUSTER",clustr);

  if( clustr<1 ) error("cannot look for a cluster larger than the largest cluster");
  if( clustr>getNumberOfNodes() ) error("cluster selected is invalid - too few atoms in system");

  // Create all tasks by copying those from underlying DFS object (which is actually MultiColvar)
  for(unsigned i=0; i<getNumberOfNodes(); ++i) addTaskToList(i);
  // And now finish the setup of everything in the base
  std::vector<AtomNumber> fake_atoms; setupMultiColvarBase( fake_atoms );
  addValue(); setNotPeriodic();
}

void ClusterSize::turnOnDerivatives() {
  error("cannot calculate derivatives of number of atoms in cluster.  This quantity is not differentiable");
}

void ClusterSize::calculate() {
  // Retrieve the atoms in the largest cluster
  std::vector<unsigned> myatoms; retrieveAtomsInCluster( clustr, myatoms ); setValue( myatoms.size() );
}

}
}
