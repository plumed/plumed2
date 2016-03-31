/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015,2016 The plumed team
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

\par Examples

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
  void calculate();
///
  void performTask( const unsigned& task_index, const unsigned& current, MultiValue& myvals ) const { plumed_error(); }
///
  void turnOnDerivatives();
};

PLUMED_REGISTER_ACTION(ClusterSize,"CLUSTER_NATOMS")

void ClusterSize::registerKeywords( Keywords& keys ){
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
   for(unsigned i=0;i<getNumberOfNodes();++i) addTaskToList(i);
   // And now finish the setup of everything in the base
   std::vector<AtomNumber> fake_atoms; setupMultiColvarBase( fake_atoms );
   addValue(); setNotPeriodic();
}

void ClusterSize::turnOnDerivatives(){
   error("cannot calculate derivatives of number of atoms in cluster.  This quantity is not differentiable");
}

void ClusterSize::calculate(){
   // Retrieve the atoms in the largest cluster
   std::vector<unsigned> myatoms; retrieveAtomsInCluster( clustr, myatoms ); setValue( myatoms.size() );
}

}
}
