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
#include "DFSBase.h"
#include "core/ActionRegister.h"

//+PLUMEDOC MCOLVARF DFSCLUSTERDIAMETER
/*
Retrieve the size of a cluster.  This quantity is NOT differentiable.

This action uses the DFS clustering algorithm described in \ref DFSCLUSTERING to find a set of connected components
based on the configuration of the atoms in your system.  Once again this can be used to find crystalline nuclei or 
bubble of atoms.  Once these connected components  you can then find the sizes of the connected components by 
measuring the distance between the two most widely separated atoms in the connected component.  This is what is 
done by this action.

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace adjmat {

class DFSClusterDiameter : public DFSBase {
private:
/// The cluster we are looking for
  unsigned clustr;
public:
/// Create manual
  static void registerKeywords( Keywords& keys );
/// Constructor
  explicit DFSClusterDiameter(const ActionOptions&);
///
  void calculate();
///
  void performTask( const unsigned& task_index, const unsigned& current, MultiValue& myvals ) const ;
///
  void turnOnDerivatives();
};

PLUMED_REGISTER_ACTION(DFSClusterDiameter,"DFSCLUSTERDIAMETER")

void DFSClusterDiameter::registerKeywords( Keywords& keys ){
  DFSBase::registerKeywords( keys );
  keys.add("compulsory","CLUSTER","1","which cluster would you like to look at 1 is the largest cluster, 2 is the second largest, 3 is the the third largest and so on.");
  keys.remove("LOWMEM"); keys.use("HIGHMEM");
}

DFSClusterDiameter::DFSClusterDiameter(const ActionOptions&ao):
Action(ao),
DFSBase(ao)
{
   // Find out which cluster we want
   parse("CLUSTER",clustr);

   if( clustr<1 ) error("cannot look for a cluster larger than the largest cluster");
   if( clustr>getNumberOfNodes() ) error("cluster selected is invalid - too few atoms in system");

   // Create the task list
   for(unsigned  i=1;i<getNumberOfNodes();++i){
       for(unsigned j=0;j<i;++j) addTaskToList( i*getNumberOfNodes() + j );
   }
   // Now create a higest vessel
   addVessel("HIGHEST", "", -1); readVesselKeywords();
}

void DFSClusterDiameter::turnOnDerivatives(){
   error("cannot calculate derivatives of cluster radius.  This quantity is not differentiable");
}

void DFSClusterDiameter::calculate(){
   // Do the clustring
   performClustering();
   // Retrieve the atoms in the largest cluster
   std::vector<unsigned> myatoms; retrieveAtomsInCluster( clustr, myatoms );
   // Activate the relevant tasks
   deactivateAllTasks(); std::vector<unsigned>  active_tasks( getFullNumberOfTasks(), 0 );
   for(unsigned i=1;i<myatoms.size();++i){
       for(unsigned j=0;j<i;++j) active_tasks[ myatoms[i]*getNumberOfNodes() + myatoms[j] ] = 1;  
   }
   activateTheseTasks( active_tasks );
   // Now do the calculation 
   runAllTasks();
}

void DFSClusterDiameter::performTask( const unsigned& task_index, const unsigned& current, MultiValue& myvals ) const { 
  unsigned iatom=current/getNumberOfNodes(), jatom = current - iatom*getNumberOfNodes();
  Vector distance=getSeparation( getPosition(iatom), getPosition(jatom) );
  double dd = distance.modulo(), inv = 1.0/dd ; myvals.setValue( 1, dd ); 
  if( !doNotCalculateDerivatives() ){
      myvals.addDerivative( 1, 3*iatom + 0, -inv*distance[0] );
      myvals.addDerivative( 1, 3*iatom + 1, -inv*distance[1] );
      myvals.addDerivative( 1, 3*iatom + 2, -inv*distance[2] );
      myvals.addDerivative( 1, 3*jatom + 0, +inv*distance[0] );
      myvals.addDerivative( 1, 3*jatom + 1, +inv*distance[1] );
      myvals.addDerivative( 1, 3*jatom + 2, +inv*distance[2] );
      Tensor vir = -inv*Tensor(distance,distance);
      unsigned vbase = myvals.getNumberOfDerivatives() - 9;
      for(unsigned i=0;i<3;++i){
          for(unsigned j=0;j<3;++j) myvals.addDerivative( 1, vbase+3*i+j, vir(i,j) );
      }
  }
}

}
}
