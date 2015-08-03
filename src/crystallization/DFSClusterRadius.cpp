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
namespace crystallization {

class DFSClusterDiameter : public DFSClustering {
private:
/// The cluster we are looking for
  unsigned clustr;
public:
/// Create manual
  static void registerKeywords( Keywords& keys );
/// Constructor
  explicit DFSClusterDiameter(const ActionOptions&);
///
  void doCalculationOnCluster();
///
  void turnOnDerivatives();
};

PLUMED_REGISTER_ACTION(DFSClusterDiameter,"DFSCLUSTERDIAMETER")

void DFSClusterDiameter::registerKeywords( Keywords& keys ){
  DFSClustering::registerKeywords( keys );
  keys.add("compulsory","CLUSTER","1","which cluster would you like to look at 1 is the largest cluster, 2 is the second largest, 3 is the the third largest and so on.");
  keys.use("WTOL"); keys.use("USE_ORIENTATION");
  keys.remove("LOWMEM"); keys.use("HIGHMEM");
}

DFSClusterDiameter::DFSClusterDiameter(const ActionOptions&ao):
Action(ao),
DFSClustering(ao)
{
   // Find out which cluster we want
   parse("CLUSTER",clustr);

   if( clustr<1 ) error("cannot look for a cluster larger than the largest cluster");
   if( clustr>getFullNumberOfBaseTasks() ) error("cluster selected is invalid - too few atoms in system");

   addValue(); setNotPeriodic();
}

void DFSClusterDiameter::turnOnDerivatives(){
   error("cannot calculate derivatives of cluster radius.  This quantity is not differentiable");
}

void DFSClusterDiameter::doCalculationOnCluster(){
   std::vector<unsigned> myatoms; retrieveAtomsInCluster( clustr, myatoms );
   unsigned size=comm.Get_size(), rank=comm.Get_rank();

   double maxdist=0;
   for(unsigned i=1;i<myatoms.size();++i){
       unsigned iatom = myatoms[i];
       Vector ipos=getPositionOfAtomForLinkCells( iatom ); 
       for(unsigned j=0;j<i;++j){
           unsigned jatom = myatoms[j];
           Vector jpos=getPositionOfAtomForLinkCells( jatom );
           Vector distance=getSeparation( ipos, jpos );
           double dmod=distance.modulo2();
           if( dmod>maxdist ) maxdist=dmod;
       }
   }
   setValue( sqrt( maxdist ) );
}

}
}
