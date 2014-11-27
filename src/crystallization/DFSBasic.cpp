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
#include "DFSClustering.h"
#include "core/ActionRegister.h"

//+PLUMEDOC MCOLVARF DFSCLUSTERING 
/*
Find average properites of atoms in a cluster

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace crystallization {

class DFSBasic : public DFSClustering {
private:
/// The cluster we are looking for
  unsigned clustr;
public:
/// Create manual
  static void registerKeywords( Keywords& keys );
/// Constructor
  DFSBasic(const ActionOptions&);
///
  void doCalculationOnCluster();
};

PLUMED_REGISTER_ACTION(DFSBasic,"DFSCLUSTERING")

void DFSBasic::registerKeywords( Keywords& keys ){
  DFSClustering::registerKeywords( keys );
  keys.add("compulsory","CLUSTER","1","which cluster would you like to look at 1 is the largest cluster, 2 is the second largest, 3 is the the third largest and so on.");
  keys.use("WTOL"); keys.use("USE_ORIENTATION");
  keys.use("MEAN"); keys.use("MORE_THAN"); keys.use("LESS_THAN");
  if( keys.reserved("VMEAN") ) keys.use("VMEAN");
  if( keys.reserved("VSUM") ) keys.use("VSUM");
  keys.use("BETWEEN"); keys.use("HISTOGRAM"); keys.use("MOMENTS");
  keys.use("MIN"); keys.use("MAX"); keys.use("SUM"); keys.remove("LOWMEM"); keys.use("HIGHMEM");
}

DFSBasic::DFSBasic(const ActionOptions&ao):
Action(ao),
DFSClustering(ao)
{
   // Find out which cluster we want
   parse("CLUSTER",clustr);

   if( clustr<1 ) error("cannot look for a cluster larger than the largest cluster");
   if( clustr>getFullNumberOfBaseTasks() ) error("cluster selected is invalid - too few atoms in system");

   // Setup the various things this will calculate
   readVesselKeywords();
}

void DFSBasic::doCalculationOnCluster(){
   std::vector<unsigned> myatoms; retrieveAtomsInCluster( clustr, myatoms );
   unsigned size=comm.Get_size(), rank=comm.Get_rank();

   // Now calculate properties of the largest cluster 
   ActionWithVessel::doJobsRequiredBeforeTaskList();  // Note we loose adjacency data by doing this
   // Get size for buffer
   unsigned bsize=0; std::vector<double> buffer( getSizeOfBuffer( bsize ), 0.0 );
   std::vector<double> vals( getNumberOfQuantities() ); std::vector<unsigned> der_index;
   MultiValue myvals( getNumberOfQuantities(), getNumberOfDerivatives() );
   MultiValue bvals( getNumberOfQuantities(), getNumberOfDerivatives() );
   // Get rid of bogus derivatives
   clearDerivatives(); getAdjacencyVessel()->setFinishedTrue();
   for(unsigned j=rank;j<myatoms.size();j+=size){
       // Note loop above over array containing atoms so this is load balanced
       unsigned i=myatoms[j];
       // Need to copy values from base action
       getVectorForTask( i, false, vals );
       for(unsigned i=0;i<vals.size();++i) myvals.setValue( i, vals[i] );
       if( !doNotCalculateDerivatives() ) getVectorDerivatives( i, false, myvals );
       // Run calculate all vessels
       calculateAllVessels( i, myvals, bvals, buffer, der_index );
       myvals.clearAll();
   }
   // MPI Gather everything
   if( buffer.size()>0 ) comm.Sum( buffer );
   finishComputations( buffer );
}

}
}
