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

//+PLUMEDOC MCOLVARF DFSMAXCLUSTER
/*
Find average properites of atoms in a cluster

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace crystallization {

class DFSMaxCluster : public DFSClustering {
private:
/// The value of beta
  double beta;
public:
/// Create manual
  static void registerKeywords( Keywords& keys );
/// Constructor
  DFSMaxCluster(const ActionOptions&);
///
  void doCalculationOnCluster();
};

PLUMED_REGISTER_ACTION(DFSMaxCluster,"DFSMAXCLUSTER")

void DFSMaxCluster::registerKeywords( Keywords& keys ){
  DFSClustering::registerKeywords( keys );
  keys.add("compulsory","BETA","the value of beta to be used in calculating the smooth maximum");
  keys.use("WTOL"); keys.use("USE_ORIENTATION");
  keys.remove("LOWMEM"); keys.use("HIGHMEM");
}

DFSMaxCluster::DFSMaxCluster(const ActionOptions&ao):
Action(ao),
DFSClustering(ao)
{
   // Find out the value of beta
   parse("BETA",beta);
   addValueWithDerivatives(); setNotPeriodic();
}

void DFSMaxCluster::doCalculationOnCluster(){
   unsigned size=comm.Get_size(), rank=comm.Get_rank(); 

   std::vector<double> tder( getNumberOfDerivatives() ), fder( 1+getNumberOfDerivatives(), 0.0 );
   MultiValue myvals( getNumberOfQuantities(), getNumberOfDerivatives() );
   std::vector<unsigned> myatoms; std::vector<double> vals( getNumberOfQuantities() );

   fder[0]=0;  // Value and derivatives are accumulated in one array so that there is only one MPI call 
   for(unsigned iclust=0;iclust<getNumberOfClusters();++iclust){
       retrieveAtomsInCluster( iclust+1, myatoms );

       double tval=0; tder.assign( tder.size(), 0.0 );
       for(unsigned j=rank;j<myatoms.size();j+=size){ 
           unsigned i=myatoms[j];
           getVectorForTask( i, false, vals );
           tval += vals[0]*vals[1];
           if( !doNotCalculateDerivatives() ){ 
               getVectorDerivatives( i, false, myvals );
               for(unsigned k=0;k<myvals.getNumberActive();++k){
                   unsigned kat=myvals.getActiveIndex(k);
                   tder[kat]+=vals[0]*myvals.getDerivative(1,kat) + vals[1]*myvals.getDerivative(0,kat);
               }
               myvals.clearAll();
           }
       }

       // Accumulate value and derivatives
       fder[0] += exp( tval/beta ); double pref = exp(tval/beta) / beta;
       if( !doNotCalculateDerivatives() ){
           for(unsigned k=0;k<getNumberOfDerivatives();++k) fder[1+k]+=pref*tder[k];
       }
   }
   comm.Sum( fder );

   // And finish the derivatives
   setValue( beta*std::log( fder[0] ) ); Value* myval=getPntrToValue();
   if( !doNotCalculateDerivatives() ){
      double pref = beta/fder[0];
      for(unsigned k=0;k<getNumberOfDerivatives();++k) myval->addDerivative( k, pref*fder[1+k] );
   }
}

}
}
