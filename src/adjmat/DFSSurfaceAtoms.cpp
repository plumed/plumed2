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
#include "AdjacencyMatrixVessel.h"
#include "core/ActionRegister.h"

//+PLUMEDOC MATRIXF DFSCLUSTERING_WITHSURFACE 
/*
Find the various connected components in an adjacency matrix and then output average
properties of the atoms in those connected components.

\par Examples


*/
//+ENDPLUMEDOC

namespace PLMD {
namespace adjmat {

class DFSBasicWithSurface : public DFSBase {
private:
/// The cluster we are looking for
  unsigned clustr;
/// The cutoff for surface atoms
  double rcut_surf2;
public:
/// Create manual
  static void registerKeywords( Keywords& keys );
/// Constructor
  explicit DFSBasicWithSurface(const ActionOptions&);
/// Do the calculation
  void calculate();
///
  void retrieveAtomsInCluster( const unsigned& clust, std::vector<unsigned>& myatoms ) const ;
/// We can use ActionWithVessel to run all the calculation
  void performTask( const unsigned& , const unsigned& , MultiValue& ) const {}
};

PLUMED_REGISTER_ACTION(DFSBasicWithSurface,"DFSCLUSTERING_WITHSURFACE")

void DFSBasicWithSurface::registerKeywords( Keywords& keys ){
  DFSBase::registerKeywords( keys );
  keys.add("compulsory","CLUSTER","1","which cluster would you like to look at 1 is the largest cluster, 2 is the second largest, 3 is the the third largest and so on.");
  keys.add("compulsory","RCUT_SURF","you also have the option to find the atoms on the surface of the cluster.  An atom must be within this distance of one of the atoms "
                                  "of the cluster in order to be considered a surface atom");
}

DFSBasicWithSurface::DFSBasicWithSurface(const ActionOptions&ao):
Action(ao),
DFSBase(ao)
{
   // Find out which cluster we want
   parse("CLUSTER",clustr);

   if( clustr<1 ) error("cannot look for a cluster larger than the largest cluster");
   if( clustr>getNumberOfNodes() ) error("cluster selected is invalid - too few atoms in system");

   // Setup switching function for surface atoms
   double rcut_surf; parse("RCUT_SURF",rcut_surf);
   if( rcut_surf>0 ) log.printf("  counting surface atoms that are within %f of the cluster atoms \n",rcut_surf);
   rcut_surf2=rcut_surf*rcut_surf;

   // Create the task list
   for(unsigned i=0;i<getNumberOfNodes();++i) addTaskToList(i);

   // Setup the various things this will calculate
   readVesselKeywords();

   // Create the value
   addValue(); setNotPeriodic();
}

void DFSBasicWithSurface::retrieveAtomsInCluster( const unsigned& clust, std::vector<unsigned>& myatoms ) const {
  std::vector<unsigned> tmpat; DFSBase::retrieveAtomsInCluster( clust, tmpat );
  // Find the atoms in the the clusters
  std::vector<bool> atoms( getNumberOfNodes(), false ); 
  for(unsigned i=0;i<tmpat.size();++i){
      for(unsigned j=0;j<getNumberOfNodes();++j){
         double dist2=getSeparation( getPosition(tmpat[i]), getPosition(j)).modulo2();
         if( dist2<rcut_surf2 ){ atoms[j]=true; }
      }
  }
  unsigned nsurf_at=0; 
  for(unsigned j=0;j<getNumberOfNodes();++j){
     if( atoms[j] ) nsurf_at++; 
  }
  myatoms.resize( nsurf_at + tmpat.size() );
  for(unsigned i=0;i<tmpat.size();++i) myatoms[i]=tmpat[i];
  unsigned nn=tmpat.size();
  for(unsigned j=0;j<getNumberOfNodes();++j){
      if( atoms[j] ){ myatoms[nn]=j; nn++; }
  }
  plumed_assert( nn==myatoms.size() );
}

void DFSBasicWithSurface::calculate(){
   // Do the clustring
   performClustering();
   // Retrieve the atoms in the largest cluster
   std::vector<unsigned> myatoms; retrieveAtomsInCluster( clustr, myatoms );
   // The size of the cluster
   setValue( myatoms.size() );
}

}
}
