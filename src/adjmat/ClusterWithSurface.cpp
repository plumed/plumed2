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
#include "ClusteringBase.h"
#include "AdjacencyMatrixVessel.h"
#include "AdjacencyMatrixBase.h"
#include "core/ActionRegister.h"

//+PLUMEDOC MATRIXF CLUSTER_WITHSURFACE 
/*
Find the various connected components in an adjacency matrix and then output average
properties of the atoms in those connected components.

\par Examples


*/
//+ENDPLUMEDOC

namespace PLMD {
namespace adjmat {

class ClusterWithSurface : public ClusteringBase {
private:
/// The clusters that we are adding surface atoms to
  ClusteringBase* myclusters;
/// The cutoff for surface atoms
  double rcut_surf2;
public:
/// Create manual
  static void registerKeywords( Keywords& keys );
/// Constructor
  explicit ClusterWithSurface(const ActionOptions&);
///
  unsigned getNumberOfDerivatives();
///
  unsigned getNumberOfNodes() const ;
///
  AtomNumber getAbsoluteIndexOfCentralAtom(const unsigned& i) const ;
///
  void retrieveAtomsInCluster( const unsigned& clust, std::vector<unsigned>& myatoms ) const ;
///
  void getInputData( const unsigned& ind, const bool& normed, const multicolvar::AtomValuePack& myatoms, std::vector<double>& orient0 ) const ;
///
  MultiValue& getInputDerivatives( const unsigned& ind, const bool& normed, const multicolvar::AtomValuePack& myatoms ) const ;
///
  unsigned getNumberOfQuantities() const ;
/// Do the calculation
  void performClustering(){};
///
  double  getCutoffForConnection() const ;  
///
  Vector getPositionOfAtomForLinkCells( const unsigned& taskIndex ) const ;
};

PLUMED_REGISTER_ACTION(ClusterWithSurface,"CLUSTER_WITHSURFACE")

void ClusterWithSurface::registerKeywords( Keywords& keys ){
  ClusteringBase::registerKeywords( keys );
  keys.remove("MATRIX");
  keys.add("compulsory","CLUSTERS","the label of the action that does the clustering");
  keys.add("compulsory","RCUT_SURF","you also have the option to find the atoms on the surface of the cluster.  An atom must be within this distance of one of the atoms "
                                  "of the cluster in order to be considered a surface atom");
}

ClusterWithSurface::ClusterWithSurface(const ActionOptions&ao):
Action(ao),
ClusteringBase(ao)
{
   std::vector<AtomNumber> fake_atoms;
   if( !parseMultiColvarAtomList("CLUSTERS",-1,fake_atoms ) ) error("unable to find CLUSTERS input");
   if( mybasemulticolvars.size()!=1 ) error("should be exactly one multicolvar input");

   // Retrieve the adjacency matrix of interest
   atom_lab.resize(0); myclusters = dynamic_cast<ClusteringBase*>( mybasemulticolvars[0] ); 
   if( !myclusters ) error( mybasemulticolvars[0]->getLabel() + " does not calculate clusters");

   // Setup switching function for surface atoms
   double rcut_surf; parse("RCUT_SURF",rcut_surf);
   if( rcut_surf>0 ) log.printf("  counting surface atoms that are within %f of the cluster atoms \n",rcut_surf);
   rcut_surf2=rcut_surf*rcut_surf;

   // And now finish the setup of everything in the base
   setupMultiColvarBase( fake_atoms ); 
}

unsigned ClusterWithSurface::getNumberOfDerivatives(){
  return myclusters->getNumberOfDerivatives();
}

unsigned ClusterWithSurface::getNumberOfNodes() const {
  return myclusters->getNumberOfNodes();
}

AtomNumber ClusterWithSurface::getAbsoluteIndexOfCentralAtom(const unsigned& i) const {
  return myclusters->getAbsoluteIndexOfCentralAtom(i);
}

void ClusterWithSurface::getInputData( const unsigned& ind, const bool& normed, const multicolvar::AtomValuePack& myatoms, std::vector<double>& orient0 ) const {
  myclusters->getInputData( ind, normed, myatoms, orient0 );
}

MultiValue& ClusterWithSurface::getInputDerivatives( const unsigned& ind, const bool& normed, const multicolvar::AtomValuePack& myatoms ) const {
  return myclusters->getInputDerivatives( ind, normed, myatoms );
}

unsigned ClusterWithSurface::getNumberOfQuantities() const {
  return myclusters->getNumberOfQuantities();
}

double  ClusterWithSurface::getCutoffForConnection() const {
  double tcut = myclusters->getCutoffForConnection();
  if( tcut>sqrt(rcut_surf2) ) return tcut;
  return sqrt(rcut_surf2);
}

void ClusterWithSurface::retrieveAtomsInCluster( const unsigned& clust, std::vector<unsigned>& myatoms ) const {
  std::vector<unsigned> tmpat; myclusters->retrieveAtomsInCluster( clust, tmpat );

  // Prevent double counting
  std::vector<bool> incluster( getNumberOfNodes(), false );
  for(unsigned i=0;i<tmpat.size();++i) incluster[tmpat[i]]=true;

  // Find the atoms in the the clusters
  std::vector<bool> surface_atom( getNumberOfNodes(), false ); 
  for(unsigned i=0;i<tmpat.size();++i){
      for(unsigned j=0;j<getNumberOfNodes();++j){
         if( incluster[j] ) continue;
         double dist2=getSeparation( getPosition(tmpat[i]), getPosition(j) ).modulo2();
         if( dist2<rcut_surf2 ){ surface_atom[j]=true; }
      }
  }
  unsigned nsurf_at=0; 
  for(unsigned j=0;j<getNumberOfNodes();++j){
     if( surface_atom[j] ) nsurf_at++; 
  }
  myatoms.resize( nsurf_at + tmpat.size() );
  for(unsigned i=0;i<tmpat.size();++i) myatoms[i]=tmpat[i];
  unsigned nn=tmpat.size();
  for(unsigned j=0;j<getNumberOfNodes();++j){
      if( surface_atom[j] ){ myatoms[nn]=j; nn++; }
  }
  plumed_assert( nn==myatoms.size() );
}

Vector ClusterWithSurface::getPositionOfAtomForLinkCells( const unsigned& iatom ) const {
  return myclusters->getPositionOfAtomForLinkCells( iatom );
}

}
}
