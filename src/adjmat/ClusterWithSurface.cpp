/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2019 The plumed team
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
Take a connected component that was found using a clustering algorithm and create a new cluster that contains those atoms that are in the cluster together with those atoms that are within a certain cutoff of the cluster.

As discussed in the section of the manual on \ref contactmatrix a useful tool for developing complex collective variables is the notion of the
so called adjacency matrix.  An adjacency matrix is an \f$N \times N\f$ matrix in which the \f$i\f$th, \f$j\f$th element tells you whether
or not the \f$i\f$th and \f$j\f$th atoms/molecules from a set of \f$N\f$ atoms/molecules are adjacent or not.  When analyzing these matrix
we can treat them as a graph and find connected components using some clustering algorithm.  This action is used in tandem with this form of analysis
and takes one of the connected components that was found during this analysis and creates a new cluster that includes all the atoms within the
connected component that was found together that were within a certain cutoff distance of the atoms in the connected component.  This form of analysis
has been used successfully in the forward flux sampling simulations described in this paper \cite gab-ice-kaolinite

\par Examples

The following input uses PLUMED to calculate a adjacency matrix that connects a pair of atoms if they both have a coordination number that is less
than 13.5 and if they are within 0.38 nm of each other.  Depth first search clustering is used to find the connected components in this matrix.  The
number of atoms with indices that are between 1 and 1996 and that are either in the second largest cluster or that are within within 0.3 nm of one of the
atoms within the the second largest cluster are then counted and this number of atoms is output to a file called size.  In addition the indices of the atoms
that were counted are output to a file called dfs2.dat.

\plumedfile
c1: COORDINATIONNUMBER SPECIES=1-1996 SWITCH={CUBIC D_0=0.34 D_MAX=0.38}
cf: MFILTER_LESS DATA=c1 SWITCH={CUBIC D_0=13 D_MAX=13.5}
mat: CONTACT_MATRIX ATOMS=cf SWITCH={CUBIC D_0=0.34 D_MAX=0.38}
dfs: DFSCLUSTERING MATRIX=mat
clust2a: CLUSTER_WITHSURFACE CLUSTERS=dfs RCUT_SURF=0.3
size2a: CLUSTER_NATOMS CLUSTERS=clust2a CLUSTER=2
PRINT ARG=size2a FILE=size FMT=%8.4f
OUTPUT_CLUSTER CLUSTERS=clust2a CLUSTER=2 FILE=dfs2.dat
\endplumedfile


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
  unsigned getNumberOfDerivatives() override;
///
  unsigned getNumberOfNodes() const override;
///
  AtomNumber getAbsoluteIndexOfCentralAtom(const unsigned& i) const override;
///
  void retrieveAtomsInCluster( const unsigned& clust, std::vector<unsigned>& myatoms ) const override;
///
  void getInputData( const unsigned& ind, const bool& normed, const multicolvar::AtomValuePack& myatoms, std::vector<double>& orient0 ) const override;
///
  MultiValue& getInputDerivatives( const unsigned& ind, const bool& normed, const multicolvar::AtomValuePack& myatoms ) const override;
///
  unsigned getNumberOfQuantities() const override;
/// Do the calculation
  void performClustering() override {};
///
  double  getCutoffForConnection() const override;
///
  Vector getPositionOfAtomForLinkCells( const unsigned& taskIndex ) const override;
};

PLUMED_REGISTER_ACTION(ClusterWithSurface,"CLUSTER_WITHSURFACE")

void ClusterWithSurface::registerKeywords( Keywords& keys ) {
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

unsigned ClusterWithSurface::getNumberOfDerivatives() {
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
  for(unsigned i=0; i<tmpat.size(); ++i) incluster[tmpat[i]]=true;

  // Find the atoms in the the clusters
  std::vector<bool> surface_atom( getNumberOfNodes(), false );
  for(unsigned i=0; i<tmpat.size(); ++i) {
    for(unsigned j=0; j<getNumberOfNodes(); ++j) {
      if( incluster[j] ) continue;
      double dist2=getSeparation( getPosition(tmpat[i]), getPosition(j) ).modulo2();
      if( dist2<rcut_surf2 ) { surface_atom[j]=true; }
    }
  }
  unsigned nsurf_at=0;
  for(unsigned j=0; j<getNumberOfNodes(); ++j) {
    if( surface_atom[j] ) nsurf_at++;
  }
  myatoms.resize( nsurf_at + tmpat.size() );
  for(unsigned i=0; i<tmpat.size(); ++i) myatoms[i]=tmpat[i];
  unsigned nn=tmpat.size();
  for(unsigned j=0; j<getNumberOfNodes(); ++j) {
    if( surface_atom[j] ) { myatoms[nn]=j; nn++; }
  }
  plumed_assert( nn==myatoms.size() );
}

Vector ClusterWithSurface::getPositionOfAtomForLinkCells( const unsigned& iatom ) const {
  return myclusters->getPositionOfAtomForLinkCells( iatom );
}

}
}
