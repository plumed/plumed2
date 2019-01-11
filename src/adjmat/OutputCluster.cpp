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
#include "tools/OFile.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/ActionPilot.h"
#include "core/ActionRegister.h"

//+PLUMEDOC CONCOMP OUTPUT_CLUSTER
/*
Output the indices of the atoms in one of the clusters identified by a clustering object

This action provides one way of getting output from a \ref DFSCLUSTERING calculation.
The output in question here is either

- a file that contains a list of the atom indices that form part of one of the clusters that was identified using \ref DFSCLUSTERING
- an xyz file containing the positions of the atoms in one of the the clusters that was identified using \ref DFSCLUSTERING

Notice also that if you choose to output an xyz file you can ask PLUMED to try to reconstruct the cluster
taking the periodic boundary conditions into account by using the MAKE_WHOLE flag.

\par Examples

The input shown below identifies those atoms with a coordination number less than 13
and then constructs a contact matrix that describes the connectivity between the atoms
that satisfy this criteria.  The DFS algorithm is then used to find the connected components
in this matrix and the indices of the atoms in the largest connected component are then output
to a file.

\plumedfile
c1: COORDINATIONNUMBER SPECIES=1-1996 SWITCH={CUBIC D_0=0.34 D_MAX=0.38}
cf: MFILTER_LESS DATA=c1 SWITCH={CUBIC D_0=13 D_MAX=13.5}
mat: CONTACT_MATRIX ATOMS=cf SWITCH={CUBIC D_0=0.34 D_MAX=0.38}
dfs: DFSCLUSTERING MATRIX=mat
OUTPUT_CLUSTER CLUSTERS=dfs CLUSTER=1 FILE=dfs.dat
\endplumedfile

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace adjmat {

class OutputCluster : public ActionPilot {
private:
  bool makewhole, output_xyz;
  OFile ofile;
  ClusteringBase* myclusters;
  double rcut2;
  unsigned clustr, maxdepth, maxgoes;
  std::vector<bool> visited;
  std::vector<unsigned> myatoms;
  std::vector<Vector> atomsin;
  std::vector<unsigned> nneigh;
  Matrix<unsigned> adj_list;
  int number_of_cluster;
  std::vector< std::pair<unsigned,unsigned> > cluster_sizes;
  std::vector<unsigned> which_cluster;
  bool explore_dfs( const unsigned& index );
  void explore( const unsigned& index, const unsigned& depth );
public:
  static void registerKeywords( Keywords& keys );
  explicit OutputCluster(const ActionOptions&);
  void calculate() {}
  void apply() {}
  void update();
};

PLUMED_REGISTER_ACTION(OutputCluster,"OUTPUT_CLUSTER")

void OutputCluster::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionPilot::registerKeywords( keys );
  keys.add("compulsory","CLUSTERS","the action that performed the clustering");
  keys.add("compulsory","CLUSTER","1","which cluster would you like to look at 1 is the largest cluster, 2 is the second largest, 3 is the the third largest and so on");
  keys.add("compulsory","STRIDE","1","the frequency with which you would like to output the atoms in the cluster");
  keys.add("compulsory","FILE","the name of the file on which to output the details of the cluster");
  keys.add("compulsory","MAXDEPTH","6","maximum depth for searches over paths to reconstruct clusters for PBC");
  keys.add("compulsory","MAXGOES","200","number of times to run searches to reconstruct clusters");
  keys.addFlag("MAKE_WHOLE",false,"reconstruct the clusters and remove all periodic boundary conditions.");
}

OutputCluster::OutputCluster(const ActionOptions& ao):
  Action(ao),
  ActionPilot(ao),
  myclusters(NULL)
{
  // Setup output file
  ofile.link(*this); std::string file; parse("FILE",file);
  if( file.length()==0 ) error("output file name was not specified");
  // Search for xyz extension
  output_xyz=false;
  if( file.find(".")!=std::string::npos ) {
    std::size_t dot=file.find_first_of('.');
    if( file.substr(dot+1)=="xyz" ) output_xyz=true;
  }

  ofile.open(file); log.printf("  on file %s \n",file.c_str());
  parseFlag("MAKE_WHOLE",makewhole); parse("MAXDEPTH",maxdepth); parse("MAXGOES",maxgoes);
  if( makewhole && !output_xyz) error("MAKE_WHOLE flag is not compatible with output of non-xyz files");

  // Find what action we are taking the clusters from
  std::vector<std::string> matname(1); parse("CLUSTERS",matname[0]);
  myclusters = plumed.getActionSet().selectWithLabel<ClusteringBase*>( matname[0] );
  if( !myclusters ) error( matname[0] + " does not calculate perform a clustering of the atomic positions");
  // N.B. the +0.3 is a fudge factor.  Reconstrucing PBC doesnt work without this GAT
  addDependency( myclusters ); double rcut=myclusters->getCutoffForConnection() + 0.3; rcut2=rcut*rcut;

  // Read in the cluster we are calculating
  parse("CLUSTER",clustr);
  if( clustr<1 ) error("cannot look for a cluster larger than the largest cluster");
  if( clustr>myclusters->getNumberOfNodes() ) error("cluster selected is invalid - too few atoms in system");
  log.printf("  outputting atoms in %u th largest cluster found by %s \n",clustr,matname[0].c_str() );
}

void OutputCluster::update() {
  myclusters->retrieveAtomsInCluster( clustr, myatoms );
  if( output_xyz ) {
    ofile.printf("%u \n",static_cast<unsigned>(myatoms.size()));
    ofile.printf("atoms in %u th largest cluster \n",clustr );
    if( makewhole ) {
      // Retrieve the atom positions
      atomsin.resize( myatoms.size() );
      for(unsigned i=0; i<myatoms.size(); ++i) atomsin[i]=myclusters->getPositionOfAtomForLinkCells( myatoms[i] );
      // Build a connectivity matrix neglecting the pbc
      nneigh.resize( myatoms.size(), 0 ); adj_list.resize( myatoms.size(), myatoms.size() );
      for(unsigned i=1; i<myatoms.size(); ++i) {
        for(unsigned j=0; j<i; ++j) {
          if( delta( atomsin[i], atomsin[j] ).modulo2()<=rcut2 ) { adj_list(i,nneigh[i])=j; adj_list(j,nneigh[j])=i; nneigh[i]++; nneigh[j]++; }
        }
      }
      // Use DFS to find the largest cluster not broken by PBC
      number_of_cluster=-1; visited.resize( myatoms.size(), false );
      cluster_sizes.resize( myatoms.size() ); which_cluster.resize( myatoms.size() );
      for(unsigned i=0; i<cluster_sizes.size(); ++i) { cluster_sizes[i].first=0; cluster_sizes[i].second=i; }

      for(unsigned i=0; i<myatoms.size(); ++i) {
        if( !visited[i] ) { number_of_cluster++; visited[i]=explore_dfs(i); }
      }
      std::sort( cluster_sizes.begin(), cluster_sizes.end() );

      // Now set visited so that only those atoms in largest cluster will be start points for PBCing
      visited.assign( visited.size(), false );
      for(unsigned i=0; i<myatoms.size(); ++i) {
        if( which_cluster[i]==cluster_sizes[cluster_sizes.size()-1].second ) visited[i]=true;
      }

      // Now retrieve the original connectivity matrix (including pbc)
      nneigh.assign( nneigh.size(), 0 );
      for(unsigned i=1; i<myatoms.size(); ++i) {
        for(unsigned j=0; j<i; ++j) {
          if( myclusters->areConnected( myatoms[i], myatoms[j] ) ) { adj_list(i,nneigh[i])=j; adj_list(j,nneigh[j])=i; nneigh[i]++; nneigh[j]++; }
        }
      }

      // Now find broken bonds and run iterative deepening depth first search to reconstruct
      for(unsigned jj=0; jj<maxgoes; ++jj) {

        for(unsigned j=0; j<myatoms.size(); ++j) {
          if( !visited[j] ) continue;

          for(unsigned k=0; k<nneigh[j]; ++k) {
            if( delta( atomsin[j],atomsin[adj_list(j,k)] ).modulo2()>rcut2 ) {
              visited[j]=true;
              for(unsigned depth=0; depth<=maxdepth; ++depth) explore( j, depth );
            }
          }
        }
      }
      // And print final positions
      for(unsigned i=0; i<myatoms.size(); ++i) ofile.printf( "X %f %f %f \n", atomsin[i][0], atomsin[i][1], atomsin[i][2] );
    } else {
      for(unsigned i=0; i<myatoms.size(); ++i) {
        Vector pos=myclusters->getPositionOfAtomForLinkCells( myatoms[i] );
        ofile.printf( "X %f %f %f \n", pos[0], pos[1], pos[2] );
      }
    }
  } else {
    ofile.printf("CLUSTERING RESULTS AT TIME %f : NUMBER OF ATOMS IN %u TH LARGEST CLUSTER EQUALS %u \n",getTime(),clustr,static_cast<unsigned>(myatoms.size()) );
    ofile.printf("INDICES OF ATOMS : ");
    for(unsigned i=0; i<myatoms.size(); ++i) ofile.printf("%d ",(myclusters->getAbsoluteIndexOfCentralAtom(myatoms[i])).index());
    ofile.printf("\n");
  }
}

void OutputCluster::explore( const unsigned& index, const unsigned& depth ) {
  if( depth==0 ) return ;

  for(unsigned i=0; i<nneigh[index]; ++i) {
    unsigned j=adj_list(index,i); visited[j]=true;
    Vector svec=myclusters->pbcDistance( atomsin[index], atomsin[j] );
    atomsin[j] = atomsin[index] + svec;
    explore( j, depth-1 );
  }
}

bool OutputCluster::explore_dfs( const unsigned& index ) {
  visited[index]=true;
  for(unsigned i=0; i<nneigh[index]; ++i) {
    unsigned j=adj_list(index,i);
    if( !visited[j] ) visited[j]=explore_dfs(j);
  }

  // Count the size of the cluster
  cluster_sizes[number_of_cluster].first++;
  which_cluster[index] = number_of_cluster;
  return visited[index];
}

}
}


