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
#include "ClusterAnalysisBase.h"
#include "AdjacencyMatrixVessel.h"
#include "core/ActionRegister.h"
#include "tools/SwitchingFunction.h"

//+PLUMEDOC CONCOMP CLUSTER_DISTRIBUTION
/*
Calculate functions of the distribution of properties in your connected components.

This collective variable was developed for looking at nucleation phenomena, where you are
interested in using studying the behavior of atoms in small aggregates or nuclei.  In these sorts of
problems you might be interested in the distribution of the sizes of the clusters in your system.
A detailed description of this CV can be found in \cite tribello-clustering.

\par Examples

The input provided below calculates the local q6 Steinhardt parameter on each atom.  The coordination number
that atoms with a high value for the local q6 Steinhardt parameter have with other atoms that have a high
value for the local q6 Steinhardt parameter is then computed.  A contact matrix is then computed that measures
whether atoms atoms \f$i\f$ and \f$j\f$ have a high value for this coordination number and if they are within
3.6 nm of each other.  The connected components of this matrix are then found using a depth first clustering
algorithm on the corresponding graph. The number of components in this graph that contain more than 27 atoms is then computed.
As discussed in \cite tribello-clustering this input was used to analyze the formation of a polycrystal of GeTe from amorphous
GeTe.

\plumedfile
q6: Q6 SPECIES=1-32768 SWITCH={GAUSSIAN D_0=5.29 R_0=0.01 D_MAX=5.3} LOWMEM
lq6: LOCAL_Q6 SPECIES=q6 SWITCH={GAUSSIAN D_0=5.29 R_0=0.01 D_MAX=5.3} LOWMEM
flq6: MFILTER_MORE DATA=lq6 SWITCH={GAUSSIAN D_0=0.19 R_0=0.01 D_MAX=0.2}
cc: COORDINATIONNUMBER SPECIES=flq6 SWITCH={GAUSSIAN D_0=3.59 R_0=0.01 D_MAX=3.6}
fcc: MFILTER_MORE DATA=cc SWITCH={GAUSSIAN D_0=5.99 R_0=0.01 D_MAX=6.0}
mat: CONTACT_MATRIX ATOMS=fcc SWITCH={GAUSSIAN D_0=3.59 R_0=0.01 D_MAX=3.6}
dfs: DFSCLUSTERING MATRIX=mat
nclust: CLUSTER_DISTRIBUTION CLUSTERS=dfs TRANSFORM={GAUSSIAN D_0=5.99 R_0=0.01 D_MAX=6.0} MORE_THAN={GAUSSIAN D_0=26.99 R_0=0.01 D_MAX=27}
PRINT ARG=nclust.* FILE=colvar
\endplumedfile

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace adjmat {

class ClusterDistribution : public ClusterAnalysisBase {
private:
  unsigned nderivatives;
///
  bool use_switch, inverse;
//
  SwitchingFunction sf;
public:
/// Create manual
  static void registerKeywords( Keywords& keys );
/// Constructor
  explicit ClusterDistribution(const ActionOptions&);
/// Do the calculation
  void calculate();
/// We can use ActionWithVessel to run all the calculation
  void performTask( const unsigned&, const unsigned&, MultiValue& ) const ;
};

PLUMED_REGISTER_ACTION(ClusterDistribution,"CLUSTER_DISTRIBUTION")

void ClusterDistribution::registerKeywords( Keywords& keys ) {
  ClusterAnalysisBase::registerKeywords( keys );
  keys.add("compulsory","TRANSFORM","none","the switching function to use to convert the crystallinity parameter to a number between zero and one");
  keys.addFlag("INVERSE_TRANSFORM",false,"when TRANSFORM appears alone the input symmetry functions, \\f$x\\f$ are transformed used \\f$1-s(x)\\f$ "
               "where \\f$s(x)\\f$ is a switching function.  When this option is used you instead transform using \\f$s(x)\\f$ only.");
  keys.use("MORE_THAN"); keys.use("LESS_THAN"); keys.use("BETWEEN");
  keys.use("HISTOGRAM"); keys.use("ALT_MIN"); keys.use("MIN"); keys.use("MAX");
}

ClusterDistribution::ClusterDistribution(const ActionOptions&ao):
  Action(ao),
  ClusterAnalysisBase(ao),
  nderivatives(0)
{
  use_switch=false;
  std::string input, errors; parse("TRANSFORM",input);
  if( input!="none" ) {
    use_switch=true; sf.set( input, errors );
    if( errors.length()!=0 ) error("problem reading SWITCH keyword : " + errors );
  }
  parseFlag("INVERSE_TRANSFORM",inverse);
  if( inverse && !use_switch ) error("INVERSE_TRANSFORM option was specified but no TRANSOFRM switching function was given");

  // Create all tasks by copying those from underlying DFS object (which is actually MultiColvar)
  for(unsigned i=0; i<getNumberOfNodes(); ++i) addTaskToList(i);

  // And now finish the setup of everything in the base
  std::vector<AtomNumber> fake_atoms; setupMultiColvarBase( fake_atoms );
}

void ClusterDistribution::calculate() {
  // Activate the relevant tasks
  nderivatives = getNumberOfDerivatives();
  deactivateAllTasks();
  for(unsigned i=0; i<getNumberOfClusters(); ++i) taskFlags[i]=1;
  lockContributors();
  // Now do the calculation
  runAllTasks();
}

void ClusterDistribution::performTask( const unsigned& task_index, const unsigned& current, MultiValue& myvals ) const {
  std::vector<unsigned> myatoms; retrieveAtomsInCluster( current+1, myatoms );
  // This deals with filters
  if( myatoms.size()==1 && !nodeIsActive(myatoms[0]) ) return ;

  std::vector<double> vals( getNumberOfQuantities() );
  MultiValue tvals( getNumberOfQuantities(), nderivatives );

  // And this builds everything for this particular atom
  double vv, df, tval=0;
  for(unsigned j=0; j<myatoms.size(); ++j) {
    unsigned i=myatoms[j];
    getPropertiesOfNode( i, vals );
    if( use_switch && !inverse ) {
      vv = 1.0 - sf.calculate( vals[1], df );
      tval += vals[0]*vv; df=-df*vals[1];
    } else if( use_switch ) {
      vv = sf.calculate( vals[1], df );
      tval += vals[0]*vv; df=df*vals[1];
    } else {
      tval += vals[0]*vals[1]; df=1.; vv=vals[1];
    }
    if( !doNotCalculateDerivatives() ) {
      getNodePropertyDerivatives( i, tvals );
      for(unsigned k=0; k<tvals.getNumberActive(); ++k) {
        unsigned kat=tvals.getActiveIndex(k);
        myvals.addDerivative( 1, kat, vals[0]*df*tvals.getDerivative(1,kat) + vv*tvals.getDerivative(0,kat) );
      }
      tvals.clearAll();
    }
  }
  myvals.setValue( 0, 1.0 ); myvals.setValue( 1, tval );
}

}
}
