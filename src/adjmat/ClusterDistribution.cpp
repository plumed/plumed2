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
#include "AdjacencyMatrixVessel.h"
#include "core/ActionRegister.h"
#include "tools/SwitchingFunction.h"

//+PLUMEDOC CONCOMP CLUSTER_DISTRIBUTION
/*
Calculate functions of the distribution of properties in your connected components.

\par Examples

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
  void performTask( const unsigned& , const unsigned& , MultiValue& ) const ;
};

PLUMED_REGISTER_ACTION(ClusterDistribution,"CLUSTER_DISTRIBUTION")

void ClusterDistribution::registerKeywords( Keywords& keys ){
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
   if( input!="none" ){
      use_switch=true; sf.set( input, errors );
      if( errors.length()!=0 ) error("problem reading SWITCH keyword : " + errors );
   }
   parseFlag("INVERSE_TRANSFORM",inverse);
   if( inverse && !use_switch ) error("INVERSE_TRANSFORM option was specified but no TRANSOFRM switching function was given");

   // Create all tasks by copying those from underlying DFS object (which is actually MultiColvar)
   for(unsigned i=0;i<getNumberOfNodes();++i) addTaskToList(i);

   // And now finish the setup of everything in the base
   std::vector<AtomNumber> fake_atoms; setupMultiColvarBase( fake_atoms );
}

void ClusterDistribution::calculate(){
   // Activate the relevant tasks
   nderivatives = getNumberOfDerivatives();
   deactivateAllTasks(); 
   for(unsigned i=0;i<getNumberOfClusters();++i) taskFlags[i]=1;
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
  for(unsigned j=0;j<myatoms.size();++j){
      unsigned i=myatoms[j];
      getPropertiesOfNode( i, vals );
      if( use_switch && !inverse ){
          vv = 1.0 - sf.calculate( vals[1], df );
          tval += vals[0]*vv; df=-df*vals[1];
      } else if( use_switch ){
          vv = sf.calculate( vals[1], df );
          tval += vals[0]*vv; df=df*vals[1];
      } else {
          tval += vals[0]*vals[1]; df=1.; vv=vals[1];
      }
      if( !doNotCalculateDerivatives() ){
          getNodePropertyDerivatives( i, tvals );
          for(unsigned k=0;k<tvals.getNumberActive();++k){
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
