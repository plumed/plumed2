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

//+PLUMEDOC MCOLVARF DFSMAXCLUSTER
/*
Find the connected components and determine the maximum size of the clusters in your system using a smooth function

This action uses the DFS clustering algorithm described in \ref DFSCLUSTERING to find a set of connected components
based on the configuration of the atoms in your system.  Once again this can be used to find crystalline nuclei or 
bubble of atoms.  Once these connected components have been identified this colvar attempts to determine the size
of the largest cluster of atoms using a continuous function.  This is done using the following formula
\f[
z = \beta \ln \left[ \sum_j \exp\left( \frac{\sum_i \sigma(s_{ij} }{\beta} \right) \right]
\f] 
Here the sum over \f$j\f$ runs over the set of connected components that are identified through the DFS clustering.
Meanwhile, the sum over \f$i\f$ runs over the set of symmetry functions, \f$s_{ij}\f$, that are contained within the \f$j\f$th connected
component.  The function \f$\sigma\f$ is a switching function that is specified using the TRANSFORM keyword and used to 
convert all symmetry functions to values between 0 and 1.  This is done in order to make it easier to interpret the final
quantity as the size of the largest cluster in the system.

\par Examples

The following example calculates the size of the largest cluster of the system.  FCCCUBIC parameters, which don't 
necessarily have to be equal to numbers between zero and one, are converted to numbers between zero and one by 
virtue of the switching function specified using the TRANSFORM keyword.

\verbatim
cubic1: FCCUBIC SPECIES=1-1000 SWITCH={CUBIC D_0=0.4  D_MAX=0.5} TOL=0.03 
clust: DFSMAXCLUSTER DATA=cubic1 BETA=0.5 SWITCH={CUBIC D_0=0.4   D_MAX=0.5} TRANSFORM={CUBIC D_0=0.035 D_MAX=0.045}
\endverbatim

As was described in the examples section of the page on the \ref DFSCLUSTERING action we can also use filters when 
we use DFSMAXCLUSTER.  Here again the filter ensures that clustering is only performed for those atoms that have a 
symmetry function greater than a certain parameter.  As such the clustering will only be performed for a subset of the
symmetry functions claculated by the action labelled cubic1.

\verbatim
cubic1: FCCUBIC SPECIES=1-1000 SWITCH={CUBIC D_0=0.4  D_MAX=0.5} TOL=0.03 
cf: MFILTER_MORE DATA=cubic1 SWITCH={CUBIC D_0=0.035 D_MAX=0.045}
clust: DFSMAXCLUSTER DATA=cf BETA=0.5 SWITCH={CUBIC D_0=0.4   D_MAX=0.5} TRANSFORM={CUBIC D_0=0.035 D_MAX=0.045} WTOL=0.01
\endverbatim

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace crystallization {

class DFSMaxCluster : public DFSClustering {
private:
/// The value of beta
  double beta;
///
  bool use_switch, inverse;
//
  SwitchingFunction sf;
public:
/// Create manual
  static void registerKeywords( Keywords& keys );
/// Constructor
  explicit DFSMaxCluster(const ActionOptions&);
///
  void doCalculationOnCluster();
};

PLUMED_REGISTER_ACTION(DFSMaxCluster,"DFSMAXCLUSTER")

void DFSMaxCluster::registerKeywords( Keywords& keys ){
  DFSClustering::registerKeywords( keys );
  keys.add("compulsory","BETA","the value of beta to be used in calculating the smooth maximum");
  keys.use("WTOL"); keys.use("USE_ORIENTATION");
  keys.remove("LOWMEM"); keys.use("HIGHMEM");
  keys.add("compulsory","TRANSFORM","none","the switching function to use to convert the crystallinity parameter to a number between zero and one");
  keys.addFlag("INVERSE_TRANSFORM",false,"when TRANSFORM appears alone the input symmetry functions, \\fx\\f$ are transformed used \\f$1-s(x)\\f$ "
                                         "where \\f$s(x)\\f$ is a switching function.  When this option is used you instead transform using \\f$s(x)\\f$ only.");
}

DFSMaxCluster::DFSMaxCluster(const ActionOptions&ao):
Action(ao),
DFSClustering(ao)
{
   // Find out the value of beta
   parse("BETA",beta);
   addValueWithDerivatives(); setNotPeriodic();

   use_switch=false;
   std::string input, errors; parse("TRANSFORM",input);
   if( input!="none" ){
      use_switch=true; sf.set( input, errors );
      if( errors.length()!=0 ) error("problem reading SWITCH keyword : " + errors );
   }
   parseFlag("INVERSE_TRANSFORM",inverse);
   if( inverse && !use_switch ) error("INVERSE_TRANSFORM option was specified but no TRANSOFRM switching function was given");
}

void DFSMaxCluster::doCalculationOnCluster(){
   unsigned size=comm.Get_size(), rank=comm.Get_rank(); 

   std::vector<double> tder( getNumberOfDerivatives() ), fder( 1+getNumberOfDerivatives(), 0.0 );
   MultiValue myvals( getNumberOfQuantities(), getNumberOfDerivatives() );
   std::vector<unsigned> myatoms; std::vector<double> vals( getNumberOfQuantities() );

   fder[0]=0;  // Value and derivatives are accumulated in one array so that there is only one MPI call 
   for(unsigned iclust=0;iclust<getNumberOfClusters();++iclust){
       retrieveAtomsInCluster( iclust+1, myatoms );
       // This deals with filters
       if( myatoms.size()==1 && !isCurrentlyActive(0,myatoms[0]) ) continue ;

       double vv, df, tval=0; tder.assign( tder.size(), 0.0 );
       for(unsigned j=0;j<myatoms.size();++j){ 
           unsigned i=myatoms[j];
           getVectorForTask( i, false, vals );
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
               getVectorDerivatives( i, false, myvals );
               for(unsigned k=0;k<myvals.getNumberActive();++k){
                   unsigned kat=myvals.getActiveIndex(k);
                   tder[kat]+=vals[0]*df*myvals.getDerivative(1,kat) + vv*myvals.getDerivative(0,kat);
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

   // And finish the derivatives
   setValue( beta*std::log( fder[0] ) ); Value* myval=getPntrToValue();
   if( !doNotCalculateDerivatives() ){
      double pref = beta/fder[0];
      for(unsigned k=0;k<getNumberOfDerivatives();++k) myval->addDerivative( k, pref*fder[1+k] );
   }
}

}
}
