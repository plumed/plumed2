/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015 The plumed team
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

//+PLUMEDOC MCOLVARF DFSNUMEROFCLUSTERS
/*
Find the number of clusters that have a size larger than a certain cutoff

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace crystallization {

class DFSNumberOfClusters : public DFSClustering {
private:
  SwitchingFunction thresh;
///
  bool use_switch, inverse;
//
  SwitchingFunction sf, tsf;
public:
/// Create manual
  static void registerKeywords( Keywords& keys );
/// Constructor
  explicit DFSNumberOfClusters(const ActionOptions&);
///
  void doCalculationOnCluster();
};

PLUMED_REGISTER_ACTION(DFSNumberOfClusters,"DFSNUMEROFCLUSTERS")

void DFSNumberOfClusters::registerKeywords( Keywords& keys ){
  DFSClustering::registerKeywords( keys );
  keys.add("compulsory","TRANSFORM","none","the switching function to use to convert the crystallinity parameter to a number between zero and one");
  keys.add("compulsory","THRESHOLD","a switching function that defines how large the clusters should be in order to be counted");
  keys.use("WTOL"); keys.use("USE_ORIENTATION");
  keys.remove("LOWMEM"); keys.use("HIGHMEM");
  keys.addFlag("INVERSE_TRANSFORM",false,"when TRANSFORM appears alone the input symmetry functions, \\fx\\f$ are transformed used \\f$1-s(x)\\f$ "
                                         "where \\f$s(x)\\f$ is a switching function.  When this option is used you instead transform using \\f$s(x)\\f$ only.");
}

DFSNumberOfClusters::DFSNumberOfClusters(const ActionOptions&ao):
Action(ao),
DFSClustering(ao)
{
   // Find out the value of beta
   addValueWithDerivatives(); setNotPeriodic();

   use_switch=false;
   std::string input, errors; parse("TRANSFORM",input);
   if( input!="none" ){
      use_switch=true; sf.set( input, errors );
      if( errors.length()!=0 ) error("problem reading TRANSFORM keyword : " + errors );
   }
   parseFlag("INVERSE_TRANSFORM",inverse);
   if( inverse && !use_switch ) error("INVERSE_TRANSFORM option was specified but no TRANSOFRM switching function was given");

   parse("THRESHOLD",input);
   if( input.length()==0 ) error("found no input for threshold switching function");
   tsf.set( input, errors );
   if( errors.length()!=0 ) error("problem reading THRESHOLD keyword : " + errors );
}

void DFSNumberOfClusters::doCalculationOnCluster(){
   unsigned size=comm.Get_size(), rank=comm.Get_rank(); 

   std::vector<double> tder( getNumberOfDerivatives() );
   MultiValue myvals( getNumberOfQuantities(), getNumberOfDerivatives() );
   std::vector<unsigned> myatoms; std::vector<double> vals( getNumberOfQuantities() );

   double fder=0, dtdf;  Value* myval=getPntrToValue(); 
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
       fder += 1.0 - tsf.calculate( tval, dtdf ); double pref = -tval*dtdf; 
       if( !doNotCalculateDerivatives() ){
           for(unsigned k=0;k<getNumberOfDerivatives();++k) myval->addDerivative( k, pref*tder[k] );
       }
   }
   // And finish the derivatives
   setValue( fder ); 
}

}
}
