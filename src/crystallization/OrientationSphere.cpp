/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2016 The plumed team
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
#include "VectorMultiColvar.h"
#include "OrientationSphere.h"

using namespace std;

namespace PLMD{
namespace crystallization {

void OrientationSphere::registerKeywords( Keywords& keys ){
  multicolvar::MultiColvarFunction::registerKeywords( keys );
  keys.add("compulsory","NN","6","The n parameter of the switching function ");
  keys.add("compulsory","MM","0","The m parameter of the switching function; 0 implies 2*NN");
  keys.add("compulsory","D_0","0.0","The d_0 parameter of the switching function");
  keys.add("compulsory","R_0","The r_0 parameter of the switching function");
  keys.add("optional","SWITCH","This keyword is used if you want to employ an alternative to the continuous swiching function defined above. "
                               "The following provides information on the \\ref switchingfunction that are available. "
                               "When this keyword is present you no longer need the NN, MM, D_0 and R_0 keywords.");
  // Use actionWithDistributionKeywords
  keys.use("MEAN"); keys.use("MORE_THAN"); keys.use("LESS_THAN"); 
  keys.use("MIN"); keys.use("BETWEEN"); keys.use("HISTOGRAM"); keys.use("MOMENTS");
  keys.use("LOWEST"); keys.use("HIGHEST");
}

OrientationSphere::OrientationSphere(const ActionOptions&ao):
Action(ao),
MultiColvarFunction(ao)
{
  if( getNumberOfBaseMultiColvars()>1 ) warning("not sure if orientation sphere works with more than one base multicolvar - check numerical derivatives");
  // Read in the switching function
  std::string sw, errors; parse("SWITCH",sw);
  if(sw.length()>0){
     switchingFunction.set(sw,errors);
  } else { 
     double r_0=-1.0, d_0; int nn, mm;
     parse("NN",nn); parse("MM",mm);
     parse("R_0",r_0); parse("D_0",d_0);
     if( r_0<0.0 ) error("you must set a value for R_0");
     switchingFunction.set(nn,mm,r_0,d_0);
  }
  log.printf("  degree of overlap in orientation between central molecule and those within %s\n",( switchingFunction.description() ).c_str() );
  // Set the link cell cutoff
  rcut2 = switchingFunction.get_dmax()*switchingFunction.get_dmax();
  setLinkCellCutoff( switchingFunction.get_dmax() );

  // Finish the setup of the object
  buildSymmetryFunctionLists();
}

double OrientationSphere::compute( const unsigned& tindex, multicolvar::AtomValuePack& myatoms ) const {
   double d2, sw, value=0, denom=0, dot, f_dot, dot_df, dfunc; 
   unsigned ncomponents=getBaseMultiColvar(0)->getNumberOfQuantities();
   std::vector<double> catom_orient( ncomponents ), this_orient( ncomponents ), catom_der( ncomponents ); 

   Vector catom_pos = myatoms.getPosition(0);
   getVectorForTask( myatoms.getIndex(0), true, catom_orient ); 
   multicolvar::CatomPack atom0; 
   MultiValue& myder0=getVectorDerivatives( myatoms.getIndex(0), true ); 
   if( !doNotCalculateDerivatives() ) atom0=getCentralAtomPackFromInput( myatoms.getIndex(0) );

   for(unsigned i=1;i<myatoms.getNumberOfAtoms();++i){
      Vector& distance=myatoms.getPosition(i);  
      if ( (d2=distance[0]*distance[0])<rcut2 &&
           (d2+=distance[1]*distance[1])<rcut2 &&
           (d2+=distance[2]*distance[2])<rcut2) {
 
         sw = switchingFunction.calculateSqr( d2, dfunc );  
 
         getVectorForTask( myatoms.getIndex(i), true, this_orient );
         // Calculate the dot product wrt to this position 
         dot=0; for(unsigned k=2;k<catom_orient.size();++k) dot+=catom_orient[k]*this_orient[k];  
         f_dot = transformDotProduct( dot, dot_df ); 

         if( !doNotCalculateDerivatives() ){
             // N.B. We are assuming here that the imaginary part of the dot product is zero
             for(unsigned k=2;k<catom_orient.size();++k){
                 this_orient[k]*=sw*dot_df; catom_der[k]=sw*dot_df*catom_orient[k];
             }
             MultiValue& myder1=getVectorDerivatives( myatoms.getIndex(i), true );
             mergeVectorDerivatives( 1, 2, this_orient.size(), myatoms.getIndex(0), this_orient, myder0, myatoms );  
             mergeVectorDerivatives( 1, 2, catom_der.size(), myatoms.getIndex(i), catom_der, myder1, myatoms );
             myatoms.addComDerivatives( 1, f_dot*(-dfunc)*distance, atom0 );
             addAtomDerivatives( 1, i, f_dot*(dfunc)*distance, myatoms );
             myatoms.addBoxDerivatives( 1, (-dfunc)*f_dot*Tensor(distance,distance) );
             myder1.clearAll();
              
             myatoms.addComDerivatives( -1, (-dfunc)*distance, atom0 );
             addAtomDerivatives( -1, i, (dfunc)*distance, myatoms );
             myatoms.addTemporyBoxDerivatives( (-dfunc)*Tensor(distance,distance) );

         }
         value += sw*f_dot;
         denom += sw;
      }
   }
   double rdenom, df2, pref=calculateCoordinationPrefactor( denom, df2 );
   if( fabs(denom)>epsilon ){ rdenom = 1.0 / denom; }
   else { plumed_assert(fabs(value)<epsilon); rdenom=1.0; } 
  
   // Now divide everything
   double rdenom2=rdenom*rdenom;
   updateActiveAtoms( myatoms ); MultiValue& myvals=myatoms.getUnderlyingMultiValue();
   for(unsigned i=0;i<myvals.getNumberActive();++i){
       unsigned ider=myvals.getActiveIndex(i);
       double  dgd=myvals.getTemporyDerivative(ider);
       myvals.setDerivative( 1, ider, rdenom*(pref*myvals.getDerivative(1,ider)+value*df2*dgd) - (value*pref*dgd)*rdenom2 );
   } 

   return pref*rdenom*value;
}

}
}

