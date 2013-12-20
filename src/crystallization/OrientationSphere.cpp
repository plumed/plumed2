/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013 The plumed team
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
#include "VectorMultiColvar.h"
#include "OrientationSphere.h"

using namespace std;

namespace PLMD{
namespace crystallization {

void OrientationSphere::registerKeywords( Keywords& keys ){
  multicolvar::MultiColvarFunction::registerKeywords( keys );
  keys.add("compulsory","NN","6","The n parameter of the switching function ");
  keys.add("compulsory","MM","12","The m parameter of the switching function ");
  keys.add("compulsory","D_0","0.0","The d_0 parameter of the switching function");
  keys.add("compulsory","R_0","The r_0 parameter of the switching function");
  keys.add("optional","SWITCH","This keyword is used if you want to employ an alternative to the continuous swiching function defined above. "
                               "The following provides information on the \\ref switchingfunction that are available. "
                               "When this keyword is present you no longer need the NN, MM, D_0 and R_0 keywords.");
  // Use actionWithDistributionKeywords
  keys.use("MEAN"); keys.use("MORE_THAN"); keys.use("LESS_THAN"); 
  keys.use("MIN"); keys.use("BETWEEN"); keys.use("HISTOGRAM"); keys.use("MOMENTS");
}

OrientationSphere::OrientationSphere(const ActionOptions&ao):
Action(ao),
MultiColvarFunction(ao)
{
  // Resize everything that stores a vector now that we know the 
  // number of components
  unsigned ncomponents=getBaseMultiColvar(0)->getNumberOfQuantities() - 5;
  catom_orient.resize( ncomponents ); 
  catom_der.resize( ncomponents );
  this_orient.resize( ncomponents ); 

  // Weight of this does not have derivatives
  weightHasDerivatives=false;
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

  // Finish the setup of the object
  buildSymmetryFunctionLists( true );

  // And check everything has been read in correctly
  checkRead();
}

void OrientationSphere::calculateWeight(){
  weightHasDerivatives=true;   // The weight has no derivatives really
  setElementValue(1,1.0);
} 

double OrientationSphere::compute(){
   // Make sure derivatives for central atom are only calculated once
   VectorMultiColvar* vv = dynamic_cast<VectorMultiColvar*>( getBaseMultiColvar(0) );
   vv->firstcall=true;

   weightHasDerivatives=true;   // The weight has no derivatives really
   double sw, value=0, denom=0, dot, f_dot, dot_df, dfunc; Vector distance;

   getValueForBaseTask(0, catom_orient );
   for(unsigned i=1;i<getNAtoms();++i){
      distance=getSeparation( getPositionOfCentralAtom(0), getPositionOfCentralAtom(i) );
      sw = switchingFunction.calculateSqr( distance.modulo2(), dfunc );
      if( sw>=getTolerance() ){    
         getValueForBaseTask( i, this_orient );
         // Calculate the dot product wrt to this position 
         dot=0; for(unsigned k=0;k<catom_orient.size();++k) dot+=catom_orient[k]*this_orient[k];  
         f_dot = transformDotProduct( dot, dot_df ); 
         // N.B. We are assuming here that the imaginary part of the dot product is zero
         for(unsigned k=0;k<catom_orient.size();++k){
            this_orient[k]*=sw*dot_df; catom_der[k]=sw*dot_df*catom_orient[k];
         }  

         // Set the derivatives wrt of the numerator
         addOrientationDerivatives( 0, this_orient ); 
         addOrientationDerivatives( i, catom_der );  
         addCentralAtomsDerivatives( 0, 0, f_dot*(-dfunc)*distance );
         addCentralAtomsDerivatives( i, 0, f_dot*(dfunc)*distance );
         addBoxDerivatives( f_dot*(-dfunc)*Tensor(distance,distance) );
         value += sw*f_dot;
         // Set the derivatives wrt to the numerator
         addCentralAtomsDerivatives( 0, 1, (-dfunc)*distance );
         addCentralAtomsDerivatives( i, 1, (dfunc)*distance );
         addBoxDerivativesOfWeight( (-dfunc)*Tensor(distance,distance) );
         denom += sw;
      } else {
         removeAtomRequest( i, sw );   
      }
   }
   
   // Now divide everything
   unsigned nder = getNumberOfDerivatives();
   for(unsigned i=0;i<nder;++i){
      setElementDerivative( i, getElementDerivative(i)/denom - (value*getElementDerivative(nder+i))/(denom*denom) );  
      setElementDerivative( nder + i, 0.0 );
   }
   weightHasDerivatives=false;   // Weight has no derivatives we just use the holder for weight to store some stuff
   return value / denom;
}

Vector OrientationSphere::getCentralAtom(){
   addDerivativeOfCentralAtomPos( 0, Tensor::identity() );
   return getPositionOfCentralAtom(0);
}

}
}

