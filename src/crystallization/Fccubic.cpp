/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2015 The plumed team
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
#include "multicolvar/MultiColvar.h"
#include "tools/NeighborList.h"
#include "core/ActionRegister.h"
#include "tools/SwitchingFunction.h"

#include <string>
#include <cmath>

using namespace std;

namespace PLMD{
namespace crystallization{

//+PLUMEDOC MCOLVAR FCCUBIC    
/*
*/
//+ENDPLUMEDOC


class Fccubic : public multicolvar::MultiColvar {
private:
//  double nl_cut;
  double rcut2;
  SwitchingFunction switchingFunction;
public:
  static void registerKeywords( Keywords& keys );
  Fccubic(const ActionOptions&);
// active methods:
  virtual double compute(); 
  Vector getCentralAtom();
/// Returns the number of coordinates of the field
  bool isPeriodic(){ return false; }
};

PLUMED_REGISTER_ACTION(Fccubic,"FCCUBIC")

void Fccubic::registerKeywords( Keywords& keys ){
  multicolvar::MultiColvar::registerKeywords( keys );
  keys.use("SPECIES"); keys.use("SPECIESA"); keys.use("SPECIESB");
  keys.add("compulsory","NN","6","The n parameter of the switching function ");
  keys.add("compulsory","MM","12","The m parameter of the switching function ");
  keys.add("compulsory","D_0","0.0","The d_0 parameter of the switching function");
  keys.add("compulsory","R_0","The r_0 parameter of the switching function");
  keys.add("optional","SWITCH","This keyword is used if you want to employ an alternative to the continuous swiching function defined above. "
                               "The following provides information on the \\ref switchingfunction that are available. "
                               "When this keyword is present you no longer need the NN, MM, D_0 and R_0 keywords.");
  // Use actionWithDistributionKeywords
  keys.use("MEAN"); keys.use("MORE_THAN"); keys.use("LESS_THAN"); keys.use("MAX");
  keys.use("MIN"); keys.use("BETWEEN"); keys.use("HISTOGRAM"); keys.use("MOMENTS");
}

Fccubic::Fccubic(const ActionOptions&ao):
PLUMED_MULTICOLVAR_INIT(ao)
{
  // Read in the switching function
  std::string sw, errors; parse("SWITCH",sw);
  if(sw.length()>0){
     switchingFunction.set(sw,errors);
     if( errors.length()!=0 ) error("problem reading SWITCH keyword : " + errors );
  } else { 
     double r_0=-1.0, d_0; int nn, mm;
     parse("NN",nn); parse("MM",mm);
     parse("R_0",r_0); parse("D_0",d_0);
     if( r_0<0.0 ) error("you must set a value for R_0");
     switchingFunction.set(nn,mm,r_0,d_0);
  }
  log.printf("  measure of simple cubicity around central atom.  Includes those atoms within %s\n",( switchingFunction.description() ).c_str() );
  // Set the link cell cutoff
  rcut2 = switchingFunction.get_dmax()*switchingFunction.get_dmax();
  setLinkCellCutoff( switchingFunction.get_dmax() );

  // Read in the atoms
  int natoms=2; readAtoms( natoms );
  // And setup the ActionWithVessel
  checkRead();
}

double Fccubic::compute(){
   weightHasDerivatives=true;
   double value=0, norm=0, dfunc; Vector distance;

   // Calculate the coordination number
   Vector myder, fder;
   double sw, t0, t1, t2, t3, x2, x4, y2, y4, z2, z4, r8, tmp;
   for(unsigned i=1;i<getNAtoms();++i){
      distance=getSeparation( getPosition(0), getPosition(i) );
      double d2 = distance.modulo2();
      if( d2<rcut2 ){ 
         sw = switchingFunction.calculateSqr( d2, dfunc ); 
   
         norm += sw;

         x2 = distance[0]*distance[0];
         x4 = x2*x2;

         y2 = distance[1]*distance[1];
         y4 = y2*y2;

         z2 = distance[2]*distance[2];
         z4 = z2*z2;
                 
         r8 = pow( distance.modulo2(), 4 );

         tmp = ((x4*y4)+(x4*z4)+(y4*z4))/r8;

         value += sw*tmp;

         t0 = (x2*y4+x2*z4)/r8;
         t1 = (y2*x4+y2*z4)/r8;
         t2 = (z2*x4+z2*y4)/r8;
         t3 = 2*tmp/distance.modulo2();         
 
         myder[0]=4*distance[0]*(t0-t3);
         myder[1]=4*distance[1]*(t1-t3);
         myder[2]=4*distance[2]*(t2-t3);
  
         fder = (+dfunc)*tmp*distance + sw*myder;

         addAtomsDerivatives( 0, -fder );
         addAtomsDerivatives( i, +fder );
         addBoxDerivatives( Tensor(distance,-fder) );
         addAtomsDerivativeOfWeight( 0, (-dfunc)*distance );
         addAtomsDerivativeOfWeight( i, (+dfunc)*distance );
         addBoxDerivativesOfWeight( (-dfunc)*Tensor(distance,distance) );
      }
   }
   
   setElementValue(0, value); setElementValue(1, norm ); 
   // values -> der of... value [0], weight[1], x coord [2], y, z... [more magic]
   updateActiveAtoms(); quotientRule( 0, 1, 0 ); clearDerivativesAfterTask(1);
   // Weight doesn't really have derivatives (just use the holder for convenience)
   weightHasDerivatives=false; setElementValue( 1, 1.0 );

   return value / norm; // this is equivalent to getting an "atomic" CV
}

Vector Fccubic::getCentralAtom(){
   addCentralAtomDerivatives( 0, Tensor::identity() );
   return getPosition(0);
}

}
}

