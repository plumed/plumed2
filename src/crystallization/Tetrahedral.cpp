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

//+PLUMEDOC MCOLVAR TETRAHEDRAL
/*

\par Examples

*/
//+ENDPLUMEDOC


class Tetrahedral : public multicolvar::MultiColvar {
private:
//  double nl_cut;
  double rcut2;
  SwitchingFunction switchingFunction;
public:
  static void registerKeywords( Keywords& keys );
  Tetrahedral(const ActionOptions&);
// active methods:
  virtual double compute(); 
  Vector getCentralAtom();
/// Returns the number of coordinates of the field
  bool isPeriodic(){ return false; }
};

PLUMED_REGISTER_ACTION(Tetrahedral,"TETRAHEDRAL")

void Tetrahedral::registerKeywords( Keywords& keys ){
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

Tetrahedral::Tetrahedral(const ActionOptions&ao):
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
  setLinkCellCutoff( switchingFunction.get_dmax() );
  rcut2 = switchingFunction.get_dmax()*switchingFunction.get_dmax();

  // Read in the atoms
  int natoms=2; readAtoms( natoms );
  // And setup the ActionWithVessel
  checkRead();
}

double Tetrahedral::compute(){
   weightHasDerivatives=true;
   double value=0, norm=0, dfunc; Vector distance;

   // Calculate the coordination number
   Vector myder, fder;
   double sw, sp1, sp2, sp3, sp4;
   double sp1c, sp2c, sp3c, sp4c, r3, r5, tmp;
   double d2, t1, t2, t3, t4, tt1, tt2, tt3, tt4;
   for(unsigned i=1;i<getNAtoms();++i){
      distance=getSeparation( getPosition(0), getPosition(i) );
      d2 = distance.modulo2();
      if( d2<rcut2 ){ 
         sw = switchingFunction.calculateSqr( d2, dfunc );

         sp1 = +distance[0]+distance[1]+distance[2];
         sp2 = +distance[0]-distance[1]-distance[2];
         sp3 = -distance[0]+distance[1]-distance[2];
         sp4 = -distance[0]-distance[1]+distance[2];

         sp1c = pow( sp1, 3 );
         sp2c = pow( sp2, 3 );
         sp3c = pow( sp3, 3 );
         sp4c = pow( sp4, 3 );

         r3 = pow( distance.modulo(), 3 );
         r5 = pow( distance.modulo(), 5 );

         tmp = sp1c/r3 + sp2c/r3 + sp3c/r3 + sp4c/r3;

         t1=(3*sp1c)/r5; tt1=((3*sp1*sp1)/r3);  
         t2=(3*sp2c)/r5; tt2=((3*sp2*sp2)/r3);  
         t3=(3*sp3c)/r5; tt3=((3*sp3*sp3)/r3);  
         t4=(3*sp4c)/r5; tt4=((3*sp4*sp4)/r3);  

         myder[0] = (tt1-(distance[0]*t1))  + (tt2-(distance[0]*t2))  + (-tt3-(distance[0]*t3))  + (-tt4-(distance[0]*t4));
         myder[1] = (tt1-(distance[1]*t1))  + (-tt2-(distance[1]*t2))  + (tt3-(distance[1]*t3))  + (-tt4-(distance[1]*t4));
         myder[2] = (tt1-(distance[2]*t1))  + (-tt2-(distance[2]*t2))  + (-tt3-(distance[2]*t3))  + (tt4-(distance[2]*t4));

         value += sw*tmp; fder = (+dfunc)*tmp*distance + sw*myder;
         addAtomsDerivatives( 0, -fder );
         addAtomsDerivatives( i, +fder );
         // Tens is a constructor that you build by taking the vector product of two vectors (remember the scalars!)
         addBoxDerivatives( Tensor(distance,-fder) );
 
         norm += sw;
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

Vector Tetrahedral::getCentralAtom(){
   addCentralAtomDerivatives( 0, Tensor::identity() );
   return getPosition(0);
}

}
}

