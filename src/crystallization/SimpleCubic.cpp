/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2016 The plumed team
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
#include "multicolvar/MultiColvar.h"
#include "tools/NeighborList.h"
#include "core/ActionRegister.h"
#include "tools/SwitchingFunction.h"

#include <string>
#include <cmath>

using namespace std;

namespace PLMD{
namespace crystallization{

//+PLUMEDOC MCOLVAR SIMPLECUBIC
/*
Calculate whether or not the coordination spheres of atoms are arranged as they would be in a simple
cubic structure.

\par Examples

The following input tells plumed to calculate the simple cubic parameter for the atoms 1-100 with themselves.
The mean value is then calculated.
\verbatim
SIMPLECUBIC SPECIES=1-100 R_0=1.0 MEAN
\endverbatim

The following input tells plumed to look at the ways atoms 1-100 are within 3.0 are arranged about atoms
from 101-110.  The number of simple cubic parameters that are greater than 0.8 is then output
\verbatim
SIMPLECUBIC SPECIESA=101-110 SPECIESB=1-100 R_0=3.0 MORE_THAN={RATIONAL R_0=0.8 NN=6 MM=12 D_0=0}
\endverbatim

*/
//+ENDPLUMEDOC


class SimpleCubic : public multicolvar::MultiColvar {
private:
//  double nl_cut;
  double rcut2;
  SwitchingFunction switchingFunction;
public:
  static void registerKeywords( Keywords& keys );
  explicit SimpleCubic(const ActionOptions&);
// active methods:
  virtual double compute( const unsigned& tindex, multicolvar::AtomValuePack& myatoms ) const ; 
/// Returns the number of coordinates of the field
  bool isPeriodic(){ return false; }
};

PLUMED_REGISTER_ACTION(SimpleCubic,"SIMPLECUBIC")

void SimpleCubic::registerKeywords( Keywords& keys ){
  multicolvar::MultiColvar::registerKeywords( keys );
  keys.use("SPECIES"); keys.use("SPECIESA"); keys.use("SPECIESB");
  keys.add("compulsory","NN","6","The n parameter of the switching function ");
  keys.add("compulsory","MM","0","The m parameter of the switching function; 0 implies 2*NN");
  keys.add("compulsory","D_0","0.0","The d_0 parameter of the switching function");
  keys.add("compulsory","R_0","The r_0 parameter of the switching function");
  keys.add("optional","SWITCH","This keyword is used if you want to employ an alternative to the continuous swiching function defined above. "
                               "The following provides information on the \\ref switchingfunction that are available. "
                               "When this keyword is present you no longer need the NN, MM, D_0 and R_0 keywords.");
  // Use actionWithDistributionKeywords
  keys.use("MEAN"); keys.use("MORE_THAN"); keys.use("LESS_THAN"); keys.use("MAX");
  keys.use("MIN"); keys.use("BETWEEN"); keys.use("HISTOGRAM"); keys.use("MOMENTS");
  keys.use("ALT_MIN"); keys.use("LOWEST"); keys.use("HIGHEST"); 
}

SimpleCubic::SimpleCubic(const ActionOptions&ao):
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

double SimpleCubic::compute( const unsigned& tindex, multicolvar::AtomValuePack& myatoms ) const { 
   double d2, value=0, norm=0, dfunc; 

   // Calculate the coordination number
   Vector myder, fder;
   double sw, t1, t2, t3, x2, x3, x4, y2, y3, y4, z2, z3, z4, r4, tmp;
   for(unsigned i=1;i<myatoms.getNumberOfAtoms();++i){
      Vector& distance=myatoms.getPosition(i);  // getSeparation( myatoms.getPosition(0), myatoms.getPosition(i) );
      if ( (d2=distance[0]*distance[0])<rcut2 &&
           (d2+=distance[1]*distance[1])<rcut2 &&
           (d2+=distance[2]*distance[2])<rcut2) {

         sw = switchingFunction.calculateSqr( d2, dfunc );        

         x2 = distance[0]*distance[0];
         x3 = distance[0]*x2;
         x4 = distance[0]*x3;

         y2 = distance[1]*distance[1];
         y3 = distance[1]*y2;
         y4 = distance[1]*y3;         
   
         z2 = distance[2]*distance[2];
         z3 = distance[2]*z2; 
         z4 = distance[2]*z3;

         r4 = pow( distance.modulo2(), 2 );
         tmp = ( x4 + y4 + z4 ) / r4;

         t1=(x2+y2+z2); t2=t1*t1; t3=(x4+y4+z4)/(t1*t2);
         myder[0] = 4*x3/t2-4*distance[0]*t3; 
         myder[1] = 4*y3/t2-4*distance[1]*t3; 
         myder[2] = 4*z3/t2-4*distance[2]*t3; 

         value += sw*tmp; fder = (+dfunc)*tmp*distance + sw*myder;
         addAtomDerivatives( 1, 0, -fder, myatoms );
         addAtomDerivatives( 1, i, +fder, myatoms );
         // Tens is a constructor that you build by taking the vector product of two vectors (remember the scalars!)
         myatoms.addBoxDerivatives( 1, Tensor(distance,-fder) );
 
         norm += sw;
         addAtomDerivatives( -1, 0, (-dfunc)*distance, myatoms );
         addAtomDerivatives( -1, i, (+dfunc)*distance, myatoms );
         myatoms.addTemporyBoxDerivatives( (-dfunc)*Tensor(distance,distance) );
      }
   }
   
   myatoms.setValue(1, value);  
   // values -> der of... value [0], weight[1], x coord [2], y, z... [more magic]
   updateActiveAtoms( myatoms ); myatoms.getUnderlyingMultiValue().quotientRule( 1, norm, 1 ); 

   return value / norm; // this is equivalent to getting an "atomic" CV
}

}
}

