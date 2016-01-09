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

//+PLUMEDOC MCOLVAR TETRAHEDRAL
/*
Calculate the degree to which the environment about ions has a tetrahedral order.

We can measure the degree to which the first coordination shell around any atom, \f$i\f$ is 
tetrahedrally ordered using the following function.

\f[
 s(i) = \frac{1}{\sum_j \sigma( r_{ij} )} \sum_j \sigma( r_{ij} )\left[ \frac{(x_{ij} + y_{ij} + z_{ij})^3}{r_{ij}^3} + 
                                                                        \frac{(x_{ij} - y_{ij} - z_{ij})^3}{r_{ij}^3} + 
                                                                        \frac{(-x_{ij} + y_{ij} - z_{ij})^3}{r_{ij}^3} + 
                                                                        \frac{(-x_{ij} - y_{ij} + z_{ij})^3}{r_{ij}^3} \right]  
\f]

Here \f$r_{ij}\f$ is the magnitude fo the vector connecting atom \f$i\f$ to atom \f$j\f$ and \f$x_{ij}\f$, \f$y_{ij}\f$ and \f$z_{ij}\f$
are its three components.  The function  \f$\sigma( r_{ij} )\f$ is a \ref switchingfunction that acts on the distance between 
atoms \f$i\f$ and \f$j\f$.  The parameters of this function should be set so that the function is equal to one
when atom \f$j\f$ is in the first coordination sphere of atom \f$i\f$ and is zero otherwise.

\par Examples

The following command calculates the average value of the tetrahedrality parameter for a set of 64 atoms all of the same type
and outputs this quantity to a file called colvar.

\verbatim
tt: TETRAHEDRAL SPECIES=1-64 SWITCH={RATIONAL D_0=1.3 R_0=0.2} MEAN
PRINT ARG=tt.mean FILE=colvar
\endverbatim  

The following command calculates the number of tetrahedrality parameters that are greater than 0.8 in a set of 10 atoms.
In this calculation it is assumed that there are two atom types A and B and that the first coordination sphere of the
10 atoms of type A contains atoms of type B.  The formula above is thus calculated for ten different A atoms and within 
it the sum over \f$j\f$ runs over 40 atoms of type B that could be in the first coordination sphere.

\verbatim
tt: TETRAHEDRAL SPECIESA=1-10 SPECIESB=11-40 SWITCH={RATIONAL D_0=1.3 R_0=0.2} MORE_THAN={RATIONAL R_0=0.8}
PRINT ARG=tt.* FILE=colvar
\endverbatim

*/
//+ENDPLUMEDOC


class Tetrahedral : public multicolvar::MultiColvar {
private:
//  double nl_cut;
  double rcut2;
  SwitchingFunction switchingFunction;
public:
  static void registerKeywords( Keywords& keys );
  explicit Tetrahedral(const ActionOptions&);
// active methods:
  virtual double compute( const unsigned& tindex, multicolvar::AtomValuePack& myatoms ) const ; 
/// Returns the number of coordinates of the field
  bool isPeriodic(){ return false; }
};

PLUMED_REGISTER_ACTION(Tetrahedral,"TETRAHEDRAL")

void Tetrahedral::registerKeywords( Keywords& keys ){
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

double Tetrahedral::compute( const unsigned& tindex, multicolvar::AtomValuePack& myatoms ) const {
   double value=0, norm=0, dfunc; 

   // Calculate the coordination number
   Vector myder, fder;
   double sw, sp1, sp2, sp3, sp4;
   double sp1c, sp2c, sp3c, sp4c, r3, r5, tmp;
   double d2, t1, t2, t3, t4, tt1, tt2, tt3, tt4;
   for(unsigned i=1;i<myatoms.getNumberOfAtoms();++i){
      Vector& distance=myatoms.getPosition(i);  // getSeparation( myatoms.getPosition(0), myatoms.getPosition(i) );
      if ( (d2=distance[0]*distance[0])<rcut2 &&
           (d2+=distance[1]*distance[1])<rcut2 &&
           (d2+=distance[2]*distance[2])<rcut2) {
      
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

