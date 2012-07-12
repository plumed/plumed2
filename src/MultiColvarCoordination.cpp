/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

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
#include "MultiColvar.h"
#include "NeighborList.h"
#include "ActionRegister.h"
#include "SwitchingFunction.h"

#include <string>
#include <cmath>
#include <cassert>

using namespace std;

namespace PLMD{

//+PLUMEDOC MCOLVAR COORDINATIONNUMBER
/*
Calculate the coordination numbers of atoms so that  you can then calculate functions of the distribution of
coordination numbers such as the minimum, the number less than a certain quantity and so on.   

\par Examples

The following input tells plumed to calculate the coordination numbers of atoms 1-100 with themselves.
The minimum coordination number is then calculated.
\verbatim
COORDINATIONNUMBER SPECIES=1-100 R_0=1.0 MIN=0.1
\endverbatim

The following input tells plumed to calculate how many atoms from 1-100 are within 3.0 of each of the atoms
from 101-110.  In the first 101 is the central atom, in the second 102 is the central atom and so on.  The 
number of coordination numbers more than 6 is then computed.
\verbatim
COORDINATIONNUMBER SPECIESA=101-110 SPECIESB=1-100 R_0=3.0 MORE_THAN=6.0
\endverbatim

\bug Analytical derivatives for off diagonal elements of virial are incorrect when you use GRADIENT/SUBCELL

*/
//+ENDPLUMEDOC


class MultiColvarCoordination : public MultiColvar {
private:
//  double nl_cut;
  SwitchingFunction switchingFunction;
public:
  static void registerKeywords( Keywords& keys );
  MultiColvarCoordination(const ActionOptions&);
// active methods:
  virtual double compute( const unsigned& j, const std::vector<Vector>& pos, std::vector<Vector>& deriv, Tensor& virial );
  void getCentralAtom( const std::vector<Vector>& pos, Vector& cpos, std::vector<Tensor>& deriv );
/// Returns the number of coordinates of the field
  unsigned getNumberOfFieldDerivatives();
  bool isPeriodic(const unsigned nn){ return false; }
};

PLUMED_REGISTER_ACTION(MultiColvarCoordination,"COORDINATIONNUMBER")

void MultiColvarCoordination::registerKeywords( Keywords& keys ){
  MultiColvar::registerKeywords( keys );
  ActionWithDistribution::autoParallelize( keys );
  keys.use("SPECIES"); keys.use("SPECIESA"); keys.use("SPECIESB");
  keys.add("optional","SWITCH","A switching function that defines which atoms are in the first coordination sphere. " + SwitchingFunction::documentation() );
  keys.add("optional","NN","The n parameter of the switching function ");
  keys.add("optional","MM","The m parameter of the switching function ");
  keys.add("optional","D_0","The d_0 parameter of the switching function");
  keys.add("optional","R_0","The r_0 parameter of the switching function");
  // Use actionWithDistributionKeywords
  keys.use("AVERAGE"); keys.use("MORE_THAN"); keys.use("LESS_THAN");
  keys.use("MIN"); keys.use("WITHIN"); keys.use("HISTOGRAM"); keys.use("MOMENTS");
  keys.use("SUBCELL"); keys.use("GRADIENT"); keys.use("DISTRIBUTION");
}

MultiColvarCoordination::MultiColvarCoordination(const ActionOptions&ao):
PLUMED_MULTICOLVAR_INIT(ao)
{
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
  log.printf("  coordination of central atom and those within %s\n",( switchingFunction.description() ).c_str() );

  // Read in the atoms
  int natoms; readAtoms( natoms );
  // And setup the ActionWithDistribution
  requestDistribution();

  // Create the groups for the neighbor list
  std::vector<AtomNumber> ga_lista, gb_lista; AtomNumber aa;
  aa.setIndex(0); ga_lista.push_back(aa);
  for(unsigned i=1;i<natoms;++i){ aa.setIndex(i); gb_lista.push_back(aa); }
  // And check everything has been read in correctly
  checkRead();
}

unsigned MultiColvarCoordination::getNumberOfFieldDerivatives(){
  return getNumberOfFunctionsInAction();
} 

double MultiColvarCoordination::compute( const unsigned& j, const std::vector<Vector>& pos, std::vector<Vector>& deriv, Tensor& virial ){
   double value=0, dfunc; Vector distance;

   // Calculate the coordination number
   double dd, sw;
   for(unsigned i=1;i<pos.size();++i){
      distance=getSeparation( pos[0], pos[i] ); dd=distance.modulo();
      sw = switchingFunction.calculate( distance.modulo(), dfunc );
      if( sw>=getTolerance() ){    //  nl_cut<0 ){
         value += sw;             // switchingFunction.calculate( distance.modulo(), dfunc );
         deriv[0] = deriv[0] + (-dfunc)*distance;
         deriv[i] = deriv[i] + (dfunc)*distance;
         virial = virial + (-dfunc)*Tensor(distance,distance);
      } else if( isTimeForNeighborListUpdate() ){
         removeAtomRequest( i );   
      }
   }

   return value;
}

void MultiColvarCoordination::getCentralAtom( const std::vector<Vector>& pos, Vector& cpos, std::vector<Tensor>& deriv ){
   cpos=pos[0]; deriv[0]=Tensor::identity();
}

}

