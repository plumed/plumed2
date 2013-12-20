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
#include "MultiColvarFunction.h"
#include "core/ActionRegister.h"
#include "tools/SwitchingFunction.h"

#include <string>
#include <cmath> 

//+PLUMEDOC MCOLVARF NLINKS 
/*
Calculate number of pairs of atoms/molecules that are "linked"

In its simplest guise this coordinate calculates a coordination number.  Each pair
of atoms is assumed "linked" if they are within some cutoff of each other.  In more 
complex applications each entity is a vector and this quantity measures whether
pairs of vectors are (a) within a certain cutoff and (b) if the two vectors have 
similar orientations

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace multicolvar {

class NumberOfLinks : public MultiColvarFunction {
private:
/// The values of the quantities in the dot products
  std::vector<double> orient0, orient1; 
/// The switching function that tells us if atoms are close enough together
  SwitchingFunction switchingFunction;
public:
  static void registerKeywords( Keywords& keys );
  NumberOfLinks(const ActionOptions&);
/// Do the stuff with the switching functions
  void calculateWeight();
/// Actually do the calculation
  double compute();
/// This returns the position of the central atom
  Vector getCentralAtom();  
/// Is the variable periodic
  bool isPeriodic(){ return false; }
};

PLUMED_REGISTER_ACTION(NumberOfLinks,"NLINKS")

void NumberOfLinks::registerKeywords( Keywords& keys ){
  MultiColvarFunction::registerKeywords( keys );
  keys.add("compulsory","NN","6","The n parameter of the switching function ");
  keys.add("compulsory","MM","12","The m parameter of the switching function ");
  keys.add("compulsory","D_0","0.0","The d_0 parameter of the switching function");
  keys.add("compulsory","R_0","The r_0 parameter of the switching function");
  keys.add("optional","SWITCH","This keyword is used if you want to employ an alternative to the continuous swiching function defined above. "
                               "The following provides information on the \\ref switchingfunction that are available. "
                               "When this keyword is present you no longer need the NN, MM, D_0 and R_0 keywords.");
  // Use actionWithDistributionKeywords
  keys.remove("LOWMEM"); 
  keys.addFlag("LOWMEM",false,"lower the memory requirements");
}

NumberOfLinks::NumberOfLinks(const ActionOptions& ao):
Action(ao),
MultiColvarFunction(ao)
{
  // The weight of this does have derivatives
  weightHasDerivatives=true;

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
  log.printf("  calculating number of links with atoms separation of %s\n",( switchingFunction.description() ).c_str() );
  buildAtomListWithPairs( false );

  for(unsigned i=0;i<getNumberOfBaseMultiColvars();++i){
    if( !getBaseMultiColvar(i)->hasDifferentiableOrientation() ) error("cannot use multicolvar of type " + getBaseMultiColvar(i)->getName() );
  }
  
  // Resize these ready for business
  if( getBaseMultiColvar(0)->getNumberOfQuantities()==5 ){ 
      orient0.resize( 1 ); orient1.resize( 1 ); 
  } else { 
      orient0.resize( getBaseMultiColvar(0)->getNumberOfQuantities() - 5 ); 
      orient1.resize( getBaseMultiColvar(0)->getNumberOfQuantities() - 5 );
  }

  // Create holders for the collective variable
  readVesselKeywords();
  plumed_assert( getNumberOfVessels()==0 );
  std::string input; addVessel( "SUM", input, -1 );
  readVesselKeywords();
}

void NumberOfLinks::calculateWeight(){
  Vector distance = getSeparation( getPositionOfCentralAtom(0), getPositionOfCentralAtom(1) );
  double dfunc, sw = switchingFunction.calculateSqr( distance.modulo2(), dfunc );
  addCentralAtomsDerivatives( 0, 1, (-dfunc)*distance );
  addCentralAtomsDerivatives( 1, 1, (dfunc)*distance );
  MultiColvarBase::addBoxDerivatives( 1, (-dfunc)*Tensor(distance,distance) );
  setElementValue(1,sw);
}

double NumberOfLinks::compute(){
   getValueForBaseTask( 0, orient0 ); 
   getValueForBaseTask( 1, orient1 );

   double dot=0;
   for(unsigned k=0;k<orient0.size();++k){
       dot+=orient0[k]*orient1[k];
   }
//   for(unsigned k=0;k<orient0.size();++k){
//      orient0[k]*=dot_df; orient1[k]*=dot_df;
//   }
   addOrientationDerivatives( 0, orient1 );
   addOrientationDerivatives( 1, orient0 );

   return dot;
}

Vector NumberOfLinks::getCentralAtom(){
   addDerivativeOfCentralAtomPos( 0, 0.5*Tensor::identity() );
   addDerivativeOfCentralAtomPos( 1, 0.5*Tensor::identity() );
   return 0.5*( getPositionOfCentralAtom(0) + getPositionOfCentralAtom(1) );
}

}
}
