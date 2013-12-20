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

//+PLUMEDOC MCOLVARF LOCAL_AVERAGE
/*
Calculate averages over spherical regions centered on atoms

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace multicolvar {

class LocalAverage : public MultiColvarFunction {
private:
/// Ensures we deal with vectors properly
  unsigned jstart;
/// The values of the quantities we need to differentiate
  std::vector<double> values;
/// The switching function that tells us if atoms are close enough together
  SwitchingFunction switchingFunction;
public:
  static void registerKeywords( Keywords& keys );
  LocalAverage(const ActionOptions&);
/// We have to overwrite this here
  unsigned getNumberOfQuantities();
/// Actually do the calculation
  double compute();
/// This returns the position of the central atom
  Vector getCentralAtom();
/// Is the variable periodic
  bool isPeriodic(){ return false; }   
};

PLUMED_REGISTER_ACTION(LocalAverage,"LOCAL_AVERAGE")

void LocalAverage::registerKeywords( Keywords& keys ){
  MultiColvarFunction::registerKeywords( keys );
  keys.add("compulsory","NN","6","The n parameter of the switching function ");
  keys.add("compulsory","MM","12","The m parameter of the switching function ");
  keys.add("compulsory","D_0","0.0","The d_0 parameter of the switching function");
  keys.add("compulsory","R_0","The r_0 parameter of the switching function");
  keys.add("optional","SWITCH","This keyword is used if you want to employ an alternative to the continuous swiching function defined above. "
                               "The following provides information on the \\ref switchingfunction that are available. "
                               "When this keyword is present you no longer need the NN, MM, D_0 and R_0 keywords.");
  // Use actionWithDistributionKeywords
  keys.remove("LOWMEM"); keys.use("MEAN"); keys.use("MORE_THAN"); keys.use("LESS_THAN");
  keys.use("BETWEEN"); keys.use("HISTOGRAM"); keys.use("MOMENTS");
  keys.addFlag("LOWMEM",false,"lower the memory requirements");
}

LocalAverage::LocalAverage(const ActionOptions& ao):
Action(ao),
MultiColvarFunction(ao)
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
  log.printf("  averaging over central molecule and those within %s\n",( switchingFunction.description() ).c_str() );
  buildSymmetryFunctionLists( false );

  // One component for regular multicolvar and nelements for vectormulticolvar
  if( getBaseMultiColvar(0)->getNumberOfQuantities()==5 ){ values.resize( 1 ); jstart=0; } 
  else { values.resize( getBaseMultiColvar(0)->getNumberOfQuantities() - 5 ); jstart=5; }
}

unsigned LocalAverage::getNumberOfQuantities(){
  return jstart + values.size();
}

double LocalAverage::compute(){
  weightHasDerivatives=true;  

  Vector distance; double sw, dfunc, nbond=1;

  getValueForBaseTask( 0, values );
  for(unsigned j=0;j<values.size();++j) addElementValue( jstart + j, values[j] );

  accumulateWeightedAverageAndDerivatives( 0, 1.0 );
  for(unsigned i=1;i<getNAtoms();++i){
     distance=getSeparation( getPositionOfCentralAtom(0), getPositionOfCentralAtom(i) );
     sw = switchingFunction.calculateSqr( distance.modulo2(), dfunc );
     if( sw>=getTolerance() ){
         Tensor vir(distance,distance); 
         getValueForBaseTask( i, values ); 
         accumulateWeightedAverageAndDerivatives( i, sw );
         for(unsigned j=0;j<values.size();++j){
             addElementValue( jstart + j, sw*values[j] );
             addCentralAtomsDerivatives( 0, jstart+j, (-dfunc)*values[j]*distance );
             addCentralAtomsDerivatives( i, jstart+j, (+dfunc)*values[j]*distance );
             MultiColvarBase::addBoxDerivatives( jstart+j, (-dfunc)*values[j]*vir );   // This is a complex number?
         }
         nbond += sw;
         addCentralAtomsDerivatives( 0, 1, (-dfunc)*distance );
         addCentralAtomsDerivatives( i, 1, (+dfunc)*distance );
         MultiColvarBase::addBoxDerivatives( 1, (-dfunc)*vir );
     } else {
         removeAtomRequest( i, sw );
     }
  }

  // Set the tempory weight
  setElementValue( 1, nbond ); 
  // Update all dynamic lists
  updateActiveAtoms();
  // Finish the calculation
  getBaseMultiColvar(0)->finishWeightedAverageCalculation( this );

  // Clear working derivatives
  clearDerivativesAfterTask(1);
  // Weight doesn't really have derivatives (just use the holder for convenience)
  weightHasDerivatives=false; setElementValue( 1, 1.0 );   
   
  return getElementValue(0);
}

Vector LocalAverage::getCentralAtom(){
  addDerivativeOfCentralAtomPos( 0, Tensor::identity() ); 
  return getPositionOfCentralAtom(0); 
} 

}
}
