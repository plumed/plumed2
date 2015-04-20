/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2015 The plumed team
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

As is explained in <a href="http://www.youtube.com/watch?v=iDvZmbWE5ps"> this video </a> certain multicolvars
calculate one scalar quantity or one vector for each of the atoms in the system.  For example 
\ref COORDINATIONNUMBER measures the coordination number of each of the atoms in the system and \ref Q4 measures
the 4th order Steinhardt parameter for each of the atoms in the system.  These quantities provide tell us something about
the disposition of the atoms in the first coordination sphere of each of the atoms of interest.  Lechner and Dellago \cite dellago-q6
have suggested that one can probe local order in a system by taking the average value of such symmetry functions over
the atoms within a spherical cutoff of each of these atoms in the systems.  When this is done with Steinhardt parameters
they claim this gives a coordinate that is better able to distinguish solid and liquid configurations of Lennard-Jones atoms. 

You can calculate such locally averaged quantities within plumed by using the LOCAL_AVERAGE command.  This command calculates 
the following atom-centered quantities:

\f[
s_i = \frac{ c_i + \sum_j \sigma(r_{ij})c_j }{ 1 + \sum_j \sigma(r_{ij}) } 
\f]

where the \f$c_i\f$ and \f$c_j\f$ values can be for any one of the symmetry functions that can be calculated using plumed 
multicolvars.  The function \f$\sigma( r_{ij} )\f$ is a \ref switchingfunction that acts on the distance between 
atoms \f$i\f$ and \f$j\f$.  Lechner and Dellago suggest that the parameters of this function should be set so that it the function is equal to one
when atom \f$j\f$ is in the first coordination sphere of atom \f$i\f$ and is zero otherwise.  

The \f$s_i\f$ quantities calculated using the above command can be again thought of as atom-centred symmetry functions.  They 
thus operate much like multicolvars.  You can thus calculate properties of the distribution of \f$s_i\f$ values using MEAN, LESS_THAN, HISTOGRAM
and so on.  You can also probe the value of these averaged variables in regions of the box by using the command in tandem with the 
\ref AROUND command.

\par Examples

This example input calculates the coordination numbers for all the atoms in the system.  These coordination numbers are then averaged over
spherical regions.  The number of averaged coordination numbers that are greater than 4 is then output to a file.

\verbatim
COORDINATIONNUMBER SPECIES=1-64 D_0=1.3 R_0=0.2 LABEL=d1
LOCAL_AVERAGE ARG=d1 SWITCH={RATIONAL D_0=1.3 R_0=0.2} MORE_THAN={RATIONAL R_0=4} LABEL=la
PRINT ARG=la.* FILE=colvar 
\endverbatim

This example input calculates the \f$q_4\f$ (see \ref Q4) vectors for each of the atoms in the system.  These vectors are then averaged 
component by component over a spherical region.  The average value for this quantity is then outputeed to a file.  This calculates the 
quantities that were used in the paper by Lechner and Dellago \cite dellago-q6 

\verbatim
Q4 SPECIES=1-64 SWITCH={RATIONAL D_0=1.3 R_0=0.2} LABEL=q4
LOCAL_AVERAGE ARG=q4 SWITCH={RATIONAL D_0=1.3 R_0=0.2} MEAN LABEL=la
PRINT ARG=la.* FILE=colvar
\endverbatim

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace multicolvar {

class LocalAverage : public MultiColvarFunction {
private:
/// Cutoff
  double rcut2;
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
  // One component for regular multicolvar and nelements for vectormulticolvar
  if( getBaseMultiColvar(0)->getNumberOfQuantities()==5 ){ values.resize( 1 ); jstart=0; }
  else { values.resize( getBaseMultiColvar(0)->getNumberOfQuantities() - 5 ); jstart=5; }

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
  rcut2 = switchingFunction.get_dmax()*switchingFunction.get_dmax();
  setLinkCellCutoff( switchingFunction.get_dmax() ); buildSymmetryFunctionLists();
  for(unsigned i=0;i<getNumberOfBaseMultiColvars();++i) getBaseMultiColvar(i)->doNotCalculateDirector();
}

unsigned LocalAverage::getNumberOfQuantities(){
  return jstart + values.size();
}

double LocalAverage::compute(){
  weightHasDerivatives=true;  

  Vector distance; double d2, sw, dfunc, nbond=1;

  getVectorForBaseTask( 0, values ); 
  for(unsigned j=0;j<values.size();++j) addElementValue( jstart + j, values[j] );

  accumulateWeightedAverageAndDerivatives( 0, 1.0 );
  for(unsigned i=1;i<getNAtoms();++i){
     distance=getSeparation( getPositionOfCentralAtom(0), getPositionOfCentralAtom(i) );
     d2 = distance.modulo2(); 
     if( d2<rcut2 ){
         sw = switchingFunction.calculateSqr( d2, dfunc );
         Tensor vir(distance,distance); 
         getVectorForBaseTask( i, values ); 
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
