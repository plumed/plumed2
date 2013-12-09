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
#include "MultiColvar.h"
#include "tools/Angle.h"
#include "tools/SwitchingFunction.h"
#include "core/ActionRegister.h"

#include <string>
#include <cmath>

using namespace std;

namespace PLMD{
namespace multicolvar{

//+PLUMEDOC MCOLVAR ANGLES 
/*
Calculate functions of the distribution of angles .

You can use this command to calculate functions such as:

\f[
 f(x) = \sum_{ijk} g( \theta_{ijk} )
\f] 

Alternatively you can use this command to calculate functions such as:

\f[
f(x) = \sum_{ijk} s(r_{ij})s(r_{jk}) g(\theta_{ijk})
\f]

where \f$s(r)\f$ is a \ref switchingfunction.  This second form means that you can
use this to calculate functions of the angles in the first coordination sphere of 
an atom / molecule \cite lj-recon.

\par Examples

The following example instructs plumed to find the average of two angles and to
print it to a file

\verbatim
ANGLES ATOMS=1,2,3 ATOMS=4,5,6 MEAN LABEL=a1
PRINT ARG=a1.mean FILE=colvar 
\endverbatim

The following example tells plumed to calculate all angles involving
at least one atom from GROUPA and two atoms from GROUPB in which the distances
are less than 1.0. The number of angles between \f$\frac{\pi}{4}\f$ and
\f$\frac{3\pi}{4}\f$ is then output  

\verbatim
ANGLES GROUPA=1-10 GROUPB=11-100 BETWEEN={GAUSSIAN LOWER=0.25pi UPPER=0.75pi} SWITCH={GAUSSIAN R_0=1.0} LABEL=a1
PRINT ARG=a1.between FILE=colvar
\endverbatim 

This final example instructs plumed to calculate all the angles in the first coordination
spheres of the atoms. A discretized-normalized histogram of the distribution is then output

\verbatim
ANGLE GROUP=1-38 HISTOGRAM={GAUSSIAN LOWER=0.0 UPPER=pi NBINS=20} SWTICH={GAUSSIAN R_0=1.0} LABEL=a1
PRINT ARG=a1.* FILE=colvar
\endverbatim

*/
//+ENDPLUMEDOC

class Angles : public MultiColvar {
private:
  bool use_sf;
  Vector dij, dik;
  SwitchingFunction sf1;
  SwitchingFunction sf2;
public:
  static void registerKeywords( Keywords& keys );
  Angles(const ActionOptions&);
/// Updates neighbor list
  virtual void doJobsRequiredBeforeTaskList();
  virtual double compute();
/// Returns the number of coordinates of the field
  void calculateWeight();
  bool isPeriodic(){ return false; }
  Vector getCentralAtom();
};

PLUMED_REGISTER_ACTION(Angles,"ANGLES")

void Angles::registerKeywords( Keywords& keys ){
  MultiColvar::registerKeywords( keys );
  keys.use("ATOMS"); keys.use("MEAN"); keys.use("LESS_THAN"); 
  keys.use("BETWEEN"); keys.use("HISTOGRAM"); keys.use("MORE_THAN");
  // Could also add Region here in theory
  keys.add("atoms-1","GROUP","Calculate angles for each distinct set of three atoms in the group");
  keys.add("atoms-2","GROUPA","A group of central atoms about which angles should be calculated");
  keys.add("atoms-2","GROUPB","When used in conjuction with GROUPA this keyword instructs plumed "
                              "to calculate all distinct angles involving one atom from GROUPA "
                              "and two atoms from GROUPB. The atom from GROUPA is the central atom.");
  keys.add("atoms-3","GROUPC","This must be used in conjuction with GROUPA and GROUPB.  All angles " 
                              "involving one atom from GROUPA, one atom from GROUPB and one atom from "
                              "GROUPC are calculated. The GROUPA atoms are assumed to be the central "
                              "atoms");
  keys.add("optional","SWITCH","A switching function that ensures that only angles between atoms that "
                               "are within a certain fixed cutoff are calculated. The following provides "
                               "information on the \\ref switchingfunction that are available.");
  keys.add("optional","SWITCHA","A switching function on the distance between the atoms in group A and the atoms in "
                                "group B");
  keys.add("optional","SWITCHB","A switching function on the distance between the atoms in group A and the atoms in "
                                "group B");
}

Angles::Angles(const ActionOptions&ao):
PLUMED_MULTICOLVAR_INIT(ao),
use_sf(false)
{
  std::string sfinput,errors; parse("SWITCH",sfinput);
  if( sfinput.length()>0 ){
      use_sf=true;
      weightHasDerivatives=true;
      sf1.set(sfinput,errors);
      if( errors.length()!=0 ) error("problem reading SWITCH keyword : " + errors ); 
      sf2.set(sfinput,errors);
      if( errors.length()!=0 ) error("problem reading SWITCH keyword : " + errors ); 
      log.printf("  only calculating angles for atoms separated by less than %s\n", sf1.description().c_str() );
  } else {
      parse("SWITCHA",sfinput); 
      if(sfinput.length()>0){
         use_sf=true;
         weightHasDerivatives=true;
         sf1.set(sfinput,errors);
         if( errors.length()!=0 ) error("problem reading SWITCHA keyword : " + errors ); 
         sfinput.clear(); parse("SWITCHB",sfinput);
         if(sfinput.length()==0) error("found SWITCHA keyword without SWITCHB");
         sf2.set(sfinput,errors); 
         if( errors.length()!=0 ) error("problem reading SWITCHB keyword : " + errors );
         log.printf("  only calculating angles when the distance between GROUPA and GROUPB atoms is less than %s\n", sf1.description().c_str() );
         log.printf("  only calculating angles when the distance between GROUPA and GROUPC atoms is less than %s\n", sf2.description().c_str() );
      }
  }
  // Read in the atoms
  int natoms=3; readAtoms( natoms );
  // And check everything has been read in correctly
  checkRead();
}

// This should give big speed ups during neighbor list update steps
void Angles::doJobsRequiredBeforeTaskList(){
  // Do jobs required by action with vessel
  ActionWithVessel::doJobsRequiredBeforeTaskList();
  if( !use_sf || getCurrentNumberOfActiveTasks()==ablocks[0].size() ) return ;
  // First step of update of three body neighbor list
  threeBodyNeighborList( sf1 );
} 

void Angles::calculateWeight(){
  dij=getSeparation( getPosition(0), getPosition(2) );
  dik=getSeparation( getPosition(0), getPosition(1) );
  if(!use_sf){ setWeight(1.0); return; }

  double w1, w2, dw1, dw2, wtot;
  w1=sf1.calculateSqr( dij.modulo2(), dw1 );
  w2=sf2.calculateSqr( dik.modulo2(), dw2 );
  wtot=w1*w2; dw1*=w2; dw2*=w1; 

  setWeight( wtot );
  if( wtot<getTolerance() ) return; 
  addAtomsDerivativeOfWeight( 1, dw2*dik );
  addAtomsDerivativeOfWeight( 0, -dw1*dij - dw2*dik ); 
  addAtomsDerivativeOfWeight( 2, dw1*dij );
  addBoxDerivativesOfWeight( (-dw1)*Tensor(dij,dij) + (-dw2)*Tensor(dik,dik) );
}

double Angles::compute(){
  Vector ddij,ddik; PLMD::Angle a; 
  double angle=a.compute(dij,dik,ddij,ddik);

  // And finish the calculation
  addAtomsDerivatives( 1, ddik );
  addAtomsDerivatives( 0, - ddik - ddij );
  addAtomsDerivatives( 2, ddij );
  addBoxDerivatives( -(Tensor(dij,ddij)+Tensor(dik,ddik)) );

  return angle;
}

Vector Angles::getCentralAtom(){
   addCentralAtomDerivatives( 0, Tensor::identity() );
   return getPosition(0);
}

}
}
