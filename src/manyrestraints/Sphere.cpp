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
#include "ManyRestraintsBase.h"
#include "core/ActionRegister.h"

namespace PLMD {
namespace manyrestraints {

class Sphere : public ManyRestraintsBase {
private:
  double at;
  double kappa;
  double exp;
  double eps;
  double offset;
  bool nopbc;
  Vector com;
  std::vector<double> com_deriv;
public:
  static void registerKeywords( Keywords& keys );
  Sphere( const ActionOptions& );
  bool isPeriodic(){ return false; }
  unsigned getNumberOfFunctionsInAction();
  void performTask( const unsigned& j );
  void calculate();
};

PLUMED_REGISTER_ACTION(Sphere,"SPHERICAL_RESTRAINT")

void Sphere::registerKeywords( Keywords& keys ){
  ManyRestraintsBase::registerKeywords( keys );
  keys.add("atoms","ATOMS","the atoms that are being confined to the sphere");
  keys.add("compulsory","RADIUS","the radius of the sphere");
  keys.add("compulsory","KAPPA","the force constant for the wall.  The k_i in the expression for a wall.");
  keys.add("compulsory","OFFSET","0.0","the offset for the start of the wall.  The o_i in the expression for a wall.");
  keys.add("compulsory","EXP","2.0","the powers for the walls.  The e_i in the expression for a wall.");
  keys.add("compulsory","EPS","1.0","the values for s_i in the expression for a wall");
  keys.addFlag("NOPBC",false,"turn off periodic boundary conditions");
}

Sphere::Sphere(const ActionOptions& ao):
Action(ao),
ManyRestraintsBase(ao)
{
  std::vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  com_deriv.resize( atoms.size() );

  parse("RADIUS",at);
  parse("OFFSET",offset);
  parse("EPS",eps);
  parse("EXP",exp);
  parse("KAPPA",kappa);
  parseFlag("NOPBC",nopbc);
  checkRead();

  requestAtoms( atoms ); 
  createRestraints( atoms.size() );
}

unsigned Sphere::getNumberOfFunctionsInAction(){
  return getNumberOfAtoms();
}

void Sphere::calculate(){
  // Calculate position of the center of mass
  double mass=0; com.zero();
  for(unsigned i=0;i<getNumberOfAtoms();i++) mass+=getMass(i);

  for(unsigned i=0;i<getNumberOfAtoms();i++){
    com+=(getMass(i)/mass)*getPosition(i);
    com_deriv[i]=(getMass(i)/mass);
  }
 
  // Now run the full set of tasks
  runAllTasks();
}

void Sphere::performTask( const unsigned& j ){
  Vector distance;

  if(!nopbc){
    distance=pbcDistance(com,getPosition(current));
  } else {
    distance=delta(com,getPosition(current));
  }

  double value=distance.modulo();
  double uscale = (value - at + offset)/eps;
  if( uscale > 0. ){
     double invvalue= 1.0 / value ;
     double power = pow( uscale, exp );
     double f = invvalue * ( kappa / eps ) * exp * power / uscale;

     setElementValue( 0, kappa*power ); setElementValue( 1, 1.0 );
     // Add derivatives for com
     for(unsigned i=0;i<getNumberOfAtoms();++i) addAtomsDerivatives( i, -com_deriv[i]*f*distance );

     // Add derivatives for other atom 
     addAtomsDerivatives( current, f*distance );

     // Add derivatives for virial
     addBoxDerivatives( -f*Tensor(distance,distance) );

     // We need to accumulate derivatives
     return;
  }

  // We need do nothing more if this is not true
  setElementValue( 1, 0.0 );
  return;
}

}
}

