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
#include "ManyRestraintsBase.h"
#include "core/ActionRegister.h"

namespace PLMD {
namespace manyrestraints {

class UWalls : public ManyRestraintsBase {
private:
  double at;
  double kappa;
  double exp;
  double eps;
  double offset;
public:
  static void registerKeywords( Keywords& keys );
  UWalls( const ActionOptions& );
  void performTask();
};

PLUMED_REGISTER_ACTION(UWalls,"UWALLS")

void UWalls::registerKeywords( Keywords& keys ){
  ManyRestraintsBase::registerKeywords( keys );
  keys.add("compulsory","AT","the radius of the sphere");
  keys.add("compulsory","KAPPA","the force constant for the wall.  The k_i in the expression for a wall.");
  keys.add("compulsory","OFFSET","0.0","the offset for the start of the wall.  The o_i in the expression for a wall.");
  keys.add("compulsory","EXP","2.0","the powers for the walls.  The e_i in the expression for a wall.");
  keys.add("compulsory","EPS","1.0","the values for s_i in the expression for a wall");
}

UWalls::UWalls(const ActionOptions& ao):
Action(ao),
ManyRestraintsBase(ao)
{
  parse("AT",at);
  parse("OFFSET",offset);
  parse("EPS",eps);
  parse("EXP",exp);
  parse("KAPPA",kappa);
  checkRead();
}

void UWalls::performTask(){
  double value=getValue(); 
  double uscale = (value - at + offset)/eps;
  if( uscale > 0. ){
     double invvalue= 1.0 / value;
     double power = pow( uscale, exp );
     double f = ( kappa / eps ) * exp * power / uscale;

     setElementValue( 0, kappa*power ); setElementValue( 1, getWeight() );
     // Add derivatives 
     applyChainRuleForDerivatives( f );
    
     return;
  }

  // We need do nothing more if this is not true
  setElementValue( 1, 0.0 );
  return;
}

}
}

