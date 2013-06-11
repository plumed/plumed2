/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013 The plumed team
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
#include "tools/Torsion.h"
#include "core/ActionRegister.h"

#include <string>
#include <cmath>

using namespace std;

namespace PLMD{
namespace multicolvar{

//+PLUMEDOC MCOLVAR ALPHABETA 
/*

*/
//+ENDPLUMEDOC

class AlphaBeta : public MultiColvar {
private:
  double target;
public:
  static void registerKeywords( Keywords& keys );
  AlphaBeta(const ActionOptions&);
  virtual double compute( const unsigned& j );
  bool isPeriodic(){ return false; }
  Vector getCentralAtom();  
};

PLUMED_REGISTER_ACTION(AlphaBeta,"ALPHABETA")

void AlphaBeta::registerKeywords( Keywords& keys ){
  MultiColvar::registerKeywords( keys );
  keys.use("ATOMS");
  keys.add("optional","REFERENCE","a single reference value for all the dihedrals");
  //keys.add("optional","REFERENCE1","specific reference value for each dihedral");
}

AlphaBeta::AlphaBeta(const ActionOptions&ao):
PLUMED_MULTICOLVAR_INIT(ao)
{
  // Read in the atoms
  int natoms=4; readAtoms( natoms );
  parse("REFERENCE",target);
  // And setup the ActionWithVessel
  std::string fake_input;
  addVessel( "SUM", fake_input, -1 );  // -1 here means that this value will be named getLabel()
  readVesselKeywords();
  // And check everything has been read in correctly
  checkRead();
}

double AlphaBeta::compute( const unsigned& j ){
  Vector d0,d1,d2;
  d0=getSeparation(getPosition(1),getPosition(0));
  d1=getSeparation(getPosition(2),getPosition(1));
  d2=getSeparation(getPosition(3),getPosition(2));

  Vector dd0,dd1,dd2;
  PLMD::Torsion t;
  double value  = t.compute(d0,d1,d2,dd0,dd1,dd2);
  double svalue = -0.5*sin(value-target);
  double cvalue = 1.+cos(value-target);

  dd0 *= -svalue;
  dd1 *= -svalue;
  dd2 *= -svalue;
  value = 0.5*cvalue;

  addAtomsDerivatives(0,dd0);
  addAtomsDerivatives(1,dd1-dd0);
  addAtomsDerivatives(2,dd2-dd1);
  addAtomsDerivatives(3,-dd2);

  addBoxDerivatives  (-(extProduct(d0,dd0)+extProduct(d1,dd1)+extProduct(d2,dd2)));

  return value;
}

Vector AlphaBeta::getCentralAtom(){
   addCentralAtomDerivatives( 0, 0.25*Tensor::identity() );
   addCentralAtomDerivatives( 1, 0.25*Tensor::identity() );
   addCentralAtomDerivatives( 2, 0.25*Tensor::identity() );
   addCentralAtomDerivatives( 3, 0.25*Tensor::identity() );
   return 0.25*( getPosition(0) + getPosition(1) + getPosition(2) + getPosition(3) );
}

}
}
