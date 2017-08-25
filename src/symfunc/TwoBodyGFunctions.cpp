/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2017 The plumed team
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
#include "SymmetryFunctionBase.h"
#include "multicolvar/MultiColvarBase.h"
#include "core/ActionRegister.h"
#include <string>
#include <cmath>

using namespace std;

namespace PLMD {
namespace symfunc {


class TwoBodyGFunctions : public SymmetryFunctionBase {
private:
  double center2, nu2;
  double kappa3;
public:
  static void registerKeywords( Keywords& keys );
  explicit TwoBodyGFunctions(const ActionOptions&);
  void compute( const double& weight, const Vector& vec, MultiValue& myvals ) const ;
};

PLUMED_REGISTER_ACTION(TwoBodyGFunctions,"GSYMFUNC_TWOBODY")

void TwoBodyGFunctions::registerKeywords( Keywords& keys ) {
  SymmetryFunctionBase::registerKeywords( keys ); 
  keys.add("compulsory","CENTER2","position of the gaussian center in the G2 symmetry function");
  keys.add("compulsory","NU2","value of the width parameter for the gaussian in the G2 symmetry function");
  keys.add("compulsory","KAPPA3","value of kappa parameter in the G3 symmetry function");
  keys.addOutputComponent("g1","default","the value of the G1 symmetry function");
  keys.addOutputComponent("g2","default","the value of the G2 symmetry function");
  keys.addOutputComponent("g3","default","the value of the G3 symmetry function");
}

TwoBodyGFunctions::TwoBodyGFunctions(const ActionOptions&ao):
  Action(ao),
  SymmetryFunctionBase(ao)
{
  addComponentWithDerivatives( "g1" );
  parse("NU2",nu2); parse("CENTER2",center2); addComponentWithDerivatives( "g2" );
  log.printf("  for g2 function gaussian center is as %f and width paramter is %f \n", center2, nu2 );
  parse("KAPPA3",kappa3); addComponentWithDerivatives( "g3" );
  log.printf("  for g3 function value of kappa parameter is %f \n",kappa3 );
}

void TwoBodyGFunctions::compute( const double& val, const Vector& distance, MultiValue& myvals ) const {
  double dlen = distance.modulo();
  // Compute G1
  addToValue( 0, val, myvals ); addWeightDerivative( 0, 1.0, myvals );
  // Compute G2
  double diff = dlen - center2, ee = exp( - nu2*diff*diff );
  addToValue( 1, ee*val, myvals ); addWeightDerivative( 1, ee, myvals ); 
  addVectorDerivatives( 1, -(val/dlen)*2*nu2*diff*ee*distance, myvals );
  // Compute G3
  double cc = cos( kappa3*dlen ), ss = -kappa3*sin( kappa3*dlen );
  addToValue( 2, cc*val, myvals ); addWeightDerivative( 2, cc, myvals );
  addVectorDerivatives( 2, (val/dlen)*ss*distance, myvals );
}

}
}
