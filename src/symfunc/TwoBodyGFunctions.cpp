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
  int g1_ind, g2_ind, g3_ind;
public:
  static void registerKeywords( Keywords& keys );
  explicit TwoBodyGFunctions(const ActionOptions&);
  void compute( const double& weight, const Vector& vec, MultiValue& myvals ) const ;
};

PLUMED_REGISTER_ACTION(TwoBodyGFunctions,"GSYMFUNC_TWOBODY")

void TwoBodyGFunctions::registerKeywords( Keywords& keys ) {
  SymmetryFunctionBase::registerKeywords( keys ); 
  keys.addFlag("NO_G1",false,"do not compute the G1 symmetry function");
  keys.add("compulsory","CENTER2","position of the gaussian center in the G2 symmetry function");
  keys.add("optional","NU2","value of the width parameter for the gaussian in the G2 symmetry function");
  keys.add("optional","KAPPA3","value of kappa parameter in the G3 symmetry function");
  keys.addOutputComponent("g1","NO_G1","the value of the G1 symmetry function - this is computed by default unless the flag NO_G1 is given");
  keys.addOutputComponent("g2","CENTER2","the value of the G2 symmetry function");
  keys.addOutputComponent("g3","KAPPA3","the value of the G3 symmetry function");
}

TwoBodyGFunctions::TwoBodyGFunctions(const ActionOptions&ao):
  Action(ao),
  SymmetryFunctionBase(ao),
  g1_ind(-1),g2_ind(-1),g3_ind(-1)
{
  bool nog1; parseFlag("NO_G1",nog1 ); unsigned oflag=0;
  if( !nog1 ){ addComponentWithDerivatives( "g1" ); g1_ind=oflag; oflag++; }
  nu2=0.0; parse("NU2",nu2);
  if( nu2!=0.0 ){
     parse("CENTER2",center2); 
     addComponentWithDerivatives( "g2" );
     g2_ind=oflag; oflag++;
     log.printf("  for g2 function gaussian center is as %f and width paramter is %f \n", center2, nu2 );  
  }
  kappa3=0.0; parse("KAPPA3",kappa3); 
  if( kappa3!=0.0 ) {
      addComponentWithDerivatives( "g3" );
      g3_ind=oflag; oflag++;
      log.printf("  for g3 function value of kappa parameter is %f \n",kappa3 );
  }
}

void TwoBodyGFunctions::compute( const double& val, const Vector& distance, MultiValue& myvals ) const {
  double dlen = distance.modulo();
  // Compute G1
  if( g1_ind>-1 ) { addToValue( g1_ind, val, myvals ); addWeightDerivative( g1_ind, 1.0, myvals ); }
  // Compute G2
  if( g2_ind>-1 ) {
     double diff = dlen - center2, ee = exp( - nu2*diff*diff );
     addToValue( g2_ind, ee*val, myvals ); addWeightDerivative( g2_ind, ee, myvals ); 
     addVectorDerivatives( g2_ind, -(val/dlen)*2*nu2*diff*ee*distance, myvals );
  }
  // Compute G3
  if( g3_ind>-1 ) {
      double cc = cos( kappa3*dlen ), ss = -kappa3*sin( kappa3*dlen );
      addToValue( g3_ind, cc*val, myvals ); addWeightDerivative( g3_ind, cc, myvals );
      addVectorDerivatives( g3_ind, (val/dlen)*ss*distance, myvals );
  }
}

}
}
