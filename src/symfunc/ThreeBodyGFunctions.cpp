/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2017 The plumed team
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
#include "core/ActionRegister.h"
#include "tools/SwitchingFunction.h"

namespace PLMD {
namespace symfunc {

class ThreeBodyGFunctions : public SymmetryFunctionBase {
private:
  SwitchingFunction sf4;
  double lambda4, nu4, zeta4, g4_prefactor;
  double lambda5, nu5, zeta5, g5_prefactor;
  double lambda6, zeta6, g6_prefactor;
  double nu7, alpha7;
public:
  static void registerKeywords( Keywords& keys );
  explicit ThreeBodyGFunctions(const ActionOptions&);
  void compute( const double& weight, const Vector& vec, MultiValue& myvals ) const { plumed_error(); }
  void computeSymmetryFunction( const unsigned& current, MultiValue& myvals ) const ;
};

PLUMED_REGISTER_ACTION(ThreeBodyGFunctions,"GSYMFUNC_THREEBODY")

void ThreeBodyGFunctions::registerKeywords( Keywords& keys ) {
  SymmetryFunctionBase::registerKeywords( keys ); keys.remove("ONESHOT");
  keys.add("compulsory","LAMBDA4","the lambda parameter in the G4 symmetry function");
  keys.add("compulsory","NU4","the NU parameteter in the G4 symmetry function");
  keys.add("compulsory","ZETA4","the zeta parameter in the G4 symmetry function");
  keys.add("compulsory","SWITCH4","the switching function that acts on the distance between atom j and atom k in the G4 symmetry function");
  keys.add("compulsory","LAMBDA5","the lambda parameter in the G5 symmetry function");
  keys.add("compulsory","NU5","the NU parameteter in the G5 symmetry function");
  keys.add("compulsory","ZETA5","the zeta parameter in the G5 symmetry function");
  keys.add("compulsory","LAMBDA6","the lambda parameter in the G6 symmetry function");
  keys.add("compulsory","ZETA6","the zeta parameter in the G6 symmetry function");
  keys.add("compulsory","NU7","the nu parameter in the G7 symmetry function");
  keys.add("compulsory","ALPHA7","the alpha parameter in the G7 symmetry function");
  keys.addOutputComponent("g4","default","the value of the G4 symmetry function");
  keys.addOutputComponent("g5","default","the value of the G5 symmetry function");
  keys.addOutputComponent("g6","default","the value of the G6 symmetry function");
  keys.addOutputComponent("g7","default","the value of the G7 symmetry function");
}

ThreeBodyGFunctions::ThreeBodyGFunctions(const ActionOptions&ao):
Action(ao),
SymmetryFunctionBase(ao)
{
   // Read input for G4
   std::string swstr, errors; parse("SWITCH4",swstr); sf4.set(swstr, errors );
   if( errors.length()!=0 ) error("problem reading switching function description " + errors );
   parse("LAMBDA4",lambda4); parse("NU4",nu4); parse("ZETA4",zeta4);
   addComponentWithDerivatives("g4"); g4_prefactor = 1.0 / pow(2.0,zeta4);
   log.printf("  parameters for g4 are lambda=%f, nu=%f and zeta=%f\n",lambda4, nu4, zeta4);
   log.printf("  using switching function distance between atom j and atom k with cutoff %s \n",sf4.description().c_str() );
   // Read input for G5
   parse("LAMBDA5",lambda5); parse("NU5",nu5); parse("ZETA5",zeta5);
   addComponentWithDerivatives("g5"); g5_prefactor = 1.0 / pow(2.0,zeta5);
   log.printf("  parameters for g5 are lambda=%f, nu=%f and zeta=%f\n",lambda5, nu5, zeta5); 
   // Read input for G6
   parse("LAMBDA6",lambda6); parse("ZETA6",zeta6);
   addComponentWithDerivatives("g6"); g6_prefactor = 1.0 / pow(2.0,zeta6);
   log.printf("  parameters for g6 are lambda=%f and zeta=%f\n",lambda6, zeta6);
   // Read input for G7
   parse("NU7",nu7); parse("ALPHA7",alpha7); 
   addComponentWithDerivatives("g7");
   log.printf("  parameters for g7 are nu=%f and alpha=%f\n",nu7, alpha7);
   checkRead();
}

void ThreeBodyGFunctions::computeSymmetryFunction( const unsigned& current, MultiValue& myvals ) const {
   Vector diri, dirj; double mod2i, mod2j, weighti, weightj;
   unsigned matind = getPntrToArgument(0)->getPositionInMatrixStash();
   unsigned matind_x = getPntrToArgument(1)->getPositionInMatrixStash();
   unsigned matind_y = getPntrToArgument(2)->getPositionInMatrixStash();
   unsigned matind_z = getPntrToArgument(3)->getPositionInMatrixStash();  
   for(unsigned i=1;i<myvals.getNumberOfStashedMatrixElements(matind);++i) {
       unsigned iind = myvals.getStashedMatrixIndex(matind,i);
       weighti = myvals.getStashedMatrixElement( matind, iind );
       if( weighti<epsilon ) continue ;
       diri[0] = myvals.getStashedMatrixElement( matind_x, iind );
       diri[1] = myvals.getStashedMatrixElement( matind_y, iind );
       diri[2] = myvals.getStashedMatrixElement( matind_z, iind );
       mod2i = diri.modulo2(); diri /= sqrt(mod2i);
       for(unsigned j=0;j<i;++j) {
           unsigned jind = myvals.getStashedMatrixIndex(matind,j);
           weightj = myvals.getStashedMatrixElement( matind, jind );
           if( weightj<epsilon ) continue ;
           dirj[0] = myvals.getStashedMatrixElement( matind_x, jind );
           dirj[1] = myvals.getStashedMatrixElement( matind_y, jind );
           dirj[2] = myvals.getStashedMatrixElement( matind_z, jind );
           mod2j = dirj.modulo2(); dirj /= sqrt(mod2j);
           Vector dirij = diri - dirj; double mod2ij = dirij.modulo2();
           double dder_ij, switchik = sf4.calculateSqr( mod2ij, dder_ij );
           // Compute cosine of angle 
           double cosa = dotProduct( diri, dirj );
           // Compute product of weights
           double weightij = weighti*weightj;
           // Compute G4
           double expg4 = exp( -nu4*( mod2i + mod2j + mod2ij ) );
           double g4p = pow( (1 + lambda4*cosa), zeta4 );
           addToValue( 0, g4_prefactor*g4p*expg4*weightij*switchik, myvals );
           // Compute G5
           double expg5 = exp( -nu4*( mod2i + mod2j ) );
           double g5p = pow( (1 + lambda5*cosa), zeta5 );
           addToValue( 1, g5_prefactor*g5p*expg5*weightij, myvals ); 
           // Compute G6
           double g6p = pow( (1 + lambda6*cosa), zeta6 );
           addToValue( 2, g6_prefactor*g6p*weightij, myvals );  
           // Compute G7
           double ang = acos( cosa ); double sina = sin( nu7*( ang - alpha7 ) );
           addToValue( 3, 0.5*sina*weightij, myvals );
       }
   }
}

}
}
