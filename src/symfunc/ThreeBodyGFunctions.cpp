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
#include "tools/Angle.h"

namespace PLMD {
namespace symfunc {

class ThreeBodyGFunctions : public SymmetryFunctionBase {
private:
  SwitchingFunction sf4;
  int g4_ind, g5_ind, g6_ind, g7_ind;
  double lambda4, nu4, zeta4, g4_prefactor;
  double lambda5, nu5, zeta5, g5_prefactor;
  double lambda6, zeta6, g6_prefactor;
  double nu7, alpha7, g7_prefactor;
public:
  static void registerKeywords( Keywords& keys );
  explicit ThreeBodyGFunctions(const ActionOptions&);
  void compute( const double& weight, const Vector& vec, MultiValue& myvals ) const { plumed_error(); }
  void computeSymmetryFunction( const unsigned& current, MultiValue& myvals ) const ;
};

PLUMED_REGISTER_ACTION(ThreeBodyGFunctions,"GSYMFUNC_THREEBODY")

void ThreeBodyGFunctions::registerKeywords( Keywords& keys ) {
  SymmetryFunctionBase::registerKeywords( keys ); keys.remove("ONESHOT");
  keys.add("optional","LAMBDA4","the lambda parameter in the G4 symmetry function");
  keys.add("compulsory","NU4","the NU parameteter in the G4 symmetry function");
  keys.add("compulsory","ZETA4","the zeta parameter in the G4 symmetry function");
  keys.add("compulsory","SWITCH4","the switching function that acts on the distance between atom j and atom k in the G4 symmetry function");
  keys.add("optional","LAMBDA5","the lambda parameter in the G5 symmetry function");
  keys.add("compulsory","NU5","the NU parameteter in the G5 symmetry function");
  keys.add("compulsory","ZETA5","the zeta parameter in the G5 symmetry function");
  keys.add("optional","LAMBDA6","the lambda parameter in the G6 symmetry function");
  keys.add("compulsory","ZETA6","the zeta parameter in the G6 symmetry function");
  keys.add("optional","NU7","the nu parameter in the G7 symmetry function");
  keys.add("compulsory","ALPHA7","the alpha parameter in the G7 symmetry function");
  keys.addOutputComponent("g4","LAMBDA4","the value of the G4 symmetry function");
  keys.addOutputComponent("g5","LAMBDA5","the value of the G5 symmetry function");
  keys.addOutputComponent("g6","LAMBDA6","the value of the G6 symmetry function");
  keys.addOutputComponent("g7","NU7","the value of the G7 symmetry function");
}

ThreeBodyGFunctions::ThreeBodyGFunctions(const ActionOptions&ao):
Action(ao),
SymmetryFunctionBase(ao),
g4_ind(-1), g5_ind(-1), g6_ind(-1), g7_ind(-1)
{
   unsigned vindex=0;
   // Read input for G4
   lambda4=0.0; parse("LAMBDA4",lambda4); 
   if( lambda4!=0 ) {
       parse("NU4",nu4); parse("ZETA4",zeta4);
       addComponentWithDerivatives("g4"); g4_prefactor = 2.0 / pow(2.0,zeta4);
       log.printf("  parameters for g4 are lambda=%f, nu=%f and zeta=%f\n",lambda4, nu4, zeta4);
       std::string swstr, errors; parse("SWITCH4",swstr); sf4.set(swstr, errors );
       if( errors.length()!=0 ) error("problem reading switching function description " + errors ); 
       log.printf("  using switching function distance between atom j and atom k with cutoff %s \n",sf4.description().c_str() );
       g4_ind=vindex; vindex++;
   }
   // Read input for G5
   lambda5=0.0; parse("LAMBDA5",lambda5); 
   if( lambda5!=0 ){
       parse("NU5",nu5); parse("ZETA5",zeta5);
       addComponentWithDerivatives("g5"); g5_prefactor = 2.0 / pow(2.0,zeta5);
       log.printf("  parameters for g5 are lambda=%f, nu=%f and zeta=%f\n",lambda5, nu5, zeta5); 
       g5_ind=vindex; vindex++;
   }
   // Read input for G6
   lambda6=0.0; parse("LAMBDA6",lambda6); 
   if( lambda6!=0 ) {
       parse("ZETA6",zeta6);
       addComponentWithDerivatives("g6"); g6_prefactor = 2.0 / pow(2.0,zeta6);
       log.printf("  parameters for g6 are lambda=%f and zeta=%f\n",lambda6, zeta6);
       g6_ind=vindex; vindex++;
   }
   // Read input for G7
   nu7=0.0; parse("NU7",nu7); 
   if( nu7!=0 ) {
       parse("ALPHA7",alpha7); 
       addComponentWithDerivatives("g7"); g7_prefactor = 2.0;
       log.printf("  parameters for g7 are nu=%f and alpha=%f\n",nu7, alpha7);
       g7_ind=vindex; vindex++; 
   }
   checkRead();
}

void ThreeBodyGFunctions::computeSymmetryFunction( const unsigned& current, MultiValue& myvals ) const {
   Vector dd_i, dd_j, disti, distj;  double mod2i, mod2j, weighti, weightj; Angle angle;
   unsigned matind = getPntrToArgument(0)->getPositionInMatrixStash();
   unsigned matind_x = getPntrToArgument(1)->getPositionInMatrixStash();
   unsigned matind_y = getPntrToArgument(2)->getPositionInMatrixStash();
   unsigned matind_z = getPntrToArgument(3)->getPositionInMatrixStash();  
   for(unsigned i=1;i<myvals.getNumberOfStashedMatrixElements(matind);++i) {
       unsigned iind = myvals.getStashedMatrixIndex(matind,i);
       weighti = myvals.getStashedMatrixElement( matind, iind );
       if( weighti<epsilon ) continue ;
       disti[0] = myvals.getStashedMatrixElement( matind_x, iind );
       disti[1] = myvals.getStashedMatrixElement( matind_y, iind );
       disti[2] = myvals.getStashedMatrixElement( matind_z, iind );
       mod2i = disti.modulo2(); 
       for(unsigned j=0;j<i;++j) {
           unsigned jind = myvals.getStashedMatrixIndex(matind,j);
           weightj = myvals.getStashedMatrixElement( matind, jind );
           if( weightj<epsilon ) continue ;
           distj[0] = myvals.getStashedMatrixElement( matind_x, jind );
           distj[1] = myvals.getStashedMatrixElement( matind_y, jind );
           distj[2] = myvals.getStashedMatrixElement( matind_z, jind );
           mod2j = distj.modulo2(); 
           Vector dirij = disti - distj; double mod2ij = dirij.modulo2();
           // Compute angle betwee bonds
           double ang = angle.compute( disti, distj, dd_i, dd_j );
           // And cosine and sine of angle
           double cosa = cos(ang), sina = sin(ang); 
           // Compute product of weights
           double weightij = weighti*weightj;
           // Compute G4
           if( g4_ind>-1 ){
               double dder_ij, switchij = sf4.calculateSqr( mod2ij, dder_ij );
               double expg4 = exp( -nu4*( mod2i + mod2j + mod2ij ) );
               double g4v = (1 + lambda4*cosa);
               double g4p = pow( g4v, zeta4-1 );
               double g4_nonweight = g4_prefactor*g4p*g4v*expg4*switchij;
               double g4_vder1 = -g4_prefactor*weightij*zeta4*g4p*lambda4*sina*expg4*switchij; 
               double g4_mder = g4_prefactor*weightij*g4p*g4v*expg4;
               double g4_vder2 = -g4_mder*2*nu4*switchij; 
               double g4_vder3 = g4_mder*(dder_ij -2*nu4*switchij);
               addToValue( g4_ind, g4_nonweight*weightij, myvals );
               myvals.setSymfuncTemporyIndex( iind ); 
               addWeightDerivative( g4_ind, g4_nonweight*weightj, myvals );
               addVectorDerivatives(g4_ind, g4_vder1*dd_i, myvals ); 
               addVectorDerivatives(g4_ind, g4_vder2*disti, myvals );
               addVectorDerivatives(g4_ind, g4_vder3*dirij, myvals );
               myvals.setSymfuncTemporyIndex( jind );
               addWeightDerivative( g4_ind, g4_nonweight*weighti, myvals );
               addVectorDerivatives(g4_ind, g4_vder1*dd_j, myvals );
               addVectorDerivatives(g4_ind, g4_vder2*distj, myvals );
               addVectorDerivatives(g4_ind, -g4_vder3*dirij, myvals );
           }
           // Compute G5
           if( g5_ind>-1 ) {
               double expg5 = exp( -nu5*( mod2i + mod2j ) );
               double g5v = (1 + lambda5*cosa);
               double g5p = pow( g5v, zeta5-1 ); 
               double g5_nonweight = g5_prefactor*g5p*g5v*expg5;
               double g5_vder1 = -g5_prefactor*weightij*expg5*zeta5*g5p*lambda5*sina;
               double g5_vder2 = -g5_prefactor*weightij*g5p*g5v*2*expg5*nu5;
               addToValue( g5_ind, g5_nonweight*weightij, myvals ); 
               myvals.setSymfuncTemporyIndex( iind ); 
               addWeightDerivative( g5_ind, g5_nonweight*weightj, myvals );
               addVectorDerivatives( g5_ind, g5_vder1*dd_i, myvals ); 
               addVectorDerivatives(g5_ind, g5_vder2*disti, myvals );
               myvals.setSymfuncTemporyIndex( jind ); 
               addWeightDerivative( g5_ind, g5_nonweight*weighti, myvals );
               addVectorDerivatives( g5_ind, g5_vder1*dd_j, myvals ); 
               addVectorDerivatives(g5_ind, g5_vder2*distj, myvals );
           }
           // Compute G6
           if( g6_ind>-1 ) {
               double g6v = (1 + lambda6*cosa);
               double g6p = pow( g6v, zeta6-1 );
               double g6_nonweight = g6_prefactor*g6p*g6v;
               double g6_vder = -g6_prefactor*weightij*zeta6*g6p*lambda6*sina;
               addToValue( g6_ind, g6_nonweight*weightij, myvals );  
               myvals.setSymfuncTemporyIndex( iind ); 
               addWeightDerivative( g6_ind, g6_nonweight*weightj, myvals );
               addVectorDerivatives(g6_ind, g6_vder*dd_i, myvals );
               myvals.setSymfuncTemporyIndex( jind ); 
               addWeightDerivative( g6_ind, g6_nonweight*weighti, myvals );
               addVectorDerivatives(g6_ind, g6_vder*dd_j, myvals );
           }
           // Compute G7
           if( g7_ind>-1 ) {
               double g7_nonweight = 0.5*g7_prefactor*sin( nu7*(ang-alpha7) );
               double g7_vder = g7_prefactor*weightij*0.5*nu7*cos( nu7*(ang-alpha7) );
               addToValue( g7_ind, g7_nonweight*weightij, myvals );
               myvals.setSymfuncTemporyIndex( iind ); 
               addWeightDerivative( g7_ind, g7_nonweight*weightj, myvals );
               addVectorDerivatives( g7_ind, g7_vder*dd_i, myvals );
               myvals.setSymfuncTemporyIndex( jind ); 
               addWeightDerivative( g7_ind, g7_nonweight*weighti, myvals );
               addVectorDerivatives( g7_ind, g7_vder*dd_j, myvals );
           }
       }
   }
}

}
}
