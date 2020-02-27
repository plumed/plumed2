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
  std::vector<unsigned> ftypes;
  std::vector<double> lambda, nu, zeta, prefactor, alpha;
public:
  static void registerKeywords( Keywords& keys );
  explicit ThreeBodyGFunctions(const ActionOptions&);
  void compute( const double& weight, const Vector& vec, MultiValue& myvals ) const { plumed_error(); }
  void computeSymmetryFunction( const unsigned& current, MultiValue& myvals ) const ;
};

PLUMED_REGISTER_ACTION(ThreeBodyGFunctions,"GSYMFUNC_THREEBODY")

void ThreeBodyGFunctions::registerKeywords( Keywords& keys ) {
  SymmetryFunctionBase::registerKeywords( keys ); keys.remove("ONESHOT"); keys.remove("USECOLS");
  keys.add("numbered","FUNCTION","the parameters of the function you would like to compute");
  keys.add("compulsory","SWITCH","the switching function that acts on the distance between atom j and atom k in the G4 symmetry function");
  ActionWithValue::useCustomisableComponents( keys );
}

ThreeBodyGFunctions::ThreeBodyGFunctions(const ActionOptions&ao):
  Action(ao),
  SymmetryFunctionBase(ao)
{
  bool hasg4=false;
  for(int i=1;; i++) {
    std::string mystr, stype, lab, num; Tools::convert(i,num);
    if( !parseNumbered("FUNCTION",i,mystr ) ) break;
    std::vector<std::string> data=Tools::getWords(mystr);
    if( !Tools::parse(data,"LABEL",lab ) ) error("found no LABEL in FUNCTION" + num + " specification");
    if( !Tools::parse(data,"TYPE",stype) ) error("found no TYPE in FUNCTION" + num + " specification");
    lambda.push_back(0); nu.push_back(0); zeta.push_back(0); prefactor.push_back(0); alpha.push_back(0);
    addComponentWithDerivatives( lab );
    if( (stype=="g4" || stype=="g5" || stype=="g6") && !Tools::parse(data,"LAMBDA",lambda[i-1]) ) error("found no LAMBDA in FUNCTION" + num + " specification");
    if( (stype=="g4" || stype=="g5" || stype=="g6") && !Tools::parse(data,"ZETA",zeta[i-1]) ) error("found no ZETA in FUNCTION" + num + " specification");
    if( (stype=="g4" || stype=="g5" || stype=="g7") && !Tools::parse(data,"NU",nu[i-1]) ) error("found no NU in FUNCTION" + num + " specification");
    if( stype=="g7" && !Tools::parse(data,"ALPHA",alpha[i-1]) ) error("found no ALPHA in FUNCTION" + num + " specification");
    if( stype=="g4" ) {
      ftypes.push_back(4); prefactor[i-1] = 2.0 / pow(2.0,zeta[i-1]); hasg4=true;
      log.printf("  component labelled %s is of type g4 with parameters lambda=%f, nu=%f and zeta=%f\n", lab.c_str(), lambda[i-1], nu[i-1], zeta[i-1]);
    } else if( stype=="g5" ) {
      ftypes.push_back(5); prefactor[i-1] = 2.0 / pow(2.0,zeta[i-1]);
      log.printf("  component labelled %s is of type g5 with parameters lambda=%f, nu=%f and zeta=%f\n", lab.c_str(), lambda[i-1], nu[i-1], zeta[i-1]);
    } else if( stype=="g6" ) {
      ftypes.push_back(6); prefactor[i-1] = 2.0 / pow(2.0,zeta[i-1]);
      log.printf("  component labelled %s is of type g6 with parameters lambda=%f and zeta=%f\n", lab.c_str(), lambda[i-1], zeta[i-1]);
    } else if( stype=="g7" ) {
      ftypes.push_back(7); prefactor[i-1] = 2.0;
      log.printf("  component labelled %s is of type g7 with parameters nu=%f and alpha=%f\n", lab.c_str(), nu[i-1], alpha[i-1]);
    } else if( stype=="g8" ) {
      ftypes.push_back(8); prefactor[i-1] = 1.0;
      log.printf("  component labelled %s is of type g8 \n", lab.c_str() );
    }
  }
  if( hasg4 ) {
    std::string swstr, errors; parse("SWITCH",swstr); sf4.set(swstr, errors );
    if( errors.length()!=0 ) error("problem reading switching function description " + errors );
    log.printf("  using switching function distance between atom j and atom k with cutoff %s \n",sf4.description().c_str() );
  }
  checkRead();
}

void ThreeBodyGFunctions::computeSymmetryFunction( const unsigned& current, MultiValue& myvals ) const {
  bool hasg4=false; for(unsigned j=0; j<ftypes.size(); ++j) { if( ftypes[j]==4 ) { hasg4=true; break; } }

  Vector dd_i, dd_j, disti, distj;  double mod2i, mod2j, weighti, weightj; Angle angle;
  unsigned matind = getPntrToArgument(0)->getPositionInMatrixStash();
  unsigned matind_x = getPntrToArgument(1)->getPositionInMatrixStash();
  unsigned matind_y = getPntrToArgument(2)->getPositionInMatrixStash();
  unsigned matind_z = getPntrToArgument(3)->getPositionInMatrixStash();
  for(unsigned i=1; i<myvals.getNumberOfStashedMatrixElements(matind); ++i) {
    unsigned iind = myvals.getStashedMatrixIndex(matind,i);
    weighti = myvals.getStashedMatrixElement( matind, iind );
    if( weighti<epsilon ) continue ;
    disti[0] = myvals.getStashedMatrixElement( matind_x, iind );
    disti[1] = myvals.getStashedMatrixElement( matind_y, iind );
    disti[2] = myvals.getStashedMatrixElement( matind_z, iind );
    mod2i = disti.modulo2();
    for(unsigned j=0; j<i; ++j) {
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
      // Calculate switching function
      double dder_ij=0, switchij=0; if( hasg4 ) switchij = sf4.calculateSqr( mod2ij, dder_ij );
      // Now compute all symmetry functions
      for(unsigned n=0; n<ftypes.size(); ++n) {
        // Compute G4
        if( ftypes[n]==4 ) {
          double expg4 = exp( -nu[n]*( mod2i + mod2j + mod2ij ) );
          double g4v = (1 + lambda[n]*cosa);
          double g4p = pow( g4v, zeta[n]-1 );
          double g4_nonweight = prefactor[n]*g4p*g4v*expg4*switchij;
          double g4_vder1 = -prefactor[n]*weightij*zeta[n]*g4p*lambda[n]*sina*expg4*switchij;
          double g4_mder = prefactor[n]*weightij*g4p*g4v*expg4;
          double g4_vder2 = -g4_mder*2*nu[n]*switchij;
          double g4_vder3 = g4_mder*(dder_ij -2*nu[n]*switchij);
          addToValue( n, g4_nonweight*weightij, myvals );
          myvals.setSymfuncTemporyIndex( iind );
          addWeightDerivative( n, g4_nonweight*weightj, myvals );
          addVectorDerivatives(n, g4_vder1*dd_i, myvals );
          addVectorDerivatives(n, g4_vder2*disti, myvals );
          addVectorDerivatives(n, g4_vder3*dirij, myvals );
          myvals.setSymfuncTemporyIndex( jind );
          addWeightDerivative( n, g4_nonweight*weighti, myvals );
          addVectorDerivatives(n, g4_vder1*dd_j, myvals );
          addVectorDerivatives(n, g4_vder2*distj, myvals );
          addVectorDerivatives(n, -g4_vder3*dirij, myvals );
          // Compute G5
        } else if( ftypes[n]==5 ) {
          double expg5 = exp( -nu[n]*( mod2i + mod2j ) );
          double g5v = (1 + lambda[n]*cosa);
          double g5p = pow( g5v, zeta[n]-1 );
          double g5_nonweight = prefactor[n]*g5p*g5v*expg5;
          double g5_vder1 = -prefactor[n]*weightij*expg5*zeta[n]*g5p*lambda[n]*sina;
          double g5_vder2 = -prefactor[n]*weightij*g5p*g5v*2*expg5*nu[n];
          addToValue( n, g5_nonweight*weightij, myvals );
          myvals.setSymfuncTemporyIndex( iind );
          addWeightDerivative( n, g5_nonweight*weightj, myvals );
          addVectorDerivatives(n, g5_vder1*dd_i, myvals );
          addVectorDerivatives(n, g5_vder2*disti, myvals );
          myvals.setSymfuncTemporyIndex( jind );
          addWeightDerivative( n, g5_nonweight*weighti, myvals );
          addVectorDerivatives(n, g5_vder1*dd_j, myvals );
          addVectorDerivatives(n, g5_vder2*distj, myvals );
          // Compute G6
        } else if( ftypes[n]==6 ) {
          double g6v = (1 + lambda[n]*cosa);
          double g6p = pow( g6v, zeta[n]-1 );
          double g6_nonweight = prefactor[n]*g6p*g6v;
          double g6_vder = -prefactor[n]*weightij*zeta[n]*g6p*lambda[n]*sina;
          addToValue( n, g6_nonweight*weightij, myvals );
          myvals.setSymfuncTemporyIndex( iind );
          addWeightDerivative( n, g6_nonweight*weightj, myvals );
          addVectorDerivatives(n, g6_vder*dd_i, myvals );
          myvals.setSymfuncTemporyIndex( jind );
          addWeightDerivative( n, g6_nonweight*weighti, myvals );
          addVectorDerivatives(n, g6_vder*dd_j, myvals );
          // Compute G7
        } else if( ftypes[n]==7 ) {
          double g7_nonweight = 0.5*prefactor[n]*sin( nu[n]*(ang-alpha[n]) );
          double g7_vder = prefactor[n]*weightij*0.5*nu[n]*cos( nu[n]*(ang-alpha[n]) );
          addToValue( n, g7_nonweight*weightij, myvals );
          myvals.setSymfuncTemporyIndex( iind );
          addWeightDerivative( n, g7_nonweight*weightj, myvals );
          addVectorDerivatives( n, g7_vder*dd_i, myvals );
          myvals.setSymfuncTemporyIndex( jind );
          addWeightDerivative( n, g7_nonweight*weighti, myvals );
          addVectorDerivatives( n, g7_vder*dd_j, myvals );
        } else if( ftypes[n]==8 ) {
          double g8_vder = cos( ang ) + (1./3.);
          double g8_nonweight = g8_vder*g8_vder; 
          g8_vder = -2*g8_vder*weightij*sin( ang );
          addToValue( n, g8_nonweight*weightij, myvals );
          myvals.setSymfuncTemporyIndex( iind );
          addWeightDerivative( n, g8_nonweight*weightj, myvals );
          addVectorDerivatives( n, g8_vder*dd_i, myvals );
          myvals.setSymfuncTemporyIndex( jind );
          addWeightDerivative( n, g8_nonweight*weighti, myvals );
          addVectorDerivatives( n, g8_vder*dd_j, myvals ); 
        } else {
          plumed_merror("invalid symmetry function type");
        }
      }
    }
  }
}

}
}
