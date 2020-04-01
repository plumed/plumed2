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
#include "tools/LeptonCall.h"
#include "tools/Angle.h"

namespace PLMD {
namespace symfunc {

class ThreeBodyGFunctions : public SymmetryFunctionBase {
private:
  std::vector<LeptonCall> functions;
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
  for(int i=1;; i++) {
    std::string myfunc, mystr, lab, num; Tools::convert(i,num);
    if( !parseNumbered("FUNCTION",i,mystr ) ) break;
    std::vector<std::string> data=Tools::getWords(mystr);
    if( !Tools::parse(data,"LABEL",lab ) ) error("found no LABEL in FUNCTION" + num + " specification");
    addComponentWithDerivatives( lab );
    if( !Tools::parse(data,"FUNC",myfunc) ) error("found no FUNC in FUNCTION" + num + " specification");
    log.printf("  component labelled %s is computed using %s \n",lab.c_str(), myfunc.c_str() ); 
    functions.push_back( LeptonCall() ); std::vector<std::string> argnames(1); argnames[0]="ajik";
    if( myfunc.find("rij")!=std::string::npos ) argnames.push_back("rij"); 
    if( myfunc.find("rik")!=std::string::npos ) { 
        if( argnames.size()<2 ) error("if you have a function of rik it must also be a function of rij -- email gareth.tribello@gmail.com if this is a problem");
        argnames.push_back("rik"); 
    }
    if( myfunc.find("rjk")!=std::string::npos ) {
        if( argnames.size()<2 ) error("if you have a function of rjk it must also be a function of rij and rik -- email gareth.tribello@gmail.com if this is a problem");
        argnames.push_back("rjk"); 
    }
    functions[i-1].set( myfunc, argnames, this, true ); 
  }
  checkRead();
}

void ThreeBodyGFunctions::computeSymmetryFunction( const unsigned& current, MultiValue& myvals ) const {
  std::vector<Vector> der_i(4), der_j(4);
  std::vector<double> values(4); Angle angle; Vector disti, distj; 
  unsigned matind = getPntrToArgument(0)->getPositionInMatrixStash();
  unsigned matind_x = getPntrToArgument(1)->getPositionInMatrixStash();
  unsigned matind_y = getPntrToArgument(2)->getPositionInMatrixStash();
  unsigned matind_z = getPntrToArgument(3)->getPositionInMatrixStash();
  for(unsigned i=1; i<myvals.getNumberOfStashedMatrixElements(matind); ++i) {
    unsigned iind = myvals.getStashedMatrixIndex(matind,i);
    double weighti = myvals.getStashedMatrixElement( matind, iind );
    if( weighti<epsilon ) continue ;
    disti[0] = myvals.getStashedMatrixElement( matind_x, iind );
    disti[1] = myvals.getStashedMatrixElement( matind_y, iind );
    disti[2] = myvals.getStashedMatrixElement( matind_z, iind );
    values[1] = disti.modulo2(); der_i[1]=2*disti; der_i[2].zero();
    for(unsigned j=0; j<i; ++j) {
      unsigned jind = myvals.getStashedMatrixIndex(matind,j);
      double weightj = myvals.getStashedMatrixElement( matind, jind );
      if( weightj<epsilon ) continue ;
      distj[0] = myvals.getStashedMatrixElement( matind_x, jind );
      distj[1] = myvals.getStashedMatrixElement( matind_y, jind );
      distj[2] = myvals.getStashedMatrixElement( matind_z, jind );
      values[2] = distj.modulo2(); der_j[1].zero(); der_j[2]=2*distj;
      der_i[3] = disti - distj; values[3] = der_i[3].modulo2(); 
      der_i[3] = 2*der_i[3]; der_j[3] = -der_i[3];
      // Compute angle between bonds
      values[0] = angle.compute( disti, distj, der_i[0], der_j[0] );
      // Compute product of weights
      double weightij = weighti*weightj;
      // Now compute all symmetry functions
      for(unsigned n=0; n<functions.size(); ++n) {
        double nonweight = functions[n].evaluate( values );
        addToValue( n, nonweight*weightij, myvals );  
        for(unsigned m=0;m<functions[n].getNumberOfArguments();++m) {
            double der = weightij*functions[n].evaluateDeriv( m, values );
            myvals.setSymfuncTemporyIndex( iind );
            addVectorDerivatives(n, der*der_i[m], myvals );
            myvals.setSymfuncTemporyIndex( jind );
            addVectorDerivatives(n, der*der_j[m], myvals );
        }
        myvals.setSymfuncTemporyIndex( iind );
        addWeightDerivative( n, nonweight*weightj, myvals );
        myvals.setSymfuncTemporyIndex( jind );   
        addWeightDerivative( n, nonweight*weighti, myvals );
      }
    }
  }
}

}
}
