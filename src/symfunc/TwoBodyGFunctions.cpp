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
#include "tools/SwitchingFunction.h"
#include <string>
#include <cmath>

using namespace std;

namespace PLMD {
namespace symfunc {


class TwoBodyGFunctions : public SymmetryFunctionBase {
private:
  std::vector<SwitchingFunction> functions;
public:
  static void registerKeywords( Keywords& keys );
  explicit TwoBodyGFunctions(const ActionOptions&);
  void compute( const double& weight, const Vector& vec, MultiValue& myvals ) const ;
};

PLUMED_REGISTER_ACTION(TwoBodyGFunctions,"GSYMFUNC_TWOBODY")

void TwoBodyGFunctions::registerKeywords( Keywords& keys ) {
  SymmetryFunctionBase::registerKeywords( keys );
  keys.add("numbered","FUNCTION","the parameters of the function you would like to compute");
  ActionWithValue::useCustomisableComponents( keys );
}

TwoBodyGFunctions::TwoBodyGFunctions(const ActionOptions&ao):
  Action(ao),
  SymmetryFunctionBase(ao)
{
  for(int i=1;; i++) {
    std::string mystr, myfunc, errors, stype, lab, num; Tools::convert(i,num);
    if( !parseNumbered("FUNCTION",i,mystr ) ) break;
    std::vector<std::string> data=Tools::getWords(mystr);
    if( !Tools::parse(data,"LABEL",lab ) ) error("found no LABEL in FUNCTION" + num + " specification");
    addComponentWithDerivatives( lab ); 
    if( !Tools::parse(data,"FUNC",myfunc) ) break;
    functions.push_back( SwitchingFunction() );
    functions[i-1].set( "CUSTOM FUNC=" + myfunc + " R_0=1.0", errors );
    if( errors.length()>0 ) error("errors reading function " + errors ); 
    log.printf("  component labelled %s is %s \n",lab.c_str(), functions[i-1].description().c_str() );
  }
}

void TwoBodyGFunctions::compute( const double& val, const Vector& distance, MultiValue& myvals ) const {
  double dlen = distance.modulo();
  for(unsigned i=0; i<functions.size(); ++i) {
    double df, fval = functions[i].calculate( dlen, df ); addToValue( i, val*fval, myvals );
    addWeightDerivative( i, fval, myvals ); addVectorDerivatives( i, df*val*distance, myvals ); 
  }
}

}
}
