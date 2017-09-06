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
  std::vector<unsigned> ftypes;
  std::vector<double> center, nu, kappa;
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
  for(int i=1;;i++) {
      std::string mystr, stype, lab, num; Tools::convert(i,num);
      if( !parseNumbered("FUNCTION",i,mystr ) ) break;
      std::vector<std::string> data=Tools::getWords(mystr);
      if( !Tools::parse(data,"LABEL",lab ) ) error("found no LABEL in FUNCTION" + num + " specification");
      if( !Tools::parse(data,"TYPE",stype) ) error("found no TYPE in FUNCTION" + num + " specification");
      nu.push_back(0); center.push_back(0); kappa.push_back(0); addComponentWithDerivatives( lab );
      if( stype=="g1" ) {
          ftypes.push_back(1);
          log.printf("  component labelled %s is of type g1 \n",lab.c_str() );
      } else if( stype=="g2" ) {
          ftypes.push_back(2); 
          if( !Tools::parse(data,"CENTER",center[i-1]) ) error("found no CENTER in FUNCTION" + num + " specification");
          if( !Tools::parse(data,"NU",nu[i-1]) ) error("found no NU in FUNCTION" + num + " specification");
          log.printf("  component labelled %s is of type g2 with gaussian center at %f and width paramter of %f \n",lab.c_str(),center[i-1],nu[i-1]);
      } else if( stype=="g3" ) {
          ftypes.push_back(3);
          if( !Tools::parse(data,"KAPPA",kappa[i-1]) ) error("found no KAPPA in FUNCTION" + num + " specification");
          log.printf("  component labelled %s is of type g3 with kappa parameter equal %f \n",lab.c_str(), kappa[i-1]); 
      } else plumed_merror("found invalid type in FUNCTION" + num  + " specification should be g1, g2 or g3");
      if( !data.empty() ) {
          std::string errmsg = "found the following rogue keywords in FUNCTION" + num + " input : ";
          for(unsigned i=0; i<data.size(); ++i) errmsg = errmsg + data[i] + " ";
          error( errmsg ); 
      }
  }
}

void TwoBodyGFunctions::compute( const double& val, const Vector& distance, MultiValue& myvals ) const {
  double dlen = distance.modulo();
  for(unsigned i=0;i<ftypes.size();++i) {
      if( ftypes[i]==1 ) {
          addToValue( i, val, myvals ); addWeightDerivative( i, 1.0, myvals ); 
      } else if( ftypes[i]==2 ) {
          double diff = dlen - center[i], ee = exp( - nu[i]*diff*diff );
          addToValue( i, ee*val, myvals ); addWeightDerivative( i, ee, myvals ); 
          addVectorDerivatives( i, -(val/dlen)*2*nu[i]*diff*ee*distance, myvals );
      } else if( ftypes[i]==3 ) {
          double cc = cos( kappa[i]*dlen ), ss = -kappa[i]*sin( kappa[i]*dlen );
          addToValue( i, cc*val, myvals ); addWeightDerivative( i, cc, myvals );
          addVectorDerivatives( i, (val/dlen)*ss*distance, myvals );
      } else {
          plumed_merror("invalid gfunction type");
      }
  }
}

}
}
