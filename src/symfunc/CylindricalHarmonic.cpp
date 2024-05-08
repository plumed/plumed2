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
#include "function/FunctionTemplateBase.h"
#include "function/FunctionShortcut.h"
#include "function/FunctionOfMatrix.h"
#include "core/ActionRegister.h"

#include <complex>

namespace PLMD {
namespace symfunc {

//+PLUMEDOC MCOLVAR CYLINDRICAL_HARMONIC
/*
Calculate the cylindrical harmonic function

\par Examples


*/
//+ENDPLUMEDOC

//+PLUMEDOC MCOLVAR CYLINDRICAL_HARMONIC_MATRIX
/*
Calculate the cylindrical harmonic function from the elements in two input matrices

\par Examples


*/
//+ENDPLUMEDOC


class CylindricalHarmonic : public function::FunctionTemplateBase {
private:
  int tmom;
public:
  void registerKeywords( Keywords& keys ) override;
  void read( ActionWithArguments* action ) override;
  void setPeriodicityForOutputs( ActionWithValue* action ) override;
  void calc( const ActionWithArguments* action, const std::vector<double>& args, std::vector<double>& vals, Matrix<double>& derivatives ) const override;
};

typedef function::FunctionShortcut<CylindricalHarmonic> CyHarmShortcut;
PLUMED_REGISTER_ACTION(CyHarmShortcut,"CYLINDRICAL_HARMONIC")
typedef function::FunctionOfMatrix<CylindricalHarmonic> MatrixCyHarm;
PLUMED_REGISTER_ACTION(MatrixCyHarm,"CYLINDRICAL_HARMONIC_MATRIX")

void CylindricalHarmonic::registerKeywords( Keywords& keys ) {
  keys.add("compulsory","DEGREE","the value of the n parameter in the equation above");
  keys.addOutputComponent("rm","default","the real part of the cylindrical harmonic");
  keys.addOutputComponent("im","default","the imaginary part of the cylindrical harmonic");
}

void CylindricalHarmonic::read( ActionWithArguments* action ) {
  parse(action,"DEGREE",tmom);
  action->log.printf("  calculating %dth order cylindrical harmonic with %s and %s as input \n", tmom, action->getPntrToArgument(0)->getName().c_str(), action->getPntrToArgument(1)->getName().c_str() );
  if( action->getNumberOfArguments()==3 ) action->log.printf("  multiplying cylindrical harmonic by weight from %s \n", action->getPntrToArgument(2)->getName().c_str() );
}

void CylindricalHarmonic::setPeriodicityForOutputs( ActionWithValue* action ) {
  action->componentIsNotPeriodic("rm"); action->componentIsNotPeriodic("im");
}

void CylindricalHarmonic::calc( const ActionWithArguments* action, const std::vector<double>& args, std::vector<double>& vals, Matrix<double>& derivatives ) const {
  double dlen2 = args[0]*args[0] + args[1]*args[1]; double dlen = sqrt( dlen2 ); double dlen3 = dlen2*dlen;
  std::complex<double> com1( args[0]/dlen,args[1]/dlen ); double weight=1; if( args.size()==3 ) weight=args[2];
  std::complex<double> ppp = pow( com1, tmom-1 ), ii( 0, 1 ); double real_z = real( ppp*com1 ), imag_z = imag( ppp*com1 );
  std::complex<double> dp_x = static_cast<double>(tmom)*ppp*( (1.0/dlen)-(args[0]*args[0])/dlen3-ii*(args[0]*args[1])/dlen3 );
  std::complex<double> dp_y = static_cast<double>(tmom)*ppp*( ii*(1.0/dlen)-(args[0]*args[1])/dlen3-ii*(args[1]*args[1])/dlen3 );
  vals[0] = weight*real_z; derivatives(0,0) = weight*real(dp_x); derivatives(0,1) = weight*real(dp_y);
  vals[1] = weight*imag_z; derivatives(1,0) = weight*imag(dp_x); derivatives(1,1) = weight*imag(dp_y);
  if( args.size()==3 ) { derivatives(0,2) = real_z; derivatives(1,2) = imag_z; }
}

}
}

