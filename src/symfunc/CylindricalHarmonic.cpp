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
#include "function/Function.h"
#include "core/ActionRegister.h"

#include <complex>

using namespace std;

namespace PLMD {
namespace symfunc {

//+PLUMEDOC MCOLVAR CYLINDRICAL_HARMONIC
/*

\par Examples


*/
//+ENDPLUMEDOC


class CylindricalHarmonic : public function::Function {
private:
  int tmom;
public:
  static void registerKeywords( Keywords& keys );
  explicit CylindricalHarmonic(const ActionOptions&);
  void calculateFunction( const std::vector<double>& args, MultiValue& myvals ) const override;
};

PLUMED_REGISTER_ACTION(CylindricalHarmonic,"CYLINDRICAL_HARMONIC")

void CylindricalHarmonic::registerKeywords( Keywords& keys ) {
  function::Function::registerKeywords( keys ); keys.use("ARG");
  keys.add("compulsory","DEGREE","the value of the n parameter in the equation above");
  keys.addOutputComponent("rm","default","the real part of the cylindrical harmonic");
  keys.addOutputComponent("im","default","the imaginary part of the cylindrical harmonic");
}

CylindricalHarmonic::CylindricalHarmonic(const ActionOptions&ao):
  Action(ao),
  Function(ao)
{
  parse("DEGREE",tmom);
  log.printf("  calculating %dth order cylindrical harmonic with %s and %s as input \n", tmom, getPntrToArgument(0)->getName().c_str(), getPntrToArgument(1)->getName().c_str() );
  if( getNumberOfArguments()==3 ) log.printf("  multiplying cylindrical harmonic by weight from %s \n", getPntrToArgument(2)->getName().c_str() );
  addComponentWithDerivatives( "rm" ); addComponentWithDerivatives( "im" );
  componentIsNotPeriodic( "rm" ); componentIsNotPeriodic( "im" );
  checkRead();
}

void CylindricalHarmonic::calculateFunction( const std::vector<double>& args, MultiValue& myvals ) const {
  double dlen2 = args[0]*args[0] + args[1]*args[1]; double dlen = sqrt( dlen2 ); double dlen3 = dlen2*dlen; 
  std::complex<double> com1( args[0]/dlen,args[1]/dlen ); double weight=1; if( args.size()==3 ) weight=args[2];
  std::complex<double> ppp = pow( com1, tmom-1 ), ii( 0, 1 ); double real_z = real( ppp*com1 ), imag_z = imag( ppp*com1 );
  std::complex<double> dp_x = static_cast<double>(tmom)*ppp*( (1.0/dlen)-(args[0]*args[0])/dlen3-ii*(args[0]*args[1])/dlen3 );
  std::complex<double> dp_y = static_cast<double>(tmom)*ppp*( ii*(1.0/dlen)-(args[0]*args[1])/dlen3-ii*(args[1]*args[1])/dlen3 );
  addValue( 0, weight*real_z, myvals ); addDerivative( 0, 0, weight*real(dp_x), myvals ); addDerivative( 0, 1, weight*real(dp_y), myvals ); 
  addValue( 1, weight*imag_z, myvals ); addDerivative( 1, 0, weight*imag(dp_x), myvals ); addDerivative( 1, 1, weight*imag(dp_y), myvals ); 
  if( args.size()==3 ) { addDerivative( 0, 2, real_z, myvals ); addDerivative( 1, 2, imag_z, myvals ); }
}

}
}

