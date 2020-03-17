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

#include <complex>

using namespace std;

namespace PLMD {
namespace symfunc {

//+PLUMEDOC MCOLVAR CYLINDRICAL_HARMONIC
/*

\par Examples


*/
//+ENDPLUMEDOC


class CylindricalHarmonic : public SymmetryFunctionBase {
private:
  int tmom;
public:
  static void registerKeywords( Keywords& keys );
  explicit CylindricalHarmonic(const ActionOptions&);
  void compute( const double& val, const Vector& dir, MultiValue& myvals ) const ;
};

PLUMED_REGISTER_ACTION(CylindricalHarmonic,"CYLINDRICAL_HARMONIC")

void CylindricalHarmonic::registerKeywords( Keywords& keys ) {
  SymmetryFunctionBase::registerKeywords( keys );
  keys.add("compulsory","DEGREE","the value of the n parameter in the equation above");
  keys.addOutputComponent("rm","default","the real part of the cylindrical harmonic");
  keys.addOutputComponent("im","default","the imaginary part of the cylindrical harmonic");
}

CylindricalHarmonic::CylindricalHarmonic(const ActionOptions&ao):
  Action(ao),
  SymmetryFunctionBase(ao)
{
  parse("DEGREE",tmom);
  log.printf("  calculating %dth order cylindrical harmonic \n", tmom);
  addComponentWithDerivatives( "rm" ); addComponentWithDerivatives( "im" );
  checkRead();
}

void CylindricalHarmonic::compute( const double& val, const Vector& distance, MultiValue& myvals ) const {
  double dlen2 = distance.modulo2(); double dlen = sqrt( dlen2 ); double dlen3 = dlen2*dlen; 
  std::complex<double> com1( distance[0]/dlen,distance[1]/dlen );
  std::complex<double> ppp = pow( com1, tmom-1 ), ii( 0, 1 );
  std::complex<double> dp_x = static_cast<double>(tmom)*ppp*( (1.0/dlen)-(distance[0]*distance[0])/dlen3-ii*(distance[0]*distance[1])/dlen3 );
  std::complex<double> dp_y = static_cast<double>(tmom)*ppp*( ii*(1.0/dlen)-(distance[0]*distance[1])/dlen3-ii*(distance[1]*distance[1])/dlen3 );
  double real_z = real( ppp*com1 ); double imag_z = imag( ppp*com1 );
  Vector myrealvec, myimagvec; myrealvec.zero(); myimagvec.zero();
  myrealvec[0] = real( dp_x ); myrealvec[1] = real( dp_y );
  myimagvec[0] = imag( dp_x ); myimagvec[1] = imag( dp_y );
  addToValue( 0, val*real_z, myvals ); 
  addVectorDerivatives( 0, val*myrealvec, myvals );
  addWeightDerivative( 0, real_z, myvals );
  addToValue( 1, val*imag_z, myvals );  
  addVectorDerivatives( 1, val*myimagvec, myvals );
  addWeightDerivative( 1, imag_z, myvals );
}

}
}

