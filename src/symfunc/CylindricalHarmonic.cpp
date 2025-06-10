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
#include "function/FunctionSetup.h"
#include "function/FunctionShortcut.h"
#include "function/FunctionOfMatrix.h"
#include "core/ActionRegister.h"

#include <complex>

namespace PLMD {
namespace symfunc {

//+PLUMEDOC MCOLVAR CYLINDRICAL_HARMONIC
/*
Calculate the cylindrical harmonic function

This action allows you to the value of the following complex function.  The action outputs
two components that are the real and imaginary parts of the following function:

$$
z = w (\frac{x}{r} + \frac{y}{r} i )^n \qquad \textrm{where} \qquad r = \sqrt(x^2 + y^2}
$$

In this expression $n$ is a parameter that is specified using the DEGREE keyword. $x$ and $y$ are the input arguments and $w$ is an optional input weight, which is set equal to
one if only two arguments are provided in input.  At present, the arguments for this action must be matrices.
These arguments must all have the same shape as the two output components will also be matrices that are
calculated by applying the function above to each of the elements of the input matrix in turn.

The following intput provides an example that demonstrates how this function is used:

```plumed
d: DISTANCE_MATRIX GROUP=1-10 COMPONENTS
c: CYLINDRICAL_HARMONIC DEGREE=6 ARG=d.x,d.y
PRINT ARG=c.rm FILE=real_part
PRINT ARG=c.im FILE=imaginary_part
```

The DISTANCE_MATRIX command in the above input computes 3 $10\times10$ matrices.  Two of these $10\times10$ matrices are used in the input to the cylindrical harmonic command,
which in turn outputs two $10\times10$ matrices that contain the real and imaginary parts when the function above is applied element-wise to the above input. These two $10\times10$
matrices are then output to two separate files.

In the above example the weights for every distance is set equal to one.  The following example shows how an argument can be used to set the $w$ values to use when computing the function
above.

```plumed
s: CONTACT_MATRIX GROUP=1-10 SWITCH={RATIONAL R_0=1.0}
sc: CONTACT_MATRIX GROUP=1-10 SWITCH={RATIONAL R_0=1.0} COMPONENTS
c: CYLINDRICAL_HARMONIC DEGREE=6 ARG=sc.x,sc.y,s
PRINT ARG=c.rm FILE=real_part
PRINT ARG=c.im FILE=imaginary_part
```

*/
//+ENDPLUMEDOC


class CylindricalHarmonic {
public:
  int tmom;
  static void registerKeywords( Keywords& keys );
  static void read( CylindricalHarmonic& func, ActionWithArguments* action, function::FunctionOptions& options );
  static void calc( const CylindricalHarmonic& func, bool noderiv, const View<const double,helpers::dynamic_extent>& args, function::FunctionOutput& funcout );
  CylindricalHarmonic& operator=(const CylindricalHarmonic& m) {
    tmom=m.tmom;
    return *this;
  }
};

typedef function::FunctionShortcut<CylindricalHarmonic> CyHarmShortcut;
PLUMED_REGISTER_ACTION(CyHarmShortcut,"CYLINDRICAL_HARMONIC")
typedef function::FunctionOfMatrix<CylindricalHarmonic> MatrixCyHarm;
PLUMED_REGISTER_ACTION(MatrixCyHarm,"CYLINDRICAL_HARMONIC_MATRIX")

void CylindricalHarmonic::registerKeywords( Keywords& keys ) {
  keys.add("compulsory","DEGREE","the value of the n parameter in the equation above");
  keys.addOutputComponent("rm","default","matrix","the real part of the cylindrical harmonic");
  keys.addOutputComponent("im","default","matrix","the imaginary part of the cylindrical harmonic");
}

void CylindricalHarmonic::read( CylindricalHarmonic& func, ActionWithArguments* action, function::FunctionOptions& options ) {
  action->parse("DEGREE",func.tmom);
  action->log.printf("  calculating %dth order cylindrical harmonic with %s and %s as input \n", func.tmom, action->getPntrToArgument(0)->getName().c_str(), action->getPntrToArgument(1)->getName().c_str() );
  if( action->getNumberOfArguments()==3 ) {
    action->log.printf("  multiplying cylindrical harmonic by weight from %s \n", action->getPntrToArgument(2)->getName().c_str() );
  }
  options.derivativeZeroIfValueIsZero = (action->getNumberOfArguments()==3 && (action->getPntrToArgument(2))->isDerivativeZeroWhenValueIsZero());
}

void CylindricalHarmonic::calc( const CylindricalHarmonic& func, bool noderiv, const View<const double,helpers::dynamic_extent>& args, function::FunctionOutput& funcout ) {
  double dlen2 = args[0]*args[0] + args[1]*args[1];
  double dlen = sqrt( dlen2 );
  double dlen3 = dlen2*dlen;
  std::complex<double> com1( args[0]/dlen,args[1]/dlen );
  double weight=1;
  if( args.size()==3 ) {
    weight=args[2];
  }
  std::complex<double> ppp = pow( com1, func.tmom-1 ), ii( 0, 1 );
  double real_z = real( ppp*com1 ), imag_z = imag( ppp*com1 );
  std::complex<double> dp_x = static_cast<double>(func.tmom)*ppp*( (1.0/dlen)-(args[0]*args[0])/dlen3-ii*(args[0]*args[1])/dlen3 );
  std::complex<double> dp_y = static_cast<double>(func.tmom)*ppp*( ii*(1.0/dlen)-(args[0]*args[1])/dlen3-ii*(args[1]*args[1])/dlen3 );
  funcout.values[0] = weight*real_z;
  funcout.values[1] = weight*imag_z;

  if( !noderiv ) {
    funcout.derivs[0][0] = weight*real(dp_x);
    funcout.derivs[0][1] = weight*real(dp_y);
    funcout.derivs[1][0] = weight*imag(dp_x);
    funcout.derivs[1][1] = weight*imag(dp_y);
    if( args.size()==3 ) {
      funcout.derivs[0][2] = real_z;
      funcout.derivs[1][2] = imag_z;
    }
  }
}

}
}

