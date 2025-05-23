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

//+PLUMEDOC MCOLVAR SPHERICAL_HARMONIC
/*
Calculate the values of all the spherical harmonic funtions for a particular value of l.

This action allows you to the all the [spherical harmonic](https://en.wikipedia.org/wiki/Spherical_harmonics) functions for a particular input
$L$ value.  As discussed in more detail in the article provided above the spherical harmonics this action calculates have the following general form:

$$
Y_l^m(\theta,\phi) = wNe^{im\phi} P_l^m(\cos\theta)
$$

where $N$ is a normalisation constant, $P_l^m$ is tn associated Legendre Polynomial and $w$ is an optional weight that can be passed in an input argument or simply set equal to one.
$e^{i\phi}$ and $\cos(\theta) are computed from the other input arguments, $x, y$ and $z$ as follows:

$$
e^{i\phi} = \frac{x}{r} + i\frac{y}{r} \qquad \textrm{and} \qquad \cos(\theta) = \frac{z}{r} \qquad \textrm{where} \qquad r = \sqrt{x^2 + y^2 + z^2}
$$

At present, the arguments for this action must be matrices. However, it would be easy to add functionality that would allow you to compute this function for scalar or vector input.
The arguments must all have the same shape as the two output components will also be matrices that are
calculated by applying the function above to each of the elements of the input matrix in turn.  The number of components output will be equal to $2(2L+1)$ and will contain
the real and imaginary parts of the $Y_l^m$ functions with the the $2l+1$ possible $m$ values.

The following intput provides an example that demonstrates how this function is used:

```plumed
d: DISTANCE_MATRIX GROUP=1-10 COMPONENTS
c: SPHERICAL_HARMONIC L=1 ARG=d.x,d.y,d.z
PRINT ARG=c.rm-n1 FILE=real_part_m-1
PRINT ARG=c.im-n1 FILE=imaginary_part_m-1
PRINT ARG=c.rm-0 FILE=real_part_m0
PRINT ARG=c.im-0 FILE=imaginary_part_m0
PRINT ARG=c.rm-p1 FILE=real_part_m+1
PRINT ARG=c.im-p1 FILE=imaginary_part_m+1
```

The DISTANCE_MATRIX command in the above input computes 3 $10\times10$ matrices.  These 3 $10\times10$ matrices are used in the input to the sphierical harmonic command,
which in turn outputs 6 $10\times10$ matrices that contain the real and imaginary parts when the three spherical harmonic functions with $l=1$ are applied element-wise to the above input. These six $10\times10$
matrices are then output to six separate files.

In the above example the weights for every distance is set equal to one.  The following example shows how an argument can be used to set the $w$ values to use when computing the function
above.

```plumed
s: CONTACT_MATRIX GROUP=1-10 SWITCH={RATIONAL R_0=1.0}
sc: CONTACT_MATRIX GROUP=1-10 SWITCH={RATIONAL R_0=1.0} COMPONENTS
c: SPHERICAL_HARMONIC L=1 ARG=sc.x,sc.y,sc.z,s
PRINT ARG=c.rm-n1 FILE=real_part_m-1
PRINT ARG=c.im-n1 FILE=imaginary_part_m-1
PRINT ARG=c.rm-0 FILE=real_part_m0
PRINT ARG=c.im-0 FILE=imaginary_part_m0
PRINT ARG=c.rm-p1 FILE=real_part_m+1
PRINT ARG=c.im-p1 FILE=imaginary_part_m+1
```

This function is used in the calculation of the Steinhardt order parameters, which are described in detail [here](https://www.plumed-tutorials.org/lessons/23/001/data/Steinhardt.html).

*/
//+ENDPLUMEDOC

//+PLUMEDOC MCOLVAR SPHERICAL_HARMONIC_MATRIX
/*
Calculate the values of all the spherical harmonic funtions for a particular value of l for all the elements of a set of three input matrices

\par Examples


*/
//+ENDPLUMEDOC

class SphericalHarmonic {
public:
  int tmom;
  std::vector<double> coeff_poly;
  std::vector<double> normaliz;
  static void registerKeywords( Keywords& keys );
  static unsigned factorial( const unsigned& n );
  static void read( SphericalHarmonic& func, ActionWithArguments* action, function::FunctionOptions& options );
  static void calc( const SphericalHarmonic& func, bool noderiv, const View<const double,helpers::dynamic_extent>& args, function::FunctionOutput& funcout );
  static double deriv_poly( const SphericalHarmonic& func, const unsigned& m, const double& val, double& df );
  static void addVectorDerivatives( const unsigned& ival, const Vector& der, View2D<double>& derivatives );
  SphericalHarmonic& operator=(const SphericalHarmonic& m) {
    tmom = m.tmom;
    coeff_poly = m.coeff_poly;
    normaliz = m.normaliz;
    return *this;
  }
};

typedef function::FunctionShortcut<SphericalHarmonic> SpHarmShortcut;
PLUMED_REGISTER_ACTION(SpHarmShortcut,"SPHERICAL_HARMONIC")
typedef function::FunctionOfMatrix<SphericalHarmonic> MatrixSpHarm;
PLUMED_REGISTER_ACTION(MatrixSpHarm,"SPHERICAL_HARMONIC_MATRIX")

void SphericalHarmonic::registerKeywords( Keywords& keys ) {
  keys.add("compulsory","L","the value of the angular momentum");
  keys.addOutputComponent("rm","default","matrix","the real parts of the spherical harmonic values with the m value given");
  keys.addOutputComponent("im","default","matrix","the real parts of the spherical harmonic values with the m value given");
  keys.add("hidden","MASKED_INPUT_ALLOWED","turns on that you are allowed to use masked inputs");
}

unsigned SphericalHarmonic::factorial( const unsigned& n ) {
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

void SphericalHarmonic::read( SphericalHarmonic& func, ActionWithArguments* action, function::FunctionOptions& options ) {
  action->parse("L",func.tmom);
  action->log.printf("  calculating %dth order spherical harmonic with %s, %s and %s as input \n", func.tmom, action->getPntrToArgument(0)->getName().c_str(), action->getPntrToArgument(1)->getName().c_str(), action->getPntrToArgument(2)->getName().c_str() );
  if( action->getNumberOfArguments()==4 ) {
    action->log.printf("  multiplying cylindrical harmonic by weight from %s \n", action->getPntrToArgument(3)->getName().c_str() );
  }

  func.normaliz.resize( func.tmom+1 );
  for(unsigned i=0; i<=func.tmom; ++i) {
    func.normaliz[i] = sqrt( (2*func.tmom+1)*factorial(func.tmom-i)/(4*pi*factorial(func.tmom+i)) );
    if( i%2==1 ) {
      func.normaliz[i]*=-1;
    }
  }

  func.coeff_poly.resize( func.tmom+1 );
  if( func.tmom==1 ) {
    // Legendre polynomial coefficients of order one
    func.coeff_poly[0]=0;
    func.coeff_poly[1]=1.0;
  } else if( func.tmom==2 ) {
    // Legendre polynomial coefficients of order two
    func.coeff_poly[0]=-0.5;
    func.coeff_poly[1]=0.0;
    func.coeff_poly[2]=1.5;
  } else if( func.tmom==3 ) {
    // Legendre polynomial coefficients of order three
    func.coeff_poly[0]=0.0;
    func.coeff_poly[1]=-1.5;
    func.coeff_poly[2]=0.0;
    func.coeff_poly[3]=2.5;
  } else if( func.tmom==4 ) {
    // Legendre polynomial coefficients of order four
    func.coeff_poly[0]=0.375;
    func.coeff_poly[1]=0.0;
    func.coeff_poly[2]=-3.75;
    func.coeff_poly[3]=0.0;
    func.coeff_poly[4]=4.375;
  } else if( func.tmom==5 ) {
    // Legendre polynomial coefficients of order five
    func.coeff_poly[0]=0.0;
    func.coeff_poly[1]=1.875;
    func.coeff_poly[2]=0.0;
    func.coeff_poly[3]=-8.75;
    func.coeff_poly[4]=0.0;
    func.coeff_poly[5]=7.875;
  } else if( func.tmom==6 ) {
    // Legendre polynomial coefficients of order six
    func.coeff_poly[0]=-0.3125;
    func.coeff_poly[1]=0.0;
    func.coeff_poly[2]=6.5625;
    func.coeff_poly[3]=0.0;
    func.coeff_poly[4]=-19.6875;
    func.coeff_poly[5]=0.0;
    func.coeff_poly[6]=14.4375;
  } else {
    action->error("Insert Legendre polynomial coefficients into SphericalHarmonics code");
  }
  std::string num;
  for(int i=-func.tmom; i<=func.tmom; ++i) {
    Tools::convert(fabs(i),num);
    if( i<0 ) {
      options.multipleValuesForEachRegisteredComponent.push_back( "-n" + num );
    } else if( i>0 ) {
      options.multipleValuesForEachRegisteredComponent.push_back( "-p" + num );
    } else {
      options.multipleValuesForEachRegisteredComponent.push_back( "-0");
    }
  }
  options.derivativeZeroIfValueIsZero = (action->getNumberOfArguments()==4 && (action->getPntrToArgument(3))->isDerivativeZeroWhenValueIsZero() );
}

void SphericalHarmonic::calc( const SphericalHarmonic& func, bool noderiv, const View<const double,helpers::dynamic_extent>& args, function::FunctionOutput& funcout ) {
  double weight=1;
  if( args.size()==4 ) {
    weight = args[3];
  }
  if( weight<epsilon ) {
    if( !noderiv ) {
      unsigned imbase = 2*func.tmom+1;
      for(int m=-func.tmom; m<=func.tmom; ++m) {
        funcout.derivs[func.tmom+m][0] = 0;
        funcout.derivs[func.tmom+m][1] = 0;
        funcout.derivs[func.tmom+m][2] = 0;
        funcout.derivs[func.tmom+m][3] = 0;
        funcout.derivs[imbase+func.tmom+m][0] = 0;
        funcout.derivs[imbase+func.tmom+m][1] = 0;
        funcout.derivs[imbase+func.tmom+m][2] = 0;
        funcout.derivs[imbase+func.tmom+m][3] = 0;
      }
    }
    return;
  }

  double dlen2 = args[0]*args[0]+args[1]*args[1]+args[2]*args[2];
  double dlen = sqrt( dlen2 );
  double dlen3 = dlen2*dlen;
  double dpoly_ass, poly_ass=deriv_poly( func, 0, args[2]/dlen, dpoly_ass );
  // Derivatives of z/r wrt x, y, z
  Vector dz;
  dz[0] = -( args[2] / dlen3 )*args[0];
  dz[1] = -( args[2] / dlen3 )*args[1];
  dz[2] = -( args[2] / dlen3 )*args[2] + (1.0 / dlen);
  // Accumulate for m=0
  funcout.values[func.tmom] = weight*poly_ass;
  addVectorDerivatives( func.tmom, weight*dpoly_ass*dz, funcout.derivs );
  if( args.size()==4 ) {
    funcout.derivs[func.tmom][3] = poly_ass;
  }

  // The complex number of which we have to take powers
  std::complex<double> com1( args[0]/dlen,args[1]/dlen ), dp_x, dp_y, dp_z;
  std::complex<double> powered = std::complex<double>(1.0,0.0);
  std::complex<double> ii( 0.0, 1.0 );
  Vector myrealvec, myimagvec, real_dz, imag_dz;
  // Do stuff for all other m values
  for(unsigned m=1; m<=func.tmom; ++m) {
    // Calculate Legendre Polynomial
    poly_ass=deriv_poly( func, m, args[2]/dlen, dpoly_ass );
    // Real and imaginary parts of z
    double real_z = real(com1*powered), imag_z = imag(com1*powered);

    // Calculate steinhardt parameter
    double tq6=poly_ass*real_z;   // Real part of steinhardt parameter
    double itq6=poly_ass*imag_z;  // Imaginary part of steinhardt parameter

    // Derivatives wrt ( x/r + iy )^m
    double md=static_cast<double>(m);
    dp_x = md*powered*( (1.0/dlen)-(args[0]*args[0])/dlen3-ii*(args[0]*args[1])/dlen3 );
    dp_y = md*powered*( ii*(1.0/dlen)-(args[0]*args[1])/dlen3-ii*(args[1]*args[1])/dlen3 );
    dp_z = md*powered*( -(args[0]*args[2])/dlen3-ii*(args[1]*args[2])/dlen3 );

    // Derivatives of real and imaginary parts of above
    real_dz[0] = real( dp_x );
    real_dz[1] = real( dp_y );
    real_dz[2] = real( dp_z );
    imag_dz[0] = imag( dp_x );
    imag_dz[1] = imag( dp_y );
    imag_dz[2] = imag( dp_z );

    // Complete derivative of steinhardt parameter
    myrealvec = weight*dpoly_ass*real_z*dz + weight*poly_ass*real_dz;
    myimagvec = weight*dpoly_ass*imag_z*dz + weight*poly_ass*imag_dz;

    // Real part
    funcout.values[func.tmom+m] = weight*tq6;
    // Imaginary part
    funcout.values[3*func.tmom+1+m] = weight*itq6;
    if( !noderiv ) {
      addVectorDerivatives( func.tmom+m, myrealvec, funcout.derivs );
      addVectorDerivatives( 3*func.tmom+1+m, myimagvec, funcout.derivs );
    }
    // Store -m part of vector
    double pref=pow(-1.0,m);
    // -m part of vector is just +m part multiplied by (-1.0)**m and multiplied by complex
    // conjugate of Legendre polynomial
    // Real part
    funcout.values[func.tmom-m] = pref*weight*tq6;
    // Imaginary part
    funcout.values[3*func.tmom+1-m] = -pref*weight*itq6;
    if( !noderiv ) {
      addVectorDerivatives( func.tmom-m, pref*myrealvec, funcout.derivs );
      addVectorDerivatives( 3*func.tmom+1-m, -pref*myimagvec, funcout.derivs );
      if( args.size()==4 ) {
        funcout.derivs[func.tmom+m][3]=tq6;
        funcout.derivs[3*func.tmom+1+m][3]=itq6;
        funcout.derivs[func.tmom-m][3]=pref*tq6;
        funcout.derivs[3*func.tmom+1-m][3]=-pref*itq6;
      }
    }
    // Calculate next power of complex number
    powered *= com1;
  }
}

double SphericalHarmonic::deriv_poly( const SphericalHarmonic& func, const unsigned& m, const double& val, double& df ) {
  double fact=1.0;
  for(unsigned j=1; j<=m; ++j) {
    fact=fact*j;
  }
  double res=func.coeff_poly[m]*fact;

  double pow=1.0, xi=val, dxi=1.0;
  df=0.0;
  for(int i=m+1; i<=func.tmom; ++i) {
    fact=1.0;
    for(unsigned j=i-m+1; j<=i; ++j) {
      fact=fact*j;
    }
    res=res+func.coeff_poly[i]*fact*xi;
    df = df + pow*func.coeff_poly[i]*fact*dxi;
    xi=xi*val;
    dxi=dxi*val;
    pow+=1.0;
  }
  df = df*func.normaliz[m];
  return func.normaliz[m]*res;
}

void SphericalHarmonic::addVectorDerivatives( const unsigned& ival, const Vector& der, View2D<double>& derivatives ) {
  for(unsigned j=0; j<3; ++j) {
    derivatives[ival][j] = der[j];
  }
}

}
}

