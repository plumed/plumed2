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
#ifdef __PLUMED_HAS_OPENACC
#define __PLUMED_USE_OPENACC 1
#endif //__PLUMED_HAS_OPENACC
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

class SphericalHarmonic {
  int tmom=0;
  std::vector<double> coeff_poly_v{};
  std::vector<double> normaliz_v{};
  double *coeff_poly=nullptr;
  double *normaliz=nullptr;
public:
  static void registerKeywords( Keywords& keys );
  static unsigned factorial( unsigned n );
  static void read( SphericalHarmonic& func,
                    ActionWithArguments* action,
                    function::FunctionOptions& options );
  static void calc( const SphericalHarmonic& func,
                    bool noderiv,
                    View<const double> args,
                    function::FunctionOutput funcout );
  double deriv_poly(                            unsigned m,
      double val,
      double& df ) const;
  static inline void addVectorDerivatives( unsigned ival,
      const Vector& der,
      View2D<double> derivatives );
  void update() {
    coeff_poly = coeff_poly_v.data();
    normaliz = normaliz_v.data();
  }
  SphericalHarmonic() = default;
  ~SphericalHarmonic() = default;
  SphericalHarmonic(const SphericalHarmonic&x):
    tmom(x.tmom),
    coeff_poly_v(x.coeff_poly_v),
    normaliz_v(x.normaliz_v) {
    update();
  }
  SphericalHarmonic(SphericalHarmonic&&x):
    tmom(x.tmom),
    coeff_poly_v(std::move(x.coeff_poly_v)),
    normaliz_v(std::move(x.normaliz_v)) {
    update();
  }
  SphericalHarmonic &operator=(const SphericalHarmonic&x) {
    if (this!=&x) {
      tmom=x.tmom;
      coeff_poly_v=x.coeff_poly_v;
      normaliz_v=x.normaliz_v;
      update();
    }
    return *this;
  }
  SphericalHarmonic &operator=(SphericalHarmonic&&x) {
    if (this!=&x) {
      tmom=x.tmom;
      coeff_poly_v=std::move(x.coeff_poly_v);
      normaliz_v=std::move(x.normaliz_v);
      update();
    }
    return *this;
  }
#ifdef __PLUMED_USE_OPENACC
  void toACCDevice() const {
    const auto num=tmom+1;
#pragma acc enter data copyin(this[0:1], \
  tmom, coeff_poly[0:num], normaliz[0:num])
  }
  void removeFromACCDevice() const {
    const auto num=tmom+1;
#pragma acc exit data delete(normaliz[0:num], \
  coeff_poly[0:num], tmom, this[0:1])
  }
#endif // __PLUMED_USE_OPENACC
};

typedef function::FunctionShortcut<SphericalHarmonic> SpHarmShortcut;
PLUMED_REGISTER_ACTION(SpHarmShortcut,"SPHERICAL_HARMONIC")
typedef function::FunctionOfMatrix<SphericalHarmonic> MatrixSpHarm;
PLUMED_REGISTER_ACTION(MatrixSpHarm,"SPHERICAL_HARMONIC_MATRIX")

void SphericalHarmonic::registerKeywords( Keywords& keys ) {
  keys.add("compulsory","L","the value of the angular momentum");
  keys.addOutputComponent("rm","default",
                          "matrix","the real parts of the spherical harmonic values with the m value given");
  keys.addOutputComponent("im","default",
                          "matrix","the real parts of the spherical harmonic values with the m value given");
  keys.add("hidden","MASKED_INPUT_ALLOWED",
           "turns on that you are allowed to use masked inputs");
}

unsigned SphericalHarmonic::factorial( unsigned n ) {
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

void SphericalHarmonic::read( SphericalHarmonic& func,
                              ActionWithArguments* action,
                              function::FunctionOptions& options ) {
  action->parse("L",func.tmom);
  action->log.printf("  calculating %dth order spherical harmonic with %s, %s and %s as input \n",
                     func.tmom,
                     action->getPntrToArgument(0)->getName().c_str(),
                     action->getPntrToArgument(1)->getName().c_str(),
                     action->getPntrToArgument(2)->getName().c_str() );
  if( action->getNumberOfArguments()==4 ) {
    action->log.printf("  multiplying cylindrical harmonic by weight from %s \n",
                       action->getPntrToArgument(3)->getName().c_str() );
  }

  func.normaliz_v.resize( func.tmom+1 );
  func.coeff_poly_v.resize( func.tmom+1 );
  func.update();

  for(int i=0; i<=func.tmom; ++i) {
    func.normaliz[i] = ( i%2==0 ? 1.0 : -1.0 )
                       * sqrt( (2*func.tmom+1)
                               * factorial(func.tmom-i)
                               /(4*PLMD::pi*factorial(func.tmom+i)) );
  }

  switch (func.tmom) {
  case 1:
    // Legendre polynomial coefficients of order one
    func.coeff_poly[0]=0;
    func.coeff_poly[1]=1.0;
    break;
  case 2:
    // Legendre polynomial coefficients of order two
    func.coeff_poly[0]=-0.5;
    func.coeff_poly[1]=0.0;
    func.coeff_poly[2]=1.5;
    break;
  case 3:
    // Legendre polynomial coefficients of order three
    func.coeff_poly[0]=0.0;
    func.coeff_poly[1]=-1.5;
    func.coeff_poly[2]=0.0;
    func.coeff_poly[3]=2.5;
    break;
  case 4:
    // Legendre polynomial coefficients of order four
    func.coeff_poly[0]=0.375;
    func.coeff_poly[1]=0.0;
    func.coeff_poly[2]=-3.75;
    func.coeff_poly[3]=0.0;
    func.coeff_poly[4]=4.375;
    break;
  case 5:
    // Legendre polynomial coefficients of order five
    func.coeff_poly[0]=0.0;
    func.coeff_poly[1]=1.875;
    func.coeff_poly[2]=0.0;
    func.coeff_poly[3]=-8.75;
    func.coeff_poly[4]=0.0;
    func.coeff_poly[5]=7.875;
    break;
  case 6:
    // Legendre polynomial coefficients of order six
    func.coeff_poly[0]=-0.3125;
    func.coeff_poly[1]=0.0;
    func.coeff_poly[2]=6.5625;
    func.coeff_poly[3]=0.0;
    func.coeff_poly[4]=-19.6875;
    func.coeff_poly[5]=0.0;
    func.coeff_poly[6]=14.4375;
    break;
  default:
    action->error("Insert Legendre polynomial coefficients into SphericalHarmonics code");
  }

  std::string num;
  for(int i=-func.tmom; i<=func.tmom; ++i) {
    Tools::convert(std::abs(i),num);
    if( i<0 ) {
      options.multipleValuesForEachRegisteredComponent.push_back( "-n" + num );
    } else if( i>0 ) {
      options.multipleValuesForEachRegisteredComponent.push_back( "-p" + num );
    } else {
      options.multipleValuesForEachRegisteredComponent.push_back( "-0");
    }
  }
  options.derivativeZeroIfValueIsZero = (action->getNumberOfArguments()==4
                                         && action->getPntrToArgument(3)
                                         ->isDerivativeZeroWhenValueIsZero() );
}

void SphericalHarmonic::calc( const SphericalHarmonic& func,
                              bool noderiv,
                              const View<const double> args,
                              function::FunctionOutput funcout ) {
  double weight=1;
  if( args.size()==4 ) {
    weight = args[3];
  }
  if( weight<epsilon ) {
    if( !noderiv ) {
      const unsigned imbase = 2*func.tmom+1;
      for(int m=-func.tmom; m<=func.tmom; ++m) {
        auto pos = func.tmom+m;
        funcout.derivs[pos][0] = 0.0;
        funcout.derivs[pos][1] = 0.0;
        funcout.derivs[pos][2] = 0.0;
        funcout.derivs[pos][3] = 0.0;
        funcout.derivs[pos + imbase][0] = 0.0;
        funcout.derivs[pos + imbase][1] = 0.0;
        funcout.derivs[pos + imbase][2] = 0.0;
        if( args.size()==4 ) {
          funcout.derivs[pos + imbase][3] = 0.0;
        }
      }
    }
    return;
  }

  const double dlen2 = args[0]*args[0] + args[1]*args[1] + args[2]*args[2];
  const double dlen  = sqrt( dlen2 );
  const double dleninv = 1.0/dlen;
  const double dlen3 = dlen2*dlen;
  const double dlen3inv = 1.0/dlen3;
  double dpoly_ass;
  double poly_ass=func.deriv_poly( 0, args[2]*dleninv, dpoly_ass );
  // Accumulate for m=0
  funcout.values[func.tmom] = weight*poly_ass;
  // Derivatives of z/r wrt x, y, z
  Vector dz;
  if( !noderiv ) {
    dz[0] = -( args[2] * dlen3inv )*args[0];
    dz[1] = -( args[2] * dlen3inv )*args[1];
    dz[2] = -( args[2] * dlen3inv )*args[2] + (1.0 * dleninv);
    funcout.derivs[func.tmom] = weight*dpoly_ass*dz ;
    if( args.size()==4 ) {
      funcout.derivs[func.tmom][3] = poly_ass;
    }
  }

  // The complex number of which we have to take powers
  std::complex<double> com1( args[0]*dleninv, args[1]*dleninv );
  std::complex<double> powered(1.0, 0.0);
  constexpr std::complex<double> ii(0.0, 1.0);
  std::complex<double> dp_x;
  std::complex<double> dp_y;
  std::complex<double> dp_z;
  Vector myrealvec;
  Vector myimagvec;
  Vector real_dz;
  Vector imag_dz;
  // Do stuff for all other m values
  const unsigned imbase = 3*func.tmom+1;
  double pref = -1.0;
  for(int m=1; m<=func.tmom; ++m) {
    // Calculate Legendre Polynomial
    poly_ass=func.deriv_poly( m, args[2]*dleninv, dpoly_ass );
    // Real and imaginary parts of z
    double real_z = real(com1*powered);
    double imag_z = imag(com1*powered);

    // Calculate steinhardt parameter
    double tq6 =poly_ass*real_z;  // Real part of steinhardt parameter
    double itq6=poly_ass*imag_z;  // Imaginary part of steinhardt parameter

    // Real part
    funcout.values[func.tmom + m] = weight* tq6;
    // Imaginary part
    funcout.values[imbase    + m] = weight*itq6;
    // Store -m part of vector
    // -m part of vector is just +m part multiplied by (-1.0)**m and multiplied by complex
    // conjugate of Legendre polynomial
    // Real part
    funcout.values[func.tmom - m] =  pref*funcout.values[func.tmom+m];
    // Imaginary part
    funcout.values[imbase    - m] = -pref*funcout.values[imbase + m];
    if( !noderiv ) {
      // Derivatives wrt ( x/r + iy )^m
      double md=static_cast<double>(m);
      dp_x = md*powered*( (1.0*dleninv)-(args[0]*args[0])*dlen3inv-ii*(args[0]*args[1])*dlen3inv );
      dp_y = md*powered*( ii*(1.0*dleninv)-(args[0]*args[1])*dlen3inv-ii*(args[1]*args[1])*dlen3inv );
      dp_z = md*powered*( -(args[0]*args[2])*dlen3inv-ii*(args[1]*args[2])*dlen3inv );
      // Derivatives of real and imaginary parts of above
      real_dz[0] = real( dp_x );
      real_dz[1] = real( dp_y );
      real_dz[2] = real( dp_z );
      imag_dz[0] = imag( dp_x );
      imag_dz[1] = imag( dp_y );
      imag_dz[2] = imag( dp_z );
      // Complete derivative of steinhardt parameter
      myrealvec = weight*(dpoly_ass*real_z*dz + poly_ass*real_dz);
      myimagvec = weight*(dpoly_ass*imag_z*dz + poly_ass*imag_dz);
      funcout.derivs[func.tmom+m] = myrealvec;
      funcout.derivs[imbase   +m] = myimagvec;

      funcout.derivs[func.tmom-m] = pref*myrealvec;
      funcout.derivs[imbase   -m] = -pref*myimagvec;
      if( args.size()==4 ) {
        funcout.derivs[func.tmom + m][3]=  tq6;
        funcout.derivs[imbase    + m][3]= itq6;
        funcout.derivs[func.tmom - m][3]= pref* tq6;
        funcout.derivs[imbase    - m][3]=-pref*itq6;
      }
    }
    // Calculate next power of complex number
    powered *= com1;
    //hopefully the compiler will optimize with bitflipping the sign here
    pref = -pref;
  }
}

double SphericalHarmonic::deriv_poly( const unsigned m,
                                      const double val,
                                      double& df ) const {
  double fact=1.0;
  for(unsigned j=2; j<=m; ++j) {
    fact *= j;
  }
  double res=coeff_poly[m]*fact;

  double pow=1.0;
  double xi=val;
  double dxi=1.0;
  df=0.0;
  for(int i=m+1; i<=tmom; ++i) {
    fact = 1.0;
    for(int j=i-m+1; j<=i; ++j) {
      fact *= j;
    }
    res += coeff_poly[i]*fact*xi;
    df  += pow*coeff_poly[i]*fact*dxi;
    xi  *= val;
    dxi *= val;
    pow += 1.0;
  }
  df *= normaliz[m];
  return normaliz[m]*res;
}

}
}

