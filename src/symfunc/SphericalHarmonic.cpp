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

//+PLUMEDOC MCOLVAR SPHERICAL_HARMONIC 
/*

\par Examples


*/
//+ENDPLUMEDOC


class SphericalHarmonic : public SymmetryFunctionBase {
private:
  int tmom;
  std::vector<double> coeff_poly;
  std::vector<double> normaliz;
  unsigned factorial( const unsigned& n ) const ;
  double deriv_poly( const unsigned& m, const double& val, double& df ) const ;
public:
  static void shortcutKeywords( Keywords& keys );
  static void expandShortcut( const std::string& lab, const std::vector<std::string>& words,
                              const std::map<std::string,std::string>& keys,
                              std::vector<std::vector<std::string> >& actions );
  static void registerKeywords( Keywords& keys );
  explicit SphericalHarmonic(const ActionOptions&);
  void compute( const double& val, const Vector& dir, MultiValue& myvals ) const ;
};

PLUMED_REGISTER_ACTION(SphericalHarmonic,"SPHERICAL_HARMONIC")
PLUMED_REGISTER_SHORTCUT(SphericalHarmonic,"SPHERICAL_HARMONIC")
PLUMED_REGISTER_SHORTCUT(SphericalHarmonic,"Q3")
PLUMED_REGISTER_SHORTCUT(SphericalHarmonic,"Q4")
PLUMED_REGISTER_SHORTCUT(SphericalHarmonic,"Q6")

void SphericalHarmonic::shortcutKeywords( Keywords& keys ) {
  SymmetryFunctionBase::shortcutKeywords( keys );
}

void SphericalHarmonic::expandShortcut( const std::string& lab, const std::vector<std::string>& words,
                                  const std::map<std::string,std::string>& keys,
                                  std::vector<std::vector<std::string> >& actions ) {
  SymmetryFunctionBase::expandMatrix( true, lab, words, keys, actions );
  std::vector<std::string> sph_input; sph_input.push_back( lab + ":" ); 
  sph_input.push_back("SPHERICAL_HARMONIC"); sph_input.push_back("WEIGHT=" + lab + "_mat.w");
  sph_input.push_back("VECTORS1=" + lab + "_mat.x" ); sph_input.push_back("VECTORS2=" + lab + "_mat.y" ); 
  sph_input.push_back("VECTORS3=" + lab + "_mat.z" ); int l; 
  if( words[0]=="Q3" ) { 
      sph_input.push_back("L=3"); l=3;
  } else if( words[0]=="Q4" ) {
      sph_input.push_back("L=4"); l=4;
  } else if( words[0]=="Q6" ) {
      sph_input.push_back("L=6"); l=6;
  } else {
      plumed_merror("invalid input");
  }
  actions.push_back( sph_input );

  // Input for denominator (coord)
  std::vector<std::string> d_input; d_input.push_back(lab + "_denom:"); d_input.push_back("COORDINATIONNUMBER");
  d_input.push_back("WEIGHT=" + lab + "_mat.w"); actions.push_back( d_input );

  // Divide all components by coordination numbers
  for(int i=-l;i<=l;++i){
      std::string snum; Tools::convert( i, snum ); 
      // Real part 
      std::vector<std::string> rm_input; rm_input.push_back(lab + "_rmn-[" + snum + "]:"); 
      rm_input.push_back("MATHEVAL"); rm_input.push_back("ARG1=" + lab + ".rm-[" + snum + "]"); 
      rm_input.push_back("ARG2=" + lab + "_denom"); rm_input.push_back("FUNC=x/y"); rm_input.push_back("PERIODIC=NO");
      actions.push_back(rm_input);
      // Imaginary part
      std::vector<std::string> im_input; im_input.push_back(lab + "_imn-[" + snum + "]:");
      im_input.push_back("MATHEVAL"); im_input.push_back("ARG1=" + lab + ".im-[" + snum + "]");
      im_input.push_back("ARG2=" + lab + "_denom"); im_input.push_back("FUNC=x/y"); im_input.push_back("PERIODIC=NO");
      actions.push_back(im_input);
  }

  // Now calculate the total length of the vector
  std::vector<std::string> norm_input; 
  norm_input.push_back(lab +"_norm2:"); 
  norm_input.push_back("COMBINE"); std::string powstr="POWERS=2";
  std::string snum, num; unsigned nn=1;
  for(int i=-l;i<=l;++i){ 
     Tools::convert( nn, num ); Tools::convert( i, snum ); 
     norm_input.push_back("ARG" + num + "=" + lab + "_rmn-[" + snum + "]"); 
     nn++; Tools::convert( nn, num );
     norm_input.push_back("ARG" + num + "=" + lab + "_imn-[" + snum + "]"); 
     nn++; 
     if( i==-l ) powstr += ",2"; else powstr += ",2,2";
  } 
  norm_input.push_back("PERIODIC=NO"); norm_input.push_back(powstr); actions.push_back( norm_input );
  std::vector<std::string> sqrt_input; sqrt_input.push_back(lab +"_norm:"); sqrt_input.push_back("MATHEVAL");
  sqrt_input.push_back("ARG1=" + lab + "_norm2"); sqrt_input.push_back( "FUNC=sqrt(x)" ); 
  sqrt_input.push_back("PERIODIC=NO"); actions.push_back( sqrt_input );
  multicolvar::MultiColvarBase::expandFunctions( lab, lab + "_norm", words, keys, actions );
}

void SphericalHarmonic::registerKeywords( Keywords& keys ) {
  SymmetryFunctionBase::registerKeywords( keys );
  keys.add("compulsory","L","the value of the angular momentum");
  keys.addOutputComponent("rm","default","the real parts of the spherical harmonic values with the m value given");
  keys.addOutputComponent("im","default","the real parts of the spherical harmonic values with the m value given");
}

unsigned SphericalHarmonic::factorial( const unsigned& n ) const {
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

SphericalHarmonic::SphericalHarmonic(const ActionOptions&ao):
  Action(ao),
  SymmetryFunctionBase(ao)
{
  parse("L",tmom);
  log.printf("  calculating %dth order spherical harmonics \n", tmom);  
  for(int i=-tmom;i<=tmom;++i){
      std::string num; Tools::convert(i,num);
      addComponentWithDerivatives( "rm-[" + num + "]" ); 
  }
  for(int i=-tmom;i<=tmom;++i){
      std::string num; Tools::convert(i,num);
      addComponentWithDerivatives( "im-[" + num + "]" );
  }
  normaliz.resize( tmom+1 ); 
  for(unsigned i=0;i<=tmom;++i){
      normaliz[i] = sqrt( (2*tmom+1)*factorial(tmom-i)/(4*pi*factorial(tmom+i)) );
      if( i%2==1 ) normaliz[i]*=-1; 
  }

  coeff_poly.resize( tmom+1 );
  if( tmom==3 ) {
      // Legendre polynomial coefficients of order three
      coeff_poly[0]=0.0; coeff_poly[1]=-1.5;
      coeff_poly[2]=0.0; coeff_poly[3]=2.5;
  } else if( tmom==4 ) {
      // Legendre polynomial coefficients of order four
      coeff_poly[0]=0.375; coeff_poly[1]=0.0;
      coeff_poly[2]=-3.75; coeff_poly[3]=0.0;
      coeff_poly[4]=4.375;
  } else if( tmom==6 ) {
      // Legendre polynomial coefficients of order six
      coeff_poly[0]=-0.3125; coeff_poly[1]=0.0;
      coeff_poly[2]=6.5625; coeff_poly[3]=0.0;
      coeff_poly[4]=-19.6875; coeff_poly[5]=0.0;
      coeff_poly[6]=14.4375;
  } else {
      error("Insert Legendre polynomial coefficients into SphericalHarmonics code");
  }
  checkRead();
}

void SphericalHarmonic::compute( const double& val, const Vector& distance, MultiValue& myvals ) const {
  double dlen2 = distance.modulo2(); double dlen = sqrt( dlen2 ); double dlen3 = dlen2*dlen;
  double tq6, itq6, dpoly_ass, poly_ass=deriv_poly( 0, distance[2]/dlen, dpoly_ass );  
  // Derivatives of z/r wrt x, y, z
  Vector dz = -( distance[2] / dlen3 )*distance; dz[2] += (1.0 / dlen);
  // Accumulate for m=0
  addToValue( tmom, val*poly_ass, myvals ); 
  addVectorDerivatives( tmom, val*dpoly_ass*dz, myvals ); 
  addWeightDerivative( tmom, poly_ass, myvals );  

  // The complex number of which we have to take powers
  std::complex<double> com1( distance[0]/dlen,distance[1]/dlen ), dp_x, dp_y, dp_z; double md, real_z, imag_z;
  std::complex<double> powered = std::complex<double>(1.0,0.0); std::complex<double> ii( 0.0, 1.0 ); 
  Vector myrealvec, myimagvec, real_dz, imag_dz;
  // Do stuff for all other m values
  for(unsigned m=1; m<=tmom; ++m) {
    // Calculate Legendre Polynomial
    poly_ass=deriv_poly( m, distance[2]/dlen, dpoly_ass );
    // Real and imaginary parts of z
    real_z = real(com1*powered); imag_z = imag(com1*powered );
    
    // Calculate steinhardt parameter
    tq6=poly_ass*real_z;   // Real part of steinhardt parameter
    itq6=poly_ass*imag_z;  // Imaginary part of steinhardt parameter
    
    // Derivatives wrt ( x/r + iy )^m
    md=static_cast<double>(m);
    dp_x = md*powered*( (1.0/dlen)-(distance[0]*distance[0])/dlen3-ii*(distance[0]*distance[1])/dlen3 );
    dp_y = md*powered*( ii*(1.0/dlen)-(distance[0]*distance[1])/dlen3-ii*(distance[1]*distance[1])/dlen3 );
    dp_z = md*powered*( -(distance[0]*distance[2])/dlen3-ii*(distance[1]*distance[2])/dlen3 );
    
    // Derivatives of real and imaginary parts of above
    real_dz[0] = real( dp_x ); real_dz[1] = real( dp_y ); real_dz[2] = real( dp_z );
    imag_dz[0] = imag( dp_x ); imag_dz[1] = imag( dp_y ); imag_dz[2] = imag( dp_z );

    // Complete derivative of steinhardt parameter
    myrealvec = val*dpoly_ass*real_z*dz + val*poly_ass*real_dz;
    myimagvec = val*dpoly_ass*imag_z*dz + val*poly_ass*imag_dz;

    // Real part
    addToValue( tmom+m, val*tq6, myvals ); 
    addVectorDerivatives( tmom+m, myrealvec, myvals ); 
    addWeightDerivative( tmom+m, tq6, myvals );
    // Imaginary part
    addToValue( 3*tmom+1+m, val*itq6, myvals ); 
    addVectorDerivatives( 3*tmom+1+m, myimagvec, myvals ); 
    addWeightDerivative( 3*tmom+1+m, itq6, myvals ); 
    // Store -m part of vector
    double pref=pow(-1.0,m);
    // -m part of vector is just +m part multiplied by (-1.0)**m and multiplied by complex
    // conjugate of Legendre polynomial
    // Real part
    addToValue( tmom-m, pref*val*tq6, myvals ); 
    addVectorDerivatives( tmom-m, pref*myrealvec, myvals ); 
    addWeightDerivative( tmom-m, pref*tq6, myvals ); 
    // Imaginary part
    addToValue( 3*tmom+1-m, -pref*val*itq6, myvals ); 
    addVectorDerivatives( 3*tmom+1-m, -pref*myimagvec, myvals ); 
    addWeightDerivative( 3*tmom+1-m, -pref*itq6, myvals );
    // Calculate next power of complex number
    powered *= com1;
  }
}

double SphericalHarmonic::deriv_poly( const unsigned& m, const double& val, double& df ) const {
  double fact=1.0;
  for(unsigned j=1; j<=m; ++j) fact=fact*j;
  double res=coeff_poly[m]*fact;

  double pow=1.0, xi=val, dxi=1.0; df=0.0;
  for(int i=m+1; i<=tmom; ++i) {
    double fact=1.0;
    for(unsigned j=i-m+1; j<=i; ++j) fact=fact*j;
    res=res+coeff_poly[i]*fact*xi;
    df = df + pow*coeff_poly[i]*fact*dxi;
    xi=xi*val; dxi=dxi*val; pow+=1.0;
  }
  df = df*normaliz[m];
  return normaliz[m]*res;
}

}
}

