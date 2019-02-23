/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2019 The plumed team
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
#include "Steinhardt.h"
#include "core/PlumedMain.h"
#include <complex>

namespace PLMD {
namespace crystallization {

void Steinhardt::registerKeywords( Keywords& keys ) {
  VectorMultiColvar::registerKeywords( keys );
  keys.add("compulsory","NN","12","The n parameter of the switching function ");
  keys.add("compulsory","MM","0","The m parameter of the switching function; 0 implies 2*NN");
  keys.add("compulsory","D_0","0.0","The d_0 parameter of the switching function");
  keys.add("compulsory","R_0","The r_0 parameter of the switching function");
  keys.add("optional","SWITCH","This keyword is used if you want to employ an alternative to the continuous switching function defined above. "
           "The following provides information on the \\ref switchingfunction that are available. "
           "When this keyword is present you no longer need the NN, MM, D_0 and R_0 keywords.");
  keys.use("SPECIES"); keys.use("SPECIESA"); keys.use("SPECIESB");
  keys.use("MEAN"); keys.use("LESS_THAN"); keys.use("MORE_THAN"); keys.use("VMEAN");
  keys.use("BETWEEN"); keys.use("HISTOGRAM"); keys.use("MOMENTS"); keys.use("MIN"); keys.use("ALT_MIN");
  keys.use("LOWEST"); keys.use("HIGHEST");
}

Steinhardt::Steinhardt( const ActionOptions& ao ):
  Action(ao),
  VectorMultiColvar(ao),
  tmom(0)
{
  // Read in the switching function
  std::string sw, errors; parse("SWITCH",sw);
  if(sw.length()>0) {
    switchingFunction.set(sw,errors);
  } else {
    double r_0=-1.0, d_0; int nn, mm;
    parse("NN",nn); parse("MM",mm);
    parse("R_0",r_0); parse("D_0",d_0);
    if( r_0<0.0 ) error("you must set a value for R_0");
    switchingFunction.set(nn,mm,r_0,d_0);
  }
  log.printf("  Steinhardt parameter of central atom and those within %s\n",( switchingFunction.description() ).c_str() );
  log<<"  Bibliography "<<plumed.cite("Tribello, Giberti, Sosso, Salvalaglio and Parrinello, J. Chem. Theory Comput. 13, 1317 (2017)")<<"\n";
  // Set the link cell cutoff
  setLinkCellCutoff( switchingFunction.get_dmax() );
  rcut = switchingFunction.get_dmax(); rcut2 = rcut*rcut;
  std::vector<AtomNumber> all_atoms; setupMultiColvarBase( all_atoms );
}

void Steinhardt::setAngularMomentum( const unsigned& ang ) {
  tmom=ang; setVectorDimensionality( 2*(2*ang + 1) );
}

void Steinhardt::calculateVector( multicolvar::AtomValuePack& myatoms ) const {
  double dfunc, dpoly_ass, md, tq6, itq6, real_z, imag_z;
  Vector dz, myrealvec, myimagvec, real_dz, imag_dz;
  // The square root of -1
  std::complex<double> ii( 0.0, 1.0 ), dp_x, dp_y, dp_z;

  unsigned ncomp=2*tmom+1;
  double sw, poly_ass, dlen; std::complex<double> powered;
  for(unsigned i=1; i<myatoms.getNumberOfAtoms(); ++i) {
    Vector& distance=myatoms.getPosition(i);  // getSeparation( myatoms.getPosition(0), myatoms.getPosition(i) );
    double d2;
    if ( (d2=distance[0]*distance[0])<rcut2 &&
         (d2+=distance[1]*distance[1])<rcut2 &&
         (d2+=distance[2]*distance[2])<rcut2 &&
         d2>epsilon ) {

      dlen = sqrt(d2);
      sw = switchingFunction.calculate( dlen, dfunc );
      accumulateSymmetryFunction( -1, i, sw, (+dfunc)*distance, (-dfunc)*Tensor( distance,distance ), myatoms );
      double dlen3 = d2*dlen;
      // Do stuff for m=0
      poly_ass=deriv_poly( 0, distance[2]/dlen, dpoly_ass );
      // Derivatives of z/r wrt x, y, z
      dz = -( distance[2] / dlen3 )*distance; dz[2] += (1.0 / dlen);
      // Derivative wrt to the vector connecting the two atoms
      myrealvec = (+sw)*dpoly_ass*dz + poly_ass*(+dfunc)*distance;
      // Accumulate the derivatives
      accumulateSymmetryFunction( 2 + tmom, i, sw*poly_ass, myrealvec, Tensor( -myrealvec,distance ), myatoms );

      // The complex number of which we have to take powers
      std::complex<double> com1( distance[0]/dlen,distance[1]/dlen );
      powered = std::complex<double>(1.0,0.0);

      // Do stuff for all other m values
      for(unsigned m=1; m<=tmom; ++m) {
        // Calculate Legendre Polynomial
        poly_ass=deriv_poly( m, distance[2]/dlen, dpoly_ass );
        // Calculate power of complex number
        // if(std::abs(com1)>epsilon) powered=pow(com1,m-1);
        // else if(m==1) powered=std::complex<double>(1.,0);
        // else powered = std::complex<double>(0.,0.);
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
        myrealvec = (+sw)*dpoly_ass*real_z*dz + (+dfunc)*distance*tq6 + (+sw)*poly_ass*real_dz;
        myimagvec = (+sw)*dpoly_ass*imag_z*dz + (+dfunc)*distance*itq6 + (+sw)*poly_ass*imag_dz;

        // Real part
        accumulateSymmetryFunction( 2 + tmom + m, i, sw*tq6, myrealvec, Tensor( -myrealvec,distance ), myatoms );
        // Imaginary part
        accumulateSymmetryFunction( 2+ncomp+tmom+m, i, sw*itq6, myimagvec, Tensor( -myimagvec,distance ), myatoms );
        // Store -m part of vector
        double pref=pow(-1.0,m);
        // -m part of vector is just +m part multiplied by (-1.0)**m and multiplied by complex
        // conjugate of Legendre polynomial
        // Real part
        accumulateSymmetryFunction( 2+tmom-m, i, pref*sw*tq6, pref*myrealvec, pref*Tensor( -myrealvec,distance ), myatoms );
        // Imaginary part
        accumulateSymmetryFunction( 2+ncomp+tmom-m, i, -pref*sw*itq6, -pref*myimagvec, pref*Tensor( myimagvec,distance ), myatoms );
        // Calculate next power of complex number
        powered *= com1;
      }
    }
  }

  // Normalize
  updateActiveAtoms( myatoms );
  for(unsigned i=0; i<getNumberOfComponentsInVector(); ++i) myatoms.getUnderlyingMultiValue().quotientRule( 2+i, 2+i );
}

double Steinhardt::deriv_poly( const unsigned& m, const double& val, double& df ) const {
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
