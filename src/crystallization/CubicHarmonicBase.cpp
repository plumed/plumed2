/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2019 The plumed team
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
#include "CubicHarmonicBase.h"
#include "tools/SwitchingFunction.h"

#include <string>
#include <cmath>

using namespace std;

namespace PLMD {
namespace crystallization {

void CubicHarmonicBase::registerKeywords( Keywords& keys ) {
  multicolvar::MultiColvarBase::registerKeywords( keys );
  keys.use("SPECIES"); keys.use("SPECIESA"); keys.use("SPECIESB");
  keys.add("compulsory","NN","6","The n parameter of the switching function ");
  keys.add("compulsory","MM","0","The m parameter of the switching function; 0 implies 2*NN");
  keys.add("compulsory","D_0","0.0","The d_0 parameter of the switching function");
  keys.add("compulsory","R_0","The r_0 parameter of the switching function");
  keys.add("optional","SWITCH","This keyword is used if you want to employ an alternative to the continuous switching function defined above. "
           "The following provides information on the \\ref switchingfunction that are available. "
           "When this keyword is present you no longer need the NN, MM, D_0 and R_0 keywords.");
  keys.add("compulsory","PHI","0.0","The Euler rotational angle phi");
  keys.add("compulsory","THETA","0.0","The Euler rotational angle theta");
  keys.add("compulsory","PSI","0.0","The Euler rotational angle psi");
  keys.addFlag("UNORMALIZED",false,"calculate the sum of the components of the vector rather than the mean");
  // Use actionWithDistributionKeywords
  keys.use("MEAN"); keys.use("MORE_THAN"); keys.use("LESS_THAN"); keys.use("MAX");
  keys.use("MIN"); keys.use("BETWEEN"); keys.use("HISTOGRAM"); keys.use("MOMENTS");
  keys.use("ALT_MIN"); keys.use("LOWEST"); keys.use("HIGHEST");
}

CubicHarmonicBase::CubicHarmonicBase(const ActionOptions&ao):
  Action(ao),
  MultiColvarBase(ao)
{
  // Read in the switching function
  std::string sw, errors; parse("SWITCH",sw);
  if(sw.length()>0) {
    switchingFunction.set(sw,errors);
    if( errors.length()!=0 ) error("problem reading SWITCH keyword : " + errors );
  } else {
    double r_0=-1.0, d_0; int nn, mm;
    parse("NN",nn); parse("MM",mm);
    parse("R_0",r_0); parse("D_0",d_0);
    if( r_0<0.0 ) error("you must set a value for R_0");
    switchingFunction.set(nn,mm,r_0,d_0);
  }

  double phi, theta, psi; parse("PHI",phi); parse("THETA",theta); parse("PSI",psi);
  log.printf("  creating rotation matrix with Euler angles phi=%f, theta=%f and psi=%f\n",phi,theta,psi);
  // Calculate the rotation matrix http://mathworld.wolfram.com/EulerAngles.html
  rotationmatrix[0][0]=cos(psi)*cos(phi)-cos(theta)*sin(phi)*sin(psi);
  rotationmatrix[0][1]=cos(psi)*sin(phi)+cos(theta)*cos(phi)*sin(psi);
  rotationmatrix[0][2]=sin(psi)*sin(theta);

  rotationmatrix[1][0]=-sin(psi)*cos(phi)-cos(theta)*sin(phi)*cos(psi);
  rotationmatrix[1][1]=-sin(psi)*sin(phi)+cos(theta)*cos(phi)*cos(psi);
  rotationmatrix[1][2]=cos(psi)*sin(theta);

  rotationmatrix[2][0]=sin(theta)*sin(phi);
  rotationmatrix[2][1]=-sin(theta)*cos(phi);
  rotationmatrix[2][2]=cos(theta);


  log.printf("  measure crystallinity around central atom.  Includes those atoms within %s\n",( switchingFunction.description() ).c_str() );
  parseFlag("UNORMALIZED",unormalized);
  if( unormalized ) log.printf("  output sum of vector functions \n");
  else log.printf("  output mean of vector functions \n");
  // Set the link cell cutoff
  rcut2 = switchingFunction.get_dmax()*switchingFunction.get_dmax();
  setLinkCellCutoff( switchingFunction.get_dmax() );
  // And setup the ActionWithVessel
  std::vector<AtomNumber> all_atoms; setupMultiColvarBase( all_atoms );
}

double CubicHarmonicBase::compute( const unsigned& tindex, multicolvar::AtomValuePack& myatoms ) const {
  double dfunc; Vector rotatedis;

  // Calculate the coordination number
  Vector myder, rotateder, fder; unsigned nat=myatoms.getNumberOfAtoms();

  for(unsigned i=1; i<nat; ++i) {
    Vector& distance=myatoms.getPosition(i);

    double d2;
    if ( (d2=distance[0]*distance[0])<rcut2 &&
         (d2+=distance[1]*distance[1])<rcut2 &&
         (d2+=distance[2]*distance[2])<rcut2 &&
         d2>epsilon ) {

      double sw = switchingFunction.calculateSqr( d2, dfunc );

      rotatedis[0]=rotationmatrix[0][0]*distance[0]
                   +rotationmatrix[0][1]*distance[1]
                   +rotationmatrix[0][2]*distance[2];
      rotatedis[1]=rotationmatrix[1][0]*distance[0]
                   +rotationmatrix[1][1]*distance[1]
                   +rotationmatrix[1][2]*distance[2];
      rotatedis[2]=rotationmatrix[2][0]*distance[0]
                   +rotationmatrix[2][1]*distance[1]
                   +rotationmatrix[2][2]*distance[2];

      double tmp = calculateCubicHarmonic( rotatedis, d2, rotateder );

      myder[0]=rotationmatrix[0][0]*rotateder[0]
               +rotationmatrix[1][0]*rotateder[1]
               +rotationmatrix[2][0]*rotateder[2];
      myder[1]=rotationmatrix[0][1]*rotateder[0]
               +rotationmatrix[1][1]*rotateder[1]
               +rotationmatrix[2][1]*rotateder[2];
      myder[2]=rotationmatrix[0][2]*rotateder[0]
               +rotationmatrix[1][2]*rotateder[1]
               +rotationmatrix[2][2]*rotateder[2];

      fder = (+dfunc)*tmp*distance + sw*myder;

      accumulateSymmetryFunction( 1, i, sw*tmp, fder, Tensor(distance,-fder), myatoms );
      accumulateSymmetryFunction( -1, i, sw, (+dfunc)*distance, (-dfunc)*Tensor(distance,distance), myatoms );
    }
  }
  // values -> der of... value [0], weight[1], x coord [2], y, z... [more magic]
  updateActiveAtoms( myatoms );
  if( !unormalized ) myatoms.getUnderlyingMultiValue().quotientRule( 1, 1 );
  return myatoms.getValue(1); // this is equivalent to getting an "atomic" CV
}

}
}

