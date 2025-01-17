/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2023 The plumed team
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
#include "DistanceFromContourBase.h"
#include "core/ActionRegister.h"

//+PLUMEDOC COLVAR DISTANCE_FROM_SPHERICAL_CONTOUR
/*
Calculate the perpendicular distance from a Willard-Chandler dividing surface.

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace contour {

class DistanceFromSphericalContour : public DistanceFromContourBase {
public:
  static void registerKeywords( Keywords& keys );
  explicit DistanceFromSphericalContour( const ActionOptions& );
  void calculate();
  void evaluateDerivatives( const Vector& root1, const double& root2 );
};

PLUMED_REGISTER_ACTION(DistanceFromSphericalContour,"DISTANCE_FROM_SPHERICAL_CONTOUR")

void DistanceFromSphericalContour::registerKeywords( Keywords& keys ) {
  DistanceFromContourBase::registerKeywords( keys );
  keys.addOutputComponent("dist","default","scalar","the distance between the reference atom and the nearest contour");
  keys.addOutputComponent("radius","default","scalar","the radial distance from the center of the contour to the edge");
  keys.add("atoms","ORIGIN","The position of the center of the region that the contour encloses");
}

DistanceFromSphericalContour::DistanceFromSphericalContour( const ActionOptions& ao ):
  Action(ao),
  DistanceFromContourBase(ao) {
  // Create the values
  std::vector<unsigned> shape;
  addComponentWithDerivatives("dist", shape );
  componentIsNotPeriodic("dist");
  addComponent("radius", shape );
  componentIsNotPeriodic("radius");
}

void DistanceFromSphericalContour::calculate() {
  // Check box is orthorhombic
  if( !getPbc().isOrthorombic() ) {
    error("cell box must be orthorhombic");
  }

  // Calculate the director of the vector connecting the center of the sphere to the molecule of interest
  Vector dirv = pbcDistance( getPosition(getNumberOfAtoms()-1), getPosition(getNumberOfAtoms()-2) );
  double len=dirv.modulo();
  dirv /= len;
  // Now work out which atoms need to be considered explicitly
  Vector myvec = pbcDistance( getPosition(getNumberOfAtoms()-1), getPosition(0) );
  nactive=1;
  active_list[0]=0;
  for(unsigned j=1; j<getNumberOfAtoms()-2; ++j) {
    if( getNumberOfArguments()==1 ) {
      if( getPntrToArgument(0)->get(j)<epsilon ) {
        continue;
      }
    }
    active_list[nactive]=j;
    nactive++;
    Vector distance=pbcDistance( getPosition(getNumberOfAtoms()-1), getPosition(j) );
    double dp = dotProduct( distance, dirv );
    double cp = distance.modulo2() - dp*dp;
    if( cp<rcut2 ) {
      active_list[nactive]=j;
      nactive++;
    }
  }
  // Get maximum length to fit in box
  double hbox = 0.5*getBox()(0,0);
  if( 0.5*getBox()(1,1)<hbox ) {
    hbox = 0.5*getBox()(1,1);
  }
  if( 0.5*getBox()(2,2)<hbox ) {
    hbox = 0.5*getBox()(2,2);
  }
  // Set initial guess for position of contour to position of closest molecule in region
  std::vector<double> pos1(3), dirv2(3);
  for(unsigned k=0; k<3; ++k) {
    dirv2[k]=hbox*dirv[k];
    pos1[k]=0;
  }
  // Now do a search for the contours
  findContour( dirv2, pos1 );
  // Now find the distance between the center of the sphere and the contour
  double rad = sqrt( pos1[0]*pos1[0] + pos1[1]*pos1[1] + pos1[2]*pos1[2] );
  // Set the radius
  getPntrToComponent("radius")->set( rad );
  // Set the distance between the contour and the molecule
  getPntrToComponent("dist")->set( len - rad );

  // Now calculate the derivatives
  if( !doNotCalculateDerivatives() ) {
    plumed_merror("derivatives not implemented");
  }
}

void DistanceFromSphericalContour::evaluateDerivatives( const Vector& root1, const double& root2 ) {
  plumed_error();
}

}
}
