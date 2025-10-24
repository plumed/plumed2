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

This action works similarly to [DISTANCE_FROM_CONTOUR](DISTANCE_FROM_CONTOUR.md). Within this action a field is constructed that measures the density
of the system at each point in space using:

$$
p(x,y,x) = \sum_{i=1}^N K\left[\frac{x-x_i}{\sigma_x},\frac{y-y_i}{\sigma_y},\frac{z-z_i}{\sigma_z} \right]
$$

In this expression $\sigma_x, \sigma_y$ and $\sigma_z$ are bandwidth parameters and
$K$ is one of a Gaussian kernel function.  With that field in place we can define a Willard-Chandler
surface is defined a surface of constant density in the above field $p(x,y,z)$.
In other words, we can define a set of points, $(x',y',z')$, in the box which have:

$$
p(x',y',z') = \rho
$$

where $\rho$ is some target density. In [DISTANCE_FROM_CONTOUR](DISTANCE_FROM_CONTOUR.md) we assume that this set of points lie on a manifold that
has the same topology as one or multiple planes.  Here, by contrast, we assume that this set of points lie on a manifold that has the same topology
as a sphere.  This action then returns the distance between this spherical manifold and the position of a test particle.  This distance is measured
along a vector perpendicular to the manifold.

## Examples

The following input calculates a [CONTACT_MATRIX](CONTACT_MATRIX.md) between a set of atoms in which element $i,j$ is only non-zero if atoms $i$ and $j$
are within 6 nm of each other.  We then use this matrix as input for a [DFSCLUSTERING](DFSCLUSTERING.md) action that finds the largest connected component
in the matrix.  The [CENTER](CENTER.md) of this cluster is then identified and the location of an isocontour in:

$$
p(x,y,x) = \sum_{i=1}^N \xi_i f(c_i) K\left[\frac{x-x_i}{\sigma_x},\frac{y-y_i}{\sigma_y},\frac{z-z_i}{\sigma_z} \right]
$$

is found using this action.  In this expression $\xi_i$ is 1 if atom $i$ is part of the largest cluster and zero otherwise, $c_i$ is the coordination number of atom $i$ and
$f$ is a swtiching function.  The distance between this isocontour and position of atom 513 as well as the distance between `com` (the center of the largest cluster)
and the isocontour is then output to a file called colvar.

```plumed
ones: ONES SIZE=512
# Calculate contact matrix
c1_mat: CONTACT_MATRIX GROUP=1-512 SWITCH={EXP D_0=4.0 R_0=0.5 D_MAX=6.0}
# Calculate coordination numbers
c1: MATRIX_VECTOR_PRODUCT ARG=c1_mat,ones
# Select coordination numbers that are more than 2.0
cf: MORE_THAN ARG=c1 SWITCH={RATIONAL D_0=2.0 R_0=0.1}
# Find largest cluster
dfs: DFSCLUSTERING ARG=c1_mat
clust1: CLUSTER_WEIGHTS CLUSTERS=dfs CLUSTER=1
com: CENTER ATOMS=1-512 WEIGHTS=clust1 PHASES
# Filtered coordination numbers for atoms in largest cluster
ff: CUSTOM ARG=clust1,cf FUNC=x*y PERIODIC=NO
# Now do the multicolvar surface
dd: DISTANCE_FROM_SPHERICAL_CONTOUR ARG=ff POSITIONS=1-512 ATOM=513 ORIGIN=com BANDWIDTH=1.0,1.0,1.0 CONTOUR=0.5
PRINT ARG=dd.* FILE=colvar
```

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
  keys.addDOI("10.1021/acs.jpcb.8b03661");
}

DistanceFromSphericalContour::DistanceFromSphericalContour( const ActionOptions& ao ):
  Action(ao),
  DistanceFromContourBase(ao) {
  // Create the values
  std::vector<std::size_t> shape;
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
  pbcDistance( getPosition(getNumberOfAtoms()-1), getPosition(0) );
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
