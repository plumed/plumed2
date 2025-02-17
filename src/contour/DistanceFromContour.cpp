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

//+PLUMEDOC COLVAR DISTANCE_FROM_CONTOUR
/*
Calculate the perpendicular distance from a Willard-Chandler dividing surface.

Suppose that you have calculated a multicolvar.  By doing so you have calculated a
set of colvars, \f$s_i\f$, and each of these colvars has a well defined position in
space \f$(x_i,y_i,z_i)\f$.  You can use this information to calculate a phase-field
model of the colvar density using:

\f[
p(x,y,x) = \sum_{i} s_i K\left[\frac{x-x_i}{\sigma_x},\frac{y-y_i}{\sigma_y},\frac{z-z_i}{\sigma_z} \right]
\f]

In this expression \f$\sigma_x, \sigma_y\f$ and \f$\sigma_z\f$ are bandwidth parameters and
\f$K\f$ is one of the \ref kernelfunctions.  This is what is done within \ref MULTICOLVARDENS

The Willard-Chandler surface is a surface of constant density in the above phase field \f$p(x,y,z)\f$.
In other words, it is a set of points, \f$(x',y',z')\f$, in your box which have:

\f[
p(x',y',z') = \rho
\f]

where \f$\rho\f$ is some target density.  This action calculates the distance projected on the \f$x, y\f$ or
\f$z\f$ axis between the position of some test particle and this surface of constant field density.

\par Examples

In this example atoms 2-100 are assumed to be concentrated along some part of the \f$z\f$ axis so that you
an interface between a liquid/solid and the vapor.  The quantity dc measures the distance between the
surface at which the density of 2-100 atoms is equal to 0.2 and the position of the test particle atom 1.

\plumedfile
dens: DENSITY SPECIES=2-100
dc: DISTANCE_FROM_CONTOUR DATA=dens ATOM=1 BANDWIDTH=0.5,0.5,0.5 DIR=z CONTOUR=0.2
\endplumedfile

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace contour {

class DistanceFromContour : public DistanceFromContourBase {
private:
  unsigned dir;
  double pbc_param;
  std::vector<double> pos1, pos2, dirv, dirv2;
  std::vector<unsigned> perp_dirs;
  std::vector<Vector> atom_deriv;
public:
  static void registerKeywords( Keywords& keys );
  explicit DistanceFromContour( const ActionOptions& );
  void calculate() override;
  void evaluateDerivatives( const Vector& root1, const double& root2 );
};

PLUMED_REGISTER_ACTION(DistanceFromContour,"DISTANCE_FROM_CONTOUR")

void DistanceFromContour::registerKeywords( Keywords& keys ) {
  DistanceFromContourBase::registerKeywords( keys );
  keys.addOutputComponent("dist1","default","scalar","the distance between the reference atom and the nearest contour");
  keys.addOutputComponent("dist2","default","scalar","the distance between the reference atom and the other contour");
  keys.addOutputComponent("qdist","default","scalar","the differentiable (squared) distance between the two contours (see above)");
  keys.addOutputComponent("thickness","default","scalar","the distance between the two contours on the line from the reference atom");
  keys.add("compulsory","DIR","the direction perpendicular to the contour that you are looking for");
  keys.add("compulsory","TOLERANCE","0.1","this parameter is used to manage periodic boundary conditions.  The problem "
           "here is that we can be between contours even when we are not within the membrane "
           "because of periodic boundary conditions.  When we are in the contour, however, we "
           "should have it so that the sums of the absolute values of the distances to the two "
           "contours is approximately the distance between the two contours.  There can be numerical errors in these calculations, however, so "
           "we specify a small tolerance here");
}

DistanceFromContour::DistanceFromContour( const ActionOptions& ao ):
  Action(ao),
  DistanceFromContourBase(ao),
  pos1(3,0.0),
  pos2(3,0.0),
  dirv(3,0.0),
  dirv2(3,0.0),
  perp_dirs(2),
  atom_deriv(active_list.size()) {
  // Get the direction
  std::string ldir;
  parse("DIR",ldir );
  if( ldir=="x" ) {
    dir=0;
    perp_dirs[0]=1;
    perp_dirs[1]=2;
    dirv[0]=1;
    dirv2[0]=-1;
  } else if( ldir=="y" ) {
    dir=1;
    perp_dirs[0]=0;
    perp_dirs[1]=2;
    dirv[1]=1;
    dirv2[1]=-1;
  } else if( ldir=="z" ) {
    dir=2;
    perp_dirs[0]=0;
    perp_dirs[1]=1;
    dirv[2]=1;
    dirv2[2]=-1;
  } else {
    error(ldir + " is not a valid direction use x, y or z");
  }

  // Read in the tolerance for the pbc parameter
  parse("TOLERANCE",pbc_param);

  std::vector<unsigned> shape;
  // Create the values
  addComponent("thickness", shape );
  componentIsNotPeriodic("thickness");
  addComponent("dist1", shape );
  componentIsNotPeriodic("dist1");
  addComponent("dist2", shape );
  componentIsNotPeriodic("dist2");
  addComponentWithDerivatives("qdist", shape );
  componentIsNotPeriodic("qdist");
}

void DistanceFromContour::calculate() {
  // Check box is orthorhombic
  if( !getPbc().isOrthorombic() ) {
    error("cell box must be orthorhombic");
  }

  // The nanoparticle is at the origin of our coordinate system
  pos1[0]=pos1[1]=pos1[2]=0.0;
  pos2[0]=pos2[1]=pos2[2]=0.0;

  // Set bracket as center of mass of membrane in active region
  Vector myvec = pbcDistance( getPosition(getNumberOfAtoms()-1), getPosition(0) );
  pos2[dir]=myvec[dir];
  nactive=1;
  active_list[0]=0;
  double d2, mindist = myvec.modulo2();
  for(unsigned j=1; j<getNumberOfAtoms()-1; ++j) {
    Vector distance=pbcDistance( getPosition(getNumberOfAtoms()-1), getPosition(j) );
    if( (d2=distance[perp_dirs[0]]*distance[perp_dirs[0]])<rcut2 &&
        (d2+=distance[perp_dirs[1]]*distance[perp_dirs[1]])<rcut2 ) {
      d2+=distance[dir]*distance[dir];
      if( d2<mindist && fabs(distance[dir])>epsilon ) {
        pos2[dir]=distance[dir];
        mindist = d2;
      }
      active_list[nactive]=j;
      nactive++;
    }
  }
  // pos1 position of the nanoparticle, in the first time
  // pos2 is the position of the closer atom in the membrane with respect the nanoparticle
  // fa = distance between pos1 and the contour
  // fb = distance between pos2 and the contour
  std::vector<double> faked(3);
  double fa = getDifferenceFromContour( pos1, faked );
  double fb = getDifferenceFromContour( pos2, faked );
  if( fa*fb>0 ) {
    unsigned maxtries = std::floor( ( getBox()(dir,dir) ) / bw[dir] );
    for(unsigned i=0; i<maxtries; ++i) {
      double sign=(pos2[dir]>0)? -1 : +1; // If the nanoparticle is inside the membrane push it out
      pos1[dir] += sign*bw[dir];
      fa = getDifferenceFromContour( pos1, faked );
      if( fa*fb<0 ) {
        break;
      }
      // if fa*fb is less than zero the new pos 1 is outside the contour
    }
  }
  // Set direction for contour search
  dirv[dir] = pos2[dir] - pos1[dir];
  // Bracket for second root starts in center of membrane
  double fc = getDifferenceFromContour( pos2, faked );
  if( fc*fb>0 ) {
    // first time is true, because fc=fb
    // push pos2 from its initial position inside the membrane towards the second contourn
    unsigned maxtries = std::floor( ( getBox()(dir,dir) ) / bw[dir] );
    for(unsigned i=0; i<maxtries; ++i) {
      double sign=(dirv[dir]>0)? +1 : -1;
      pos2[dir] += sign*bw[dir];
      fc = getDifferenceFromContour( pos2, faked );
      if( fc*fb<0 ) {
        break;
      }
    }
    dirv2[dir] = ( pos1[dir] + dirv[dir] ) - pos2[dir];
  }

  // Now do a search for the two contours
  findContour( dirv, pos1 );
  // Save the first value
  Vector root1;
  root1.zero();
  root1[dir] = pval[dir];
  findContour( dirv2, pos2 );
  // Calculate the separation between the two roots using PBC
  Vector root2;
  root2.zero();
  root2[dir] = pval[dir];
  Vector sep = pbcDistance( root1, root2 );
  double spacing = fabs( sep[dir] );
  plumed_assert( spacing>epsilon );
  getPntrToComponent("thickness")->set( spacing );

  // Make sure the sign is right
  double predir=(root1[dir]*root2[dir]<0)? -1 : 1;
  // This deals with periodic boundary conditions - if we are inside the membrane the sum of the absolute
  // distances from the contours should add up to the spacing.  When this is not the case we must be outside
  // the contour
  // if( predir==-1 && (fabs(root1[dir])+fabs(root2[dir]))>(spacing+pbc_param) ) predir=1;
  // Set the final value to root that is closest to the "origin" = position of atom
  if( fabs(root1[dir])<fabs(root2[dir]) ) {
    getPntrToComponent("dist1")->set( predir*fabs(root1[dir]) );
    getPntrToComponent("dist2")->set( fabs(root2[dir]) );
  } else {
    getPntrToComponent("dist1")->set( predir*fabs(root2[dir]) );
    getPntrToComponent("dist2")->set( fabs(root1[dir]) );
  }
  getPntrToComponent("qdist")->set( root2[dir]*root1[dir] );

  // Now calculate the derivatives
  if( !doNotCalculateDerivatives() ) {
    evaluateDerivatives( root1, root2[dir] );
    evaluateDerivatives( root2, root1[dir] );
  }
}

void DistanceFromContour::evaluateDerivatives( const Vector& root1, const double& root2 ) {
  if( getNumberOfArguments()>0 ) {
    plumed_merror("derivatives for phase field distance from contour have not been implemented yet");
  }

  Vector origind;
  origind.zero();
  Tensor vir;
  vir.zero();
  double sumd = 0;
  std::vector<double> pp(3), ddd(3,0);
  for(unsigned i=0; i<nactive; ++i) {
    double newval = evaluateKernel( getPosition(active_list[i]), root1, ddd );
    Vector distance = pbcDistance( getPosition(getNumberOfAtoms()-1), getPosition(active_list[i]) );

    if( getNumberOfArguments()==1 ) {
    } else {
      sumd += ddd[dir];
      for(unsigned j=0; j<3; ++j) {
        atom_deriv[i][j] = -ddd[j];
      }
      origind += -atom_deriv[i];
      vir -= Tensor(atom_deriv[i],distance);
    }
  }

  // Add derivatives to atoms involved
  Value* val=getPntrToComponent("qdist");
  double prefactor =  root2 / sumd;
  for(unsigned i=0; i<nactive; ++i) {
    val->addDerivative( 3*active_list[i] + 0, -prefactor*atom_deriv[i][0] );
    val->addDerivative( 3*active_list[i] + 1, -prefactor*atom_deriv[i][1] );
    val->addDerivative( 3*active_list[i] + 2, -prefactor*atom_deriv[i][2] );
  }

  // Add derivatives to atoms at origin
  unsigned nbase = 3*(getNumberOfAtoms()-1);
  val->addDerivative( nbase, -prefactor*origind[0] );
  nbase++;
  val->addDerivative( nbase, -prefactor*origind[1] );
  nbase++;
  val->addDerivative( nbase, -prefactor*origind[2] );
  nbase++;

  // Add derivatives to virial
  for(unsigned i=0; i<3; ++i)
    for(unsigned j=0; j<3; ++j) {
      val->addDerivative( nbase, -prefactor*vir(i,j) );
      nbase++;
    }
}

}
}
