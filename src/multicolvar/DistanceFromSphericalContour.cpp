/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2018 The plumed team
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
#include "tools/KernelFunctions.h"
#include "tools/RootFindingBase.h"

//+PLUMEDOC COLVAR DISTANCE_FROM_SPHERICAL_CONTOUR
/*
Calculate the perpendicular distance from a Willard-Chandler dividing surface.

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace multicolvar {

class DistanceFromSphericalContour : public DistanceFromContourBase {
private:
  RootFindingBase<DistanceFromSphericalContour> mymin;
public:
  static void registerKeywords( Keywords& keys );
  explicit DistanceFromSphericalContour( const ActionOptions& );
  void calculate();
  void evaluateDerivatives( const Vector root1, const double& root2 );
};

PLUMED_REGISTER_ACTION(DistanceFromSphericalContour,"DISTANCE_FROM_SPHERICAL_CONTOUR")

void DistanceFromSphericalContour::registerKeywords( Keywords& keys ) {
  DistanceFromContourBase::registerKeywords( keys ); 
  keys.addOutputComponent("dist","default","the distance between the reference atom and the nearest contour");
  keys.addOutputComponent("radius","default","the radial distance from the center of the contour to the edge");
  keys.add("atoms","ORIGIN","The position of the center of the region that the contour encloses");
}

DistanceFromSphericalContour::DistanceFromSphericalContour( const ActionOptions& ao ):
  Action(ao),
  DistanceFromContourBase(ao),
  mymin(this)
{
  // Create the values
  std::vector<unsigned> shape;
  addComponentWithDerivatives("dist", shape ); componentIsNotPeriodic("dist");
  addComponent("radius", shape ); componentIsNotPeriodic("radius");
}

void DistanceFromSphericalContour::calculate() {
  // Check box is orthorhombic
  if( !getPbc().isOrthorombic() ) error("cell box must be orthorhombic");

  // Calculate the director of the vector connecting particle to center of region
  Vector dirv = pbcDistance( getPosition(getNumberOfAtoms()-2), getPosition(getNumberOfAtoms()-1) );
  double len=dirv.modulo(); dirv /= len;
  // Now work out which atoms need to be considered explicitly
  Vector myvec = pbcDistance( getPosition(getNumberOfAtoms()-2), getPosition(0) ); 
  nactive=1; active_list[0]=0; 
  double dp = dotProduct( myvec, dirv ), mindist = fabs(dp);
  for(unsigned j=1; j<getNumberOfAtoms()-1; ++j) {
    Vector distance=pbcDistance( getPosition(getNumberOfAtoms()-2), getPosition(j) );
    double dp = dotProduct( distance, dirv ); double cp = distance.modulo2() - dp*dp;
    if( cp<rcut2 ) {
      if( fabs(dp)<mindist && fabs(dp)>epsilon ) { mindist = dp; }
      active_list[nactive]=j; nactive++;
    }
  }
  // Set initial guess for position of contour to position of closest molecule in region
  std::vector<double> pos1(3), dirv2(3); 
  for(unsigned k=0;k<3;++k){ dirv2[k]=dirv[k]; pos1[k]=mindist*dirv[k]; }
  // Now do a search for the contours
  mymin.lsearch( dirv2, pos1, &DistanceFromSphericalContour::getDifferenceFromContour );
  // Now get the final position
  Vector pos2; for(unsigned i=0;i<3;++i) pos2[i] = pval[i]->get();
  // And the dot product is the value of the CV
  double fval = dotProduct( pos2, dirv );
  // Set the radius 
  getPntrToComponent("radius")->set( len - fval );
  getPntrToComponent("dist")->set( fval );

  // Now calculate the derivatives
  // if( !doNotCalculateDerivatives() ) {
  //   evaluateDerivatives( root1, root2[dir] ); evaluateDerivatives( root2, root1[dir] );
  // }
}

void DistanceFromSphericalContour::evaluateDerivatives( const Vector root1, const double& root2 ) {
//   if( getNumberOfArguments()>0 ) plumed_merror("derivatives for phase field distance from contour have not been implemented yet");
//   for(unsigned j=0; j<3; ++j) pval[j]->set( root1[j] );
// 
//   Vector origind; origind.zero(); Tensor vir; vir.zero();
//   double sumd = 0; std::vector<double> pp(3), ddd(3,0);
//   for(unsigned i=0; i<nactive; ++i) {
//     Vector distance = pbcDistance( getPosition(getNumberOfAtoms()-1), getPosition(active_list[i]) );
//     for(unsigned j=0; j<3; ++j) pp[j] = distance[j];
// 
//     // Now create the kernel and evaluate
//     KernelFunctions kernel( pp, bw, kerneltype, "DIAGONAL", 1.0 ); 
//     std::vector<Value*> fvals; kernel.normalize( fvals );
//     double newval = kernel.evaluate( pval, ddd, true );
//     if( getNumberOfArguments()==1 ) {
//     } else {
//       sumd += ddd[dir];
//       for(unsigned j=0; j<3; ++j) atom_deriv[i][j] = -ddd[j];
//       origind += -atom_deriv[i]; vir -= Tensor(atom_deriv[i],distance);
//     }
//   }
// 
//   // Add derivatives to atoms involved
//   Value* val=getPntrToComponent("qdist"); double prefactor =  root2 / sumd;
//   for(unsigned i=0; i<nactive; ++i) {
//     val->addDerivative( 3*active_list[i] + 0, -prefactor*atom_deriv[i][0] );
//     val->addDerivative( 3*active_list[i] + 1, -prefactor*atom_deriv[i][1] );
//     val->addDerivative( 3*active_list[i] + 2, -prefactor*atom_deriv[i][2] );
//   }
// 
//   // Add derivatives to atoms at origin
//   unsigned nbase = 3*(getNumberOfAtoms()-1);
//   val->addDerivative( nbase, -prefactor*origind[0] ); nbase++;
//   val->addDerivative( nbase, -prefactor*origind[1] ); nbase++;
//   val->addDerivative( nbase, -prefactor*origind[2] ); nbase++;
// 
//   // Add derivatives to virial
//   for(unsigned i=0; i<3; ++i) for(unsigned j=0; j<3; ++j) { val->addDerivative( nbase, -prefactor*vir(i,j) ); nbase++; }
}

}
}
