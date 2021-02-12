/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2020 The plumed team
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
#include "MultiColvarBase.h"
#include "AtomValuePack.h"
#include "core/ActionRegister.h"
#include "tools/KernelFunctions.h"
#include "tools/RootFindingBase.h"
#include "vesselbase/ValueVessel.h"

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
namespace multicolvar {

class DistanceFromContour : public MultiColvarBase {
private:
  unsigned dir;
  bool derivTime;
  double rcut2;
  double contour;
  double pbc_param;
  std::string kerneltype;
  std::vector<std::unique_ptr<Value>> pval;
  std::vector<double> bw, pos1, pos2, dirv, dirv2;
  std::vector<double> forces;
  std::vector<unsigned> perp_dirs;
  vesselbase::ValueVessel* myvalue_vessel;
  vesselbase::ValueVessel* myderiv_vessel;
  RootFindingBase<DistanceFromContour> mymin;
public:
  static void registerKeywords( Keywords& keys );
  explicit DistanceFromContour( const ActionOptions& );
  bool isDensity() const override { return true; }
  void calculate() override;
  unsigned getNumberOfQuantities() const override;
  bool isPeriodic() override { return false; }
  double compute( const unsigned& tindex, AtomValuePack& myatoms ) const override;
  double getDifferenceFromContour( const std::vector<double>& x, std::vector<double>& der );
// We need an apply action as we are using an independent value
  void apply() override;
};

PLUMED_REGISTER_ACTION(DistanceFromContour,"DISTANCE_FROM_CONTOUR")

void DistanceFromContour::registerKeywords( Keywords& keys ) {
  MultiColvarBase::registerKeywords( keys );
  keys.addOutputComponent("dist1","default","the distance between the reference atom and the nearest contour");
  keys.addOutputComponent("dist2","default","the distance between the reference atom and the other contour");
  keys.addOutputComponent("qdist","default","the differentiable (squared) distance between the two contours (see above)");
  keys.addOutputComponent("thickness","default","the distance between the two contours on the line from the reference atom");
  keys.add("compulsory","DATA","The input base multicolvar which is being used to calculate the contour");
  keys.add("atoms","ATOM","The atom whose perpendicular distance we are calculating from the contour");
  keys.add("compulsory","BANDWIDTH","the bandwidths for kernel density estimation");
  keys.add("compulsory","KERNEL","gaussian","the kernel function you are using.  More details on  the kernels available "
           "in plumed plumed can be found in \\ref kernelfunctions.");
  keys.add("compulsory","DIR","the direction perpendicular to the contour that you are looking for");
  keys.add("compulsory","CONTOUR","the value we would like for the contour");
  keys.add("compulsory","TOLERANCE","0.1","this parameter is used to manage periodic boundary conditions.  The problem "
           "here is that we can be between contours even when we are not within the membrane "
           "because of periodic boundary conditions.  When we are in the contour, however, we "
           "should have it so that the sums of the absolute values of the distances to the two "
           "contours is approximately the distance between the two contours.  There can be numerical errors in these calculations, however, so "
           "we specify a small tolerance here");
}

DistanceFromContour::DistanceFromContour( const ActionOptions& ao ):
  Action(ao),
  MultiColvarBase(ao),
  derivTime(false),
  bw(3),
  pos1(3,0.0),
  pos2(3,0.0),
  dirv(3,0.0),
  dirv2(3,0.0),
  perp_dirs(2),
  mymin(this)
{
  // Read in the multicolvar/atoms
  std::vector<AtomNumber> all_atoms; parse("TOLERANCE",pbc_param);
  bool read2 = parseMultiColvarAtomList("DATA", -1, all_atoms);
  if( !read2 ) error("missing DATA keyword");
  bool read1 = parseMultiColvarAtomList("ATOM", -1, all_atoms);
  if( !read1 ) error("missing ATOM keyword");
  if( all_atoms.size()!=1 ) error("should only be one atom specified");
  // Read in the center of the binding object
  log.printf("  computing distance of atom %d from contour \n",all_atoms[0].serial() );
  setupMultiColvarBase( all_atoms ); forces.resize( 3*getNumberOfAtoms() + 9 );
  if( getNumberOfBaseMultiColvars()!=1 ) error("should only be one input multicolvar");

  // Get the direction
  std::string ldir; parse("DIR",ldir );
  if( ldir=="x" ) { dir=0; perp_dirs[0]=1; perp_dirs[1]=2; dirv[0]=1; dirv2[0]=-1; }
  else if( ldir=="y" ) { dir=1; perp_dirs[0]=0; perp_dirs[1]=2; dirv[1]=1; dirv2[1]=-1; }
  else if( ldir=="z" ) { dir=2; perp_dirs[0]=0; perp_dirs[1]=1; dirv[2]=1; dirv2[2]=-1; }
  else error(ldir + " is not a valid direction use x, y or z");

  // Read in details of phase field construction
  parseVector("BANDWIDTH",bw); parse("KERNEL",kerneltype); parse("CONTOUR",contour);
  log.printf("  searching for contour in %s direction at %f in phase field for multicolvar %s \n",ldir.c_str(), contour, mybasemulticolvars[0]->getLabel().c_str() );
  log.printf("  constructing phase field using %s kernels with bandwidth (%f, %f, %f) \n",kerneltype.c_str(), bw[0], bw[1], bw[2] );

  // Now create a task list
  for(unsigned i=0; i<mybasemulticolvars[0]->getFullNumberOfTasks(); ++i) addTaskToList(i);
  // And a cutoff
  std::vector<double> pp( bw.size(),0 );
  KernelFunctions kernel( pp, bw, kerneltype, "DIAGONAL", 1.0 );
  double rcut = kernel.getCutoff( bw[0] );
  for(unsigned j=1; j<bw.size(); ++j) {
    if( kernel.getCutoff(bw[j])>rcut ) rcut=kernel.getCutoff(bw[j]);
  }
  rcut2=rcut*rcut;
  // Create the values
  addComponent("thickness"); componentIsNotPeriodic("thickness");
  addComponent("dist1"); componentIsNotPeriodic("dist1");
  addComponent("dist2"); componentIsNotPeriodic("dist2");
  addComponentWithDerivatives("qdist"); componentIsNotPeriodic("qdist");
  // Create sum vessels
  std::string fake_input; std::string deriv_input="COMPONENT=2";
  if( mybasemulticolvars[0]->isDensity() ) {
    addVessel( "SUM", fake_input, -1 ); addVessel( "SUM", deriv_input, -1 );
  } else {
    addVessel( "MEAN", fake_input, -1 ); addVessel( "MEAN", deriv_input, -1 );
  }
  // And convert to a value vessel so we can get the final value
  myvalue_vessel = dynamic_cast<vesselbase::ValueVessel*>( getPntrToVessel(0) );
  myderiv_vessel = dynamic_cast<vesselbase::ValueVessel*>( getPntrToVessel(1) );
  plumed_assert( myvalue_vessel && myderiv_vessel ); resizeFunctions();

  // Create the vector of values that holds the position
  for(unsigned i=0; i<3; ++i) pval.emplace_back( Tools::make_unique<Value>() );
}

unsigned DistanceFromContour::getNumberOfQuantities() const {
  return 3;
}

void DistanceFromContour::calculate() {
  // Check box is orthorhombic
  if( !getPbc().isOrthorombic() ) error("cell box must be orthorhombic");

  // The nanoparticle is at the origin of our coordinate system
  pos1[0]=pos1[1]=pos1[2]=0.0; pos2[0]=pos2[1]=pos2[2]=0.0;

  // Set bracket as center of mass of membrane in active region
  deactivateAllTasks();
  Vector myvec = getSeparation( getPosition(getNumberOfAtoms()-1), getPosition(0) ); pos2[dir]=myvec[dir];
  taskFlags[0]=1; double mindist = myvec.modulo2();
  for(unsigned j=1; j<getNumberOfAtoms()-1; ++j) {
    Vector distance=getSeparation( getPosition(getNumberOfAtoms()-1), getPosition(j) );
    double d2;
    if( (d2=distance[perp_dirs[0]]*distance[perp_dirs[0]])<rcut2 &&
        (d2+=distance[perp_dirs[1]]*distance[perp_dirs[1]])<rcut2 ) {
      d2+=distance[dir]*distance[dir];
      if( d2<mindist && fabs(distance[dir])>epsilon ) { pos2[dir]=distance[dir]; mindist = d2; }
      taskFlags[j]=1;
    }
  }
  lockContributors(); derivTime=false;
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
      pos1[dir] += sign*bw[dir]; fa = getDifferenceFromContour( pos1, faked );
      if( fa*fb<0 ) break;
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
      pos2[dir] += sign*bw[dir]; fc = getDifferenceFromContour( pos2, faked );
      if( fc*fb<0 ) break;
    }
    dirv2[dir] = ( pos1[dir] + dirv[dir] ) - pos2[dir];
  }

  // Now do a search for the two contours
  mymin.lsearch( dirv, pos1, &DistanceFromContour::getDifferenceFromContour );
  // Save the first value
  Vector root1; root1.zero(); root1[dir] = pval[dir]->get();
  mymin.lsearch( dirv2, pos2, &DistanceFromContour::getDifferenceFromContour );
  // Calculate the separation between the two roots using PBC
  Vector root2; root2.zero(); root2[dir]=pval[dir]->get();
  Vector sep = getSeparation( root1, root2 ); double spacing = fabs( sep[dir] ); plumed_assert( spacing>epsilon );
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
    Value* ival=myvalue_vessel->getFinalValue(); ival->clearDerivatives();
    std::vector<double> root1v(3); for(unsigned i=0; i<3; ++i) root1v[i]=root1[i];
    derivTime=true; std::vector<double> der(3); getDifferenceFromContour( root1v, der ); double prefactor;
    if( mybasemulticolvars[0]->isDensity() ) prefactor = root2[dir] / myderiv_vessel->getOutputValue();
    else plumed_error();
    Value* val=getPntrToComponent("qdist");
    for(unsigned i=0; i<val->getNumberOfDerivatives(); ++i) val->setDerivative( i, -prefactor*ival->getDerivative(i) );
    ival->clearDerivatives();
    std::vector<double> root2v(3); for(unsigned i=0; i<3; ++i) root2v[i]=root2[i];
    getDifferenceFromContour( root2v, der );
    if( mybasemulticolvars[0]->isDensity() ) prefactor = root1[dir] / myderiv_vessel->getOutputValue();
    else plumed_error();
    for(unsigned i=0; i<val->getNumberOfDerivatives(); ++i) val->addDerivative( i, -prefactor*ival->getDerivative(i) );
  }
}

double DistanceFromContour::getDifferenceFromContour( const std::vector<double>& x, std::vector<double>& der ) {
  std::string min, max;
  for(unsigned j=0; j<3; ++j) {
    Tools::convert( -0.5*getBox()(j,j), min );
    Tools::convert( +0.5*getBox()(j,j), max );
    pval[j]->setDomain( min, max ); pval[j]->set( x[j] );
  }
  runAllTasks();
  return myvalue_vessel->getOutputValue() - contour;
}

double DistanceFromContour::compute( const unsigned& tindex, AtomValuePack& myatoms ) const {
  Vector distance = getSeparation( getPosition(getNumberOfAtoms()-1), myatoms.getPosition(0) );
  std::vector<double> pp(3), der(3,0); for(unsigned j=0; j<3; ++j) pp[j] = distance[j];

  // convert the pointer once
  auto pval_ptr=Tools::unique2raw(pval);

  // Now create the kernel and evaluate
  KernelFunctions kernel( pp, bw, kerneltype, "DIAGONAL", 1.0 );
  kernel.normalize( pval_ptr ); double newval = kernel.evaluate( pval_ptr, der, true );

  if( mybasemulticolvars[0]->isDensity() ) {
    if( !doNotCalculateDerivatives() && derivTime ) {
      MultiValue& myvals=myatoms.getUnderlyingMultiValue();
      Vector vder; unsigned basen=myvals.getNumberOfDerivatives() - 12;
      for(unsigned i=0; i<3; ++i) {
        vder[i]=der[i]; myvals.addDerivative( 1, basen+i, vder[i] );
      }
      myatoms.setValue( 2, der[dir] );
      addAtomDerivatives( 1, 0, -vder, myatoms );
      myatoms.addBoxDerivatives( 1, Tensor(vder,distance) );
    }
    myatoms.setValue( 0, 1.0 ); return newval;
  }

  // This does the stuff for averaging
  myatoms.setValue( 0, newval );

  // This gets the average if we are using a phase field
  std::vector<double> cvals( mybasemulticolvars[0]->getNumberOfQuantities() );
  mybasedata[0]->retrieveValueWithIndex( tindex, false, cvals );
  return newval*cvals[0]*cvals[1];
}

void DistanceFromContour::apply() {
  if( getPntrToComponent("qdist")->applyForce( forces ) ) setForcesOnAtoms( forces );
}

}
}
