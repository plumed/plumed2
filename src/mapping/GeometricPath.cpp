/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016,2017 The plumed team
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
#include "core/ActionWithVector.h"
#include "core/ActionRegister.h"
#include "core/ActionSet.h"
#include "PathProjectionCalculator.h"

//+PLUMEDOC COLVAR GEOMETRIC_PATH
/*
Distance along and from a path calculated using geometric formulas

The Path Collective Variables developed by Branduardi and that are described in the first paper that is cited below alow one
to compute the progress along a high-dimensional path and the distance from the high-dimensional path.  The method introduced in that
paper is implemented in the shortcut [PATH](PATH.md). This action provides an implementation of the alternative method for calculating
the position on and distance from the path that was proposed by Diaz Leines and Ensing in the second paper that is cited below.  In their
method, the progress along the path $s$ is calculated using:

$$
s = i_2 + \textrm{sign}(i_2-i_1) \frac{ \sqrt{( \mathbf{v}_1\cdot\mathbf{v}_2 )^2 - |\mathbf{v}_3|^2(|\mathbf{v}_1|^2 - |\mathbf{v}_2|^2) } }{2|\mathbf{v}_3|^2} - \frac{\mathbf{v}_1\cdot\mathbf{v}_3 - |\mathbf{v}_3|^2}{2|\mathbf{v}_3|^2}
$$

where $\mathbf{v}_1$ and $\mathbf{v}_3$ are the vectors connecting the current position to the closest and second closest node of the path,
respectfully and $i_1$ and $i_2$ are the projections of the closest and second closest frames of the path. $\mathbf{v}_2$, meanwhile, is the
vector connecting the closest frame to the second closest frame.  The distance from the path, $z$ is calculated using:

$$
z = \sqrt{ \left[ |\mathbf{v}_1|^2 - |\mathbf{v}_2| \left( \frac{ \sqrt{( \mathbf{v}_1\cdot\mathbf{v}_2 )^2 - |\mathbf{v}_3|^2(|\mathbf{v}_1|^2 - |\mathbf{v}_2|^2) } }{2|\mathbf{v}_3|^2} - \frac{\mathbf{v}_1\cdot\mathbf{v}_3 - |\mathbf{v}_3|^2}{2|\mathbf{v}_3|^2} \right) \right]^2 }
$$

The symbols here are as they were for $s$.

Please note that the shortcut [GPATH](GPATH.md) provides an interface to this action that is less flexible than using this action directly but that is easier to use.

## Examples

The example input below shows how to use this action

```plumed
#SETTINGS INPUTFILES=regtest/mapping/rt-adapt/mypath.pdb
# Calculate the instaneous position in the high dimensional space
d1: DISTANCE ATOMS=1,2 COMPONENTS
# Read in the refence coordinates from the path
ref_d1x: PDB2CONSTANT REFERENCE=regtest/mapping/rt-adapt/mypath.pdb ARG=d1.x
ref_d1y: PDB2CONSTANT REFERENCE=regtest/mapping/rt-adapt/mypath.pdb ARG=d1.y
path: CONSTANT VALUES=-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5
# Now calculate the differences between our instantaneous position and the milestones on the high dimensional space
d1_dataP: DISPLACEMENT ARG2=d1.x,d1.y ARG1=ref_d1x,ref_d1y
# We need the vectors pointing from the reference positions to the instantaneous positions here which is
# minus the vector calculated above. Remember we can only specify the vectors in ARG1 for the DISPLACEMENT command
d1_data: CUSTOM ARG=d1_dataP FUNC=-x PERIODIC=NO
# And finally the path CV
pp: GEOMETRIC_PATH ARG=d1_data METRIC={DIFFERENCE} REFERENCE=ref_d1x,ref_d1y PROPERTY=path
PRINT ARG=d1.x,d1.y,pp.* FILE=colvar
```

The curved path here is defined using a series of points in a two dimensional space.  The input above includes actions ([PDB2CONSTANT](PDB2CONSTANT.md)) that set up
constant values that hold the values of the arguments for each of the milestones on the path.  A [DISPLACEMENT](DISPLACEMENT.md) shortcut is then used to calculate the vectors
connecting each of the path milestones to the instantenous coordinates of the system in the high dimensional space.  Notice, how the METRIC keyword is used in the input
to GEOMETRIC_PATH.  This keyword is necessary as the vectors connecting adjacent milestones on the path appear in the equations for $s$ and $z$ above.  The METRIC keyword tells
this action how distances between frames should be calculated.  There are essentially two options that we have tested here.  You can use [DIFFERENCE](DIFFERENCE.md) as has been done above.  In that
case the assumption is the path milestones are vectors of coordinates and the difference $\mathbf{v}$ between two of them, $\mathbf{m}_1$ and $\mathbf{m}_2$, can be calculated as:

$$
\mathbf{v} = \mathbf{m}_1 - \mathbf{m}_2
$$

with periodic variables treated appropriately.  Alternatively, the path milestones might be the positions of atoms and you may want to calculate differences between
them using an [RMSD](RMSD.md) action as has been done in the example below.

```plumed
##SETTINGS INPUTFILES=regtest/trajectories/path_msd/all.pdb
# Positions for all the frames in the path
path: CONSTANT VALUES=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42
# Calculate the vectors connecting the instantaneous position to the path milesones
vecs: RMSD DISPLACEMENT SQUARED REFERENCE=regtest/trajectories/path_msd/all.pdb TYPE=OPTIMAL
p1b: GEOMETRIC_PATH ARG=vecs.disp REFERENCE=vecs_ref METRIC={RMSD DISPLACEMENT TYPE=OPTIMAL ALIGN=1,1,1,1,1,1,1,1,1,1,1,1,1 DISPLACE=1,1,1,1,1,1,1,1,1,1,1,1,1} PROPERTY=path METRIC_COMPONENT=disp
PRINT ARG=p1b.* FILE=colvar_b STRIDE=1
```

Here an [RMSD](RMSD.md) action with label `vecs` is used to calculate the vector connecting the instantaneous position and each of the path milestones.  These vectors thus describe
how far each atom has been displaced in going from one structure to the other in a way that neglects translation of the center of mass and rotation of the reference frame.  When we
compute the formula aboves, which include the vectors that connect various milestone configurations on our path we must compute the vector connecting the milestones using the same
approach. The input associated with the METRIC keyword is thus the same input as was used in the RMSD action that was used to calculate the vectors connecting the milestones to the
instantaneous configuration.  Furthermore, we must also use METRIC_COMPONENT here to specify that it is the `disp` component of the [RMSD](RMSD.md) action that contains the vector of
interest.

It is worth noting that the METRIC keyword in the inputs above can accept inputs for one or multiple PLUMED actions.  There is thus flexibility to implement more complicated variants on the geometric path
implemented here by using the input for different actions as input to the METRIC keyword. If you are interested in using these features you can look at the code within `src/mapping/PathProjectionCalculator.h
and `src/mapping/PathProjectionCalculator.cpp` to see how this method called a separate PLUMED instance to calculate the vectors connecting the various milesones in the path.


*/
//+ENDPLUMEDOC

namespace PLMD {
namespace mapping {

class GeometricPath : public ActionWithVector {
private:
  PathProjectionCalculator path_projector;
public:
  static void registerKeywords(Keywords& keys);
  explicit GeometricPath(const ActionOptions&);
  unsigned getNumberOfDerivatives() override ;
  void calculate() override ;
};

PLUMED_REGISTER_ACTION(GeometricPath,"GEOMETRIC_PATH")

void GeometricPath::registerKeywords(Keywords& keys) {
  ActionWithVector::registerKeywords(keys);
  PathProjectionCalculator::registerKeywords(keys);
  keys.addInputKeyword("compulsory","PROPERTY","vector","the label of a value that contains the coordinates we are projecting these points onto");
  keys.addOutputComponent("s","default","scalar","the position on the path");
  keys.addOutputComponent("z","default","scalar","the distance from the path");
  keys.addDOI("10.1063/1.2432340");
  keys.addDOI("10.1103/PhysRevLett.109.020601");
  keys.needsAction("GEOMETRIC_PATH");
  keys.needsAction("PDB2CONSTANT");
}

GeometricPath::GeometricPath(const ActionOptions&ao):
  Action(ao),
  ActionWithVector(ao),
  path_projector(this) {
  // Get the coordinates in the low dimensional space
  std::vector<std::string> pcoord;
  parseVector("PROPERTY", pcoord );
  std::vector<Value*> theprop;
  ActionWithArguments::interpretArgumentList( pcoord, plumed.getActionSet(), this, theprop );
  if( theprop.size()!=1 ) {
    error("did not find property to project on");
  }
  if( theprop[0]->getNumberOfValues()!=getPntrToArgument(0)->getShape()[0] ) {
    error("mismatch between number of frames and property of interest");
  }
  log.printf("  projecting onto : %s \n", theprop[0]->getName().c_str() );
  std::vector<Value*> args( getArguments() );
  args.push_back( theprop[0] );
  requestArguments( args );
  // Create the values to store the output
  addComponentWithDerivatives("s");
  componentIsNotPeriodic("s");
  addComponentWithDerivatives("z");
  componentIsNotPeriodic("z");
}

unsigned GeometricPath::getNumberOfDerivatives() {
  return getPntrToArgument(0)->getShape()[0]*getPntrToArgument(0)->getShape()[1] + getPntrToArgument(1)->getShape()[0];
}

void GeometricPath::calculate() {
  unsigned k=0, iclose1=0, iclose2=0;
  double v1v1=0, v3v3=0;
  unsigned nrows = getPntrToArgument(0)->getShape()[0];
  unsigned ncols = getPntrToArgument(0)->getShape()[1];
  for(unsigned i=0; i<nrows; ++i) {
    double dist = 0;
    for(unsigned j=0; j<ncols; ++j) {
      double tmp = getPntrToArgument(0)->get(k);
      dist += tmp*tmp;
      k++;
    }
    if( i==0 ) {
      v1v1 = dist;
      iclose1 = 0;
    } else if( dist<v1v1 ) {
      v3v3=v1v1;
      v1v1=dist;
      iclose2=iclose1;
      iclose1=i;
    } else if( i==1 ) {
      v3v3=dist;
      iclose2=1;
    } else if( dist<v3v3 ) {
      v3v3=dist;
      iclose2=i;
    }
  }
  // And find third closest point
  int isign = iclose1 - iclose2;
  if( isign>1 ) {
    isign=1;
  } else if( isign<-1 ) {
    isign=-1;
  }
  int iclose3 = iclose1 + isign;
  unsigned ifrom=iclose1, ito=iclose3;
  if( iclose3<0 || static_cast<unsigned>(iclose3)>=nrows ) {
    ifrom=iclose2;
    ito=iclose1;
  }

  // And calculate projection of vector connecting current point to closest frame on vector connecting nearest two frames
  std::vector<double> displace;
  path_projector.getDisplaceVector( ifrom, ito, displace );
  double v2v2=0, v1v2=0;
  k=ncols*iclose1;
  for(unsigned i=0; i<displace.size(); ++i) {
    v2v2 += displace[i]*displace[i];
    v1v2 += displace[i]*getPntrToArgument(0)->get(k+i);
  }

  // This computes s value
  double spacing = getPntrToArgument(1)->get(iclose1) - getPntrToArgument(1)->get(iclose2);
  double root = sqrt( v1v2*v1v2 - v2v2 * ( v1v1 - v3v3) );
  double dx = 0.5 * ( (root + v1v2) / v2v2 - 1.);
  double path_s = getPntrToArgument(1)->get(iclose1) + spacing * dx;
  Value* sp = getPntrToComponent(0);
  sp->set( path_s );
  if( !doNotCalculateDerivatives() ) {
    for(unsigned i=0; i<ncols; ++i) {
      sp->addDerivative( ncols*iclose1 + i, 0.5*spacing*(v1v2*displace[i]/v2v2 - getPntrToArgument(0)->get(ncols*iclose1 + i))/root + 0.5*spacing*displace[i]/v2v2 );
      sp->addDerivative( ncols*iclose2 + i, 0.5*spacing*getPntrToArgument(0)->get(ncols*iclose2 + i)/root );
    }
  }

  // This computes z value
  path_projector.getDisplaceVector( iclose2, iclose1, displace );
  double v4v4=0, proj=0;
  k=ncols*iclose1;
  for(unsigned i=0; i<displace.size(); ++i) {
    v4v4 += displace[i]*displace[i];
    proj += displace[i]*getPntrToArgument(0)->get(k+i);
  }
  double path_z = v1v1 + dx*dx*v4v4 - 2*dx*proj;
  path_z = sqrt(path_z);
  Value* zp = getPntrToComponent(1);
  zp->set( path_z );
  if( !doNotCalculateDerivatives() ) {
    for(unsigned i=0; i<ncols; ++i) {
      zp->addDerivative( ncols*iclose1 + i, (1/path_z)*(getPntrToArgument(0)->get(ncols*iclose1 + i) +
                         (v4v4*dx-proj)*sp->getDerivative(ncols*iclose1 + i)/spacing -
                         dx*displace[i]) );
      zp->addDerivative( ncols*iclose2 + i, (v4v4*dx-proj)*sp->getDerivative(ncols*iclose2 + i)/(path_z*spacing) );
    }
  }
}

}
}
