/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2020 The plumed team
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
#include "core/ActionRegister.h"
#include "tools/Pbc.h"
#include "tools/SwitchingFunction.h"
#include "ActionVolume.h"
#include "tools/HistogramBead.h"
#include "VolumeShortcut.h"

//+PLUMEDOC VOLUMES INCYLINDER
/*
This quantity can be used to calculate functions of the distribution of collective variables for the atoms that lie in a particular, user-specified part of of the cell.

This action can be used to calculate whether each of the atoms are within a particular part of the simulation box or not as illustrated by the following example:

```plumed
f: FIXEDATOM AT=0,0,0
a: INCYLINDER ...
   ATOMS=1-100 CENTER=f DIRECTION=Z
   RADIUS={TANH R_0=1.5} SIGMA=0.1
   LOWER=-1.0 UPPER=1.0
   KERNEL=gaussian
...
PRINT ARG=a FILE=colvar
```

The 100 elements of the vector `a` that is returned from the INCYLINDER action in the above input are calculated using:

$$
w(x_i,y_i,z_i) = s\left( r_{xy} \right) \int_{zl}^{zu} \textrm{d}z K\left( \frac{z - z_i}{\sigma} \right)
$$

where

$$
r_{xy} = \sqrt{ x_i^2 + y_i^2 }
$$

In these expressions $K$ is one of the kernel functions described in the documentation for the function [BETWEEN](BETWEEN.md), $\sigma$ is a bandwidth parameter and the limits
for the integrals are the values specified using the keywords `LOWER` and `UPPER`. $x_i$, $y_i$ and $z_i$ are then the components
of the vector that connects the $i$th atom that was specified using the `ATOMS` keyword to the atom that was specified using the `CENTER` keyword.  In other words,
$w(x_i,y_i,z_i)$ is 1 if the atom is within a cylinder that is centered on the $z$ axis that has an extent along the $z$ direction around the position of atom `f` that
is determined by the values specified by `LOWER` and `UPPER` keywords.  The radial extent of this cylinder is determined by the parameters of the switching function that is
specified using the `RADIUS` keyword.

If you want to caluclate and print the number of atoms in the cylinder of interest you can use an input like the one shown below:

```plumed
f: FIXEDATOM AT=0,0,0
a: INCYLINDER ...
  ATOMS=1-100 CENTER=f DIRECTION=X
  RADIUS={TANH R_0=1.5} SIGMA=0.1
  LOWER=-1.0 UPPER=1.0
...
s: SUM ARG=a PERIODIC=NO
PRINT ARG=s FILE=colvar
```

Alternatively, if you want to calculate an print the number of atoms that are not in the cylinder of interest you can use the OUTSIDE flag as shown below:

```plumed
f: FIXEDATOM AT=0,0,0
a: INCYLINDER ...
  ATOMS=1-100 CENTER=f DIRECTION=X
  RADIUS={TANH R_0=1.5} SIGMA=0.1
  LOWER=-1.0 UPPER=1.0
  OUTSIDE
...
s: SUM ARG=a PERIODIC=NO
PRINT ARG=s FILE=colvar
```

Notice that now the cylinder is centered on the `x` axis rather than the `z` axis as we have changed the input for the `DIRECTION` keyword.

!!! note ""

    You can also calculate the average values of symmetry functions in the cylinder of interest by using inputs similar to those described the documentation for the [AROUND](AROUND.md)
    action. In other words, you can swap out AROUND actions for an INCLYLINDER actions. Also as with [AROUND](AROUND.md), earlier versions of PLUMED used a different syntax for doing these
    types of calculations, which can still be used with this new version of the command.  We strongly recommend using the newer syntax but if you are interested in the
    old syntax you can find more information in the old syntax section of the documentation for [AROUND](AROUND.md).


*/
//+ENDPLUMEDOC

namespace PLMD {
namespace volumes {

class VolumeInCylinder {
public:
  bool docylinder;
  double min, max, sigma;
  HistogramBead::KernelType kerneltype;
  std::array<unsigned,3> dir;
#ifdef __PLUMED_USE_OPENACC
  SwitchingFunctionAccelerable switchingFunction;
#else
  SwitchingFunction switchingFunction;
#endif //__PLUMED_USE_OPENACC
  static void registerKeywords( Keywords& keys );
  void parseInput( ActionVolume<VolumeInCylinder>* action );
  void setupRegions( ActionVolume<VolumeInCylinder>* action, const Pbc& pbc, const std::vector<Vector>& positions ) {}
  static void parseAtoms( ActionVolume<VolumeInCylinder>* action, std::vector<AtomNumber>& atom );
  static void calculateNumberInside( const VolumeInput& input, const VolumeInCylinder& actioninput, VolumeOutput& output );
#ifdef __PLUMED_USE_OPENACC
  void toACCDevice() const {
#pragma acc enter data copyin(this[0:1], \
docylinder, min, max, sigma, kerneltype, dir[0:3])
    switchingFunction.toACCDevice();
  }
  void removeFromACCDevice() const {
    switchingFunction.removeFromACCDevice();
#pragma acc exit data delete(dir[0:3], kerneltype, sigma, max, min, \
      docylinder, this[0:1])
  }
#endif //__PLUMED_USE_OPENACC
};

typedef ActionVolume<VolumeInCylinder> Volc;
PLUMED_REGISTER_ACTION(Volc,"INCYLINDER_CALC")
char glob_cylinder[] = "INCYLINDER";
typedef VolumeShortcut<glob_cylinder> VolumeInCylinderShortcut;
PLUMED_REGISTER_ACTION(VolumeInCylinderShortcut,"INCYLINDER")

void VolumeInCylinder::registerKeywords( Keywords& keys ) {
  keys.setDisplayName("INCYLINDER");
  keys.add("atoms","CENTER","the atom whose vicinity we are interested in examining");
  keys.add("optional","SIGMA","the width of the function to be used for kernel density estimation");
  keys.add("compulsory","KERNEL","gaussian","the type of kernel function to be used");
  keys.add("compulsory","DIRECTION","the direction of the long axis of the cylinder. Must be x, y or z");
  keys.add("compulsory","RADIUS","a switching function that gives the extent of the cylinder in the plane perpendicular to the direction");
  keys.add("compulsory","LOWER","0.0","the lower boundary on the direction parallel to the long axis of the cylinder");
  keys.add("compulsory","UPPER","0.0","the upper boundary on the direction parallel to the long axis of the cylinder");
  keys.linkActionInDocs("RADIUS","LESS_THAN");
}

void VolumeInCylinder::parseInput( ActionVolume<VolumeInCylinder>* action ) {
  action->parse("SIGMA",sigma);
  std::string mykerneltype;
  action->parse("KERNEL",mykerneltype);
  kerneltype=HistogramBead::getKernelType(mykerneltype);
  std::string sdir;
  action->parse("DIRECTION",sdir);
  if( sdir=="X") {
    dir[0]=1;
    dir[1]=2;
    dir[2]=0;
  } else if( sdir=="Y") {
    dir[0]=0;
    dir[1]=2;
    dir[2]=1;
  } else if( sdir=="Z") {
    dir[0]=0;
    dir[1]=1;
    dir[2]=2;
  } else {
    action->error(sdir + "is not a valid direction.  Should be X, Y or Z");
  }
  action->log.printf("  cylinder's long axis is along %s axis\n",sdir.c_str() );

  std::string errors;
  std::string swinput;
  action->parse("RADIUS",swinput);
  if(swinput.length()==0) {
    action->error("missing RADIUS keyword");
  }
  switchingFunction.set(swinput,errors);
  if( errors.length()!=0 ) {
    action->error("problem reading RADIUS keyword : " + errors );
  }
  action->log.printf("  radius of cylinder is given by %s \n", ( switchingFunction.description() ).c_str() );

  docylinder=false;
  action->parse("LOWER",min);
  action->parse("UPPER",max);
  if( min!=0.0 ||  max!=0.0 ) {
    if( min>max ) {
      action->error("minimum of cylinder should be less than maximum");
    }
    docylinder=true;
    action->log.printf("  cylinder extends from %f to %f along the %s axis\n",min,max,sdir.c_str() );
  }
}

void VolumeInCylinder::parseAtoms( ActionVolume<VolumeInCylinder>* action, std::vector<AtomNumber>& atom ) {
  action->parseAtomList("CENTER",atom);
  if( atom.size()!=1 ) {
    action->error("should only be one atom specified");
  }
  action->log.printf("  center of cylinder is at position of atom : %d\n",atom[0].serial() );
}

void VolumeInCylinder::calculateNumberInside( const VolumeInput& input, const VolumeInCylinder& actioninput, VolumeOutput& output ) {
  // Calculate position of atom wrt to origin
  Vector fpos=input.pbc.distance( Vector(input.refpos[0][0],input.refpos[0][1],input.refpos[0][2]), Vector(input.cpos[0],input.cpos[1],input.cpos[2]) );

  double vcylinder, dcylinder;
  if( actioninput.docylinder ) {
    HistogramBead bead( actioninput.kerneltype,
                        actioninput.min, actioninput.max, actioninput.sigma );
    vcylinder=bead.calculate( fpos[actioninput.dir[2]], dcylinder );
  } else {
    vcylinder=1.0;
    dcylinder=0.0;
  }

  const double dd = fpos[actioninput.dir[0]]*fpos[actioninput.dir[0]] + fpos[actioninput.dir[1]]*fpos[actioninput.dir[1]];
  double dfunc;
  double vswitch = actioninput.switchingFunction.calculateSqr( dd, dfunc );
  output.values[0]=vswitch*vcylinder;
  output.derivatives[actioninput.dir[0]]=vcylinder*dfunc*fpos[actioninput.dir[0]];
  output.derivatives[actioninput.dir[1]]=vcylinder*dfunc*fpos[actioninput.dir[1]];
  output.derivatives[actioninput.dir[2]]=vswitch*dcylinder;
  // Add derivatives wrt to position of origin atom
  output.refders[0][0] = -output.derivatives[0];
  output.refders[0][1] = -output.derivatives[1];
  output.refders[0][2] = -output.derivatives[2];
  // Add virial contribution
  output.virial.set( 0, -Tensor(fpos,Vector(output.derivatives[0], output.derivatives[1], output.derivatives[2])) );
}

}
}
