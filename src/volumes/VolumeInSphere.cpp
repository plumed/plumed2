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
#include "VolumeShortcut.h"

//+PLUMEDOC VOLUMES INSPHERE
/*
This quantity can be used to calculate functions of the distribution of collective variables for the atoms that lie in a particular, user-specified part of of the cell.

This action can be used to calculate whether each of the atoms are within a particular part of the simulation box or not as illustrated by the following example:

```plumed
f: FIXEDATOM AT=0,0,0
a: INSPHERE ATOMS=1-100 CENTER=f RADIUS={GAUSSIAN R_0=1.5}
PRINT ARG=a FILE=colvar
```

The 100 elements of the vector `a` that is returned from the INSPHERE action in the above input are calculated using:

$$
w(x_i,y_i,z_i) = \sigma\left( \sqrt{ x_i^2 + y_i^2 + z_i^2}  \right)
$$

In this expression $x_i$, $y_i$ and $z_i$ are the components of the vector that connects the $i$th atom that was specified using the `ATOMS` keyword to the atom that was specified using the `CENTER` keyword and
$\sigma$ is one of the switching functions that is described that in the documentation for the action [LESS_THAN](LESS_THAN.md).  In other words,
$w(x_i,y_i,z_i)$ is 1 if atom $i$ is within a sphere with a radial extent that is determined by the parameters of the specified switching function
and zero otherwise.

If you want to caluclate and print the number of atoms in the sphere of interest you can use an input like the one shown below:

```plumed
f: FIXEDATOM AT=0,0,0
a: INSPHERE ATOMS=1-100 CENTER=f RADIUS={GAUSSIAN R_0=1.5}
s: SUM ARG=a PERIODIC=NO
PRINT ARG=s FILE=colvar
```

If by constrast you want to calculate and print the number of atoms that are not in the sphere of interest you OUTSIDE flag as shown below

```plumed
f: FIXEDATOM AT=0,0,0
a: INSPHERE ...
  ATOMS=1-100 CENTER=f
  RADIUS={GAUSSIAN R_0=1.5}
  OUTSIDE
...
s: SUM ARG=a PERIODIC=NO
PRINT ARG=s FILE=colvar
```

!!! note ""

    You can also calculate the average values of symmetry functions in the sphere of interest by using inputs similar to those described the documentation for the [AROUND](AROUND.md)
    action. In other words, you can swap out AROUND actions for an INSPHERE actions.  Also as with [AROUND](AROUND.md), earlier versions of PLUMED used a different syntax for doing these types of calculations, which can
    still be used with this new version of the command.  We strongly recommend using the newer syntax but if you are interested in the
    old syntax you can find more information in the old syntax section of the documentation for [AROUND](AROUND.md).

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace volumes {

struct VolumeInSphere {
#ifdef __PLUMED_USE_OPENACC
  SwitchingFunctionAccelerable switchingFunction;
#else
  SwitchingFunction switchingFunction;
#endif //__PLUMED_USE_OPENACC
  static void registerKeywords( Keywords& keys );
  void parseInput( ActionVolume<VolumeInSphere>* action );
  void setupRegions( ActionVolume<VolumeInSphere>* action, const Pbc& pbc, const std::vector<Vector>& positions ) {}
  static void parseAtoms( ActionVolume<VolumeInSphere>* action, std::vector<AtomNumber>& atom );
  static void calculateNumberInside( const VolumeInput& input, const VolumeInSphere& actioninput, VolumeOutput& output );
#ifdef __PLUMED_USE_OPENACC
  void toACCDevice() const {
#pragma acc enter data copyin(this[0:1])
    switchingFunction.toACCDevice();
  }
  void removeFromACCDevice() const {
    switchingFunction.removeFromACCDevice();
#pragma acc exit data delete(this[0:1])
  }
#endif //__PLUMED_USE_OPENACC
};

typedef ActionVolume<VolumeInSphere> Vols;
PLUMED_REGISTER_ACTION(Vols,"INSPHERE_CALC")
char glob_sphere[] = "INSPHERE";
typedef VolumeShortcut<glob_sphere> VolumeInSphereShortcut;
PLUMED_REGISTER_ACTION(VolumeInSphereShortcut,"INSPHERE")

void VolumeInSphere::registerKeywords( Keywords& keys ) {
  keys.setDisplayName("INSPHERE");
  keys.add("atoms","CENTER","the atom whose vicinity we are interested in examining");
  keys.addDeprecatedKeyword("ATOM","CENTER");
  keys.add("compulsory","RADIUS","the switching function that tells us the extent of the sphereical region of interest");
  keys.linkActionInDocs("RADIUS","LESS_THAN");
}

void VolumeInSphere::parseInput( ActionVolume<VolumeInSphere>* action ) {
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

  action->log.printf("  radius of sphere is given by %s \n",
                     switchingFunction.description().c_str() );
}

void VolumeInSphere::parseAtoms( ActionVolume<VolumeInSphere>* action, std::vector<AtomNumber>& atom ) {
  action->parseAtomList("CENTER",atom);
  if( atom.size()==0 ) {
    action->parseAtomList("ATOM",atom);
  }
  if( atom.size()!=1 ) {
    action->error("should only be one atom specified");
  }
  action->log.printf("  center of sphere is at position of atom : %d\n",atom[0].serial() );
}

void VolumeInSphere::calculateNumberInside( const VolumeInput& input,
    const VolumeInSphere& actioninput,
    VolumeOutput& output ) {
  // Calculate position of atom wrt to origin
  Vector fpos=input.pbc.distance( Vector(input.refpos[0][0],input.refpos[0][1],input.refpos[0][2]), Vector(input.cpos[0],input.cpos[1],input.cpos[2]) );
  double dfunc;
  output.values[0] = actioninput.switchingFunction.calculateSqr( fpos.modulo2(), dfunc );
  output.derivatives = dfunc*fpos;
  output.refders[0][0] = -output.derivatives[0];
  output.refders[0][1] = -output.derivatives[1];
  output.refders[0][2] = -output.derivatives[2];
  // Add a virial contribution
  output.virial.set( 0, -Tensor(fpos,Vector(output.derivatives[0], output.derivatives[1], output.derivatives[2])) );
}

}
}
