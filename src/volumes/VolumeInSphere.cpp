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

You can also calculate the average values of symmetry functions in the sphere of interest by using inputs similar to those described the documentation for the [AROUND](AROUND.md) 
action. In other words, you can swap out AROUND actions for an INSPHERE actions.  Also as with [AROUND](AROUND.md), earlier versions of PLUMED used a different syntax for doing these types of calculations, which can 
still be used with this new version of the command.  However, we strongly recommend using the newer syntax.


*/
//+ENDPLUMEDOC

//+PLUMEDOC MCOLVAR INSPHERE_CALC
/*
Calculate a vector from the input positions with elements equal to one when the positions are in a particular part of the cell and elements equal to zero otherwise

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace volumes {

class VolumeInSphere : public ActionVolume {
private:
  Vector origin;
  SwitchingFunction switchingFunction;
public:
  static void registerKeywords( Keywords& keys );
  explicit VolumeInSphere(const ActionOptions& ao);
  void setupRegions() override;
  double calculateNumberInside( const Vector& cpos, Vector& derivatives, Tensor& vir, std::vector<Vector>& refders ) const override;
};

PLUMED_REGISTER_ACTION(VolumeInSphere,"INSPHERE_CALC")
char glob_sphere[] = "INSPHERE";
typedef VolumeShortcut<glob_sphere> VolumeInSphereShortcut;
PLUMED_REGISTER_ACTION(VolumeInSphereShortcut,"INSPHERE")

void VolumeInSphere::registerKeywords( Keywords& keys ) {
  ActionVolume::registerKeywords( keys ); keys.setDisplayName("INSPHERE");
  keys.add("atoms","CENTER","the atom whose vicinity we are interested in examining");
  keys.add("atoms-2","ATOM","the atom whose vicinity we are interested in examining");
  keys.add("compulsory","RADIUS","the switching function that tells us the extent of the sphereical region of interest");
  keys.remove("SIGMA"); keys.linkActionInDocs("RADIUS","LESS_THAN");
}

VolumeInSphere::VolumeInSphere(const ActionOptions& ao):
  Action(ao),
  ActionVolume(ao)
{
  std::vector<AtomNumber> atom; parseAtomList("CENTER",atom);
  if( atom.size()==0 ) parseAtomList("ATOM",atom);
  if( atom.size()!=1 ) error("should only be one atom specified");
  log.printf("  center of sphere is at position of atom : %d\n",atom[0].serial() );

  std::string sw, errors; parse("RADIUS",sw);
  if(sw.length()==0) error("missing RADIUS keyword");
  switchingFunction.set(sw,errors);
  if( errors.length()!=0 ) error("problem reading RADIUS keyword : " + errors );
  log.printf("  radius of sphere is given by %s \n", ( switchingFunction.description() ).c_str() );

  checkRead(); requestAtoms(atom);
}

void VolumeInSphere::setupRegions() { }

double VolumeInSphere::calculateNumberInside( const Vector& cpos, Vector& derivatives, Tensor& vir, std::vector<Vector>& refders ) const {
  // Calculate position of atom wrt to origin
  Vector fpos=pbcDistance( getPosition(0), cpos );
  double dfunc, value = switchingFunction.calculateSqr( fpos.modulo2(), dfunc );
  derivatives.zero(); derivatives = dfunc*fpos; refders[0] = -derivatives;
  // Add a virial contribution
  vir -= Tensor(fpos,derivatives);
  return value;
}

}
}
