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

//+PLUMEDOC VOLUMES INCYLINDER
/*
This quantity can be used to calculate functions of the distribution of collective variables for the atoms that lie in a particular, user-specified part of of the cell.

This action can be used to calculate whether each of the atoms are within a particular part of the simulation box or not as illustrated by the following example:

```plumed
f: FIXEDATOM AT=0,0,0
a: INCYLINDER ATOMS=1-100 CENTER=f DIRECTION=Z RADIUS={TANH R_0=1.5} SIGMA=0.1 LOWER=-1.0 UPPER=1.0 
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
a: INCYLINDER ATOMS=1-100 CENTER=f DIRECTION=X RADIUS={TANH R_0=1.5} SIGMA=0.1 LOWER=-1.0 UPPER=1.0 
s: SUM ARG=a PERIODIC=NO
PRINT ARG=s FILE=colvar
```

Notice that now the cylinder is centered on the `x` axis rather than the `z` axis as we have changed the input for the `DIRECTION` keyword.  

You can also calculate the average values of symmetry functions in the cylinder of interest by using inputs similar to those described the documentation for the [AROUND](AROUND.md) 
action. In other words, you can swap out AROUND actions for an INCLYLINDER actions. Also as with [AROUND](AROUND.md), earlier versions of PLUMED used a different syntax for doing these 
types of calculations, which can still be used with this new version of the command.  However, we strongly recommend using the newer syntax.


*/
//+ENDPLUMEDOC

//+PLUMEDOC MCOLVAR INCYLINDER_CALC
/*
Calculate a vector from the input positions with elements equal to one when the positions are in a particular part of the cell and elements equal to zero otherwise

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace volumes {

class VolumeInCylinder : public ActionVolume {
private:
  bool docylinder;
  Vector origin;
  HistogramBead bead;
  std::vector<unsigned> dir;
  SwitchingFunction switchingFunction;
public:
  static void registerKeywords( Keywords& keys );
  explicit VolumeInCylinder (const ActionOptions& ao);
  void setupRegions() override;
  double calculateNumberInside( const Vector& cpos, Vector& derivatives, Tensor& vir, std::vector<Vector>& refders ) const override;
};

PLUMED_REGISTER_ACTION(VolumeInCylinder,"INCYLINDER_CALC")
char glob_cylinder[] = "INCYLINDER";
typedef VolumeShortcut<glob_cylinder> VolumeInCylinderShortcut;
PLUMED_REGISTER_ACTION(VolumeInCylinderShortcut,"INCYLINDER")

void VolumeInCylinder::registerKeywords( Keywords& keys ) {
  ActionVolume::registerKeywords( keys ); keys.setDisplayName("INCYLINDER");
  keys.add("atoms","CENTER","the atom whose vicinity we are interested in examining");
  keys.add("compulsory","DIRECTION","the direction of the long axis of the cylinder. Must be x, y or z");
  keys.add("compulsory","RADIUS","a switching function that gives the extent of the cylinder in the plane perpendicular to the direction");
  keys.add("compulsory","LOWER","0.0","the lower boundary on the direction parallel to the long axis of the cylinder");
  keys.add("compulsory","UPPER","0.0","the upper boundary on the direction parallel to the long axis of the cylinder");
  keys.reset_style("SIGMA","optional"); keys.linkActionInDocs("RADIUS","LESS_THAN");
}

VolumeInCylinder::VolumeInCylinder(const ActionOptions& ao):
  Action(ao),
  ActionVolume(ao),
  docylinder(false)
{
  std::vector<AtomNumber> atom;
  parseAtomList("CENTER",atom);
  if( atom.size()!=1 ) error("should only be one atom specified");
  log.printf("  center of cylinder is at position of atom : %d\n",atom[0].serial() );

  std::string sdir; parse("DIRECTION",sdir);
  if( sdir=="X") {dir.push_back(1); dir.push_back(2); dir.push_back(0); }
  else if( sdir=="Y") {dir.push_back(0); dir.push_back(2); dir.push_back(1); }
  else if( sdir=="Z") {dir.push_back(0); dir.push_back(1); dir.push_back(2); }
  else { error(sdir + "is not a valid direction.  Should be X, Y or Z"); }
  log.printf("  cylinder's long axis is along %s axis\n",sdir.c_str() );

  std::string sw, errors; parse("RADIUS",sw);
  if(sw.length()==0) error("missing RADIUS keyword");
  switchingFunction.set(sw,errors);
  if( errors.length()!=0 ) error("problem reading RADIUS keyword : " + errors );
  log.printf("  radius of cylinder is given by %s \n", ( switchingFunction.description() ).c_str() );

  double min, max; parse("LOWER",min); parse("UPPER",max);
  if( min!=0.0 ||  max!=0.0 ) {
    if( min>max ) error("minimum of cylinder should be less than maximum");
    docylinder=true;
    log.printf("  cylinder extends from %f to %f along the %s axis\n",min,max,sdir.c_str() );
    bead.isNotPeriodic(); bead.setKernelType( getKernelType() ); bead.set( min, max, getSigma() );
  }

  checkRead(); requestAtoms(atom);
}

void VolumeInCylinder::setupRegions() { }

double VolumeInCylinder::calculateNumberInside( const Vector& cpos, Vector& derivatives, Tensor& vir, std::vector<Vector>& refders ) const {
  // Calculate position of atom wrt to origin
  Vector fpos=pbcDistance( getPosition(0), cpos );

  double vcylinder, dcylinder;
  if( docylinder ) {
    vcylinder=bead.calculate( fpos[dir[2]], dcylinder );
  } else {
    vcylinder=1.0; dcylinder=0.0;
  }

  const double dd = fpos[dir[0]]*fpos[dir[0]] + fpos[dir[1]]*fpos[dir[1]];
  double dfunc, vswitch = switchingFunction.calculateSqr( dd, dfunc );
  derivatives.zero(); double value=vswitch*vcylinder;
  derivatives[dir[0]]=vcylinder*dfunc*fpos[dir[0]];
  derivatives[dir[1]]=vcylinder*dfunc*fpos[dir[1]];
  derivatives[dir[2]]=vswitch*dcylinder;
  // Add derivatives wrt to position of origin atom
  refders[0] = -derivatives;
  // Add virial contribution
  vir -= Tensor(fpos,derivatives);
  return value;
}

}
}
