/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2019 The plumed team
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

//+PLUMEDOC VOLUMES INCYLINDER
/*
This quantity can be used to calculate functions of the distribution of collective variables for the atoms that lie in a particular, user-specified part of of the cell.

Each of the base quantities calculated by a multicolvar can can be assigned to a particular point in three
dimensional space. For example, if we have the coordination numbers for all the atoms in the
system each coordination number can be assumed to lie on the position of the central atom.
Because each base quantity can be assigned to a particular point in space we can calculate functions of the
distribution of base quantities in a particular part of the box by using:

\f[
\overline{s}_{\tau} = \frac{ \sum_i f(s_i) \sigma(r_{xy}) }{ \sum_i \sigma(r_{xy}) }
\f]

where the sum is over the collective variables, \f$s_i\f$, each of which can be thought to be at \f$ (x_i,y_i,z_i)\f$.
The function \f$\sigma\f$ is a \ref switchingfunction that acts on the distance between the point at which the
collective is located \f$(x_i,y_i,z_i)\f$ and the position of the atom that was specified using the ORIGIN keyword
projected in the xy plane if DIRECTION=z is used.  In other words:
\f[
r_{xy} = sqrt{ ( x_i - x_0)^2 + ( y_i - y_0)^2 }
\f]
In short this function, \f$\sigma(r_{xy})\f$, measures whether or not the CV is within a cylinder that
runs along the axis specified using the DIRECTION keyword and that is centered on the position of the atom specified using
ORIGIN.

The function \f$(s_i)\f$ can be any of the usual LESS_THAN, MORE_THAN, WITHIN etc that are used in all other multicolvars.

When INCYLINDER is used with the \ref DENSITY action the number of atoms in the specified region is calculated

\par Examples

The input below can be use to calculate the average coordination numbers for those atoms that are within a cylindrical tube
of radius 1.5 nm that is centered on the position of atom 101 and that has its long axis parallel to the z-axis.

\plumedfile
c1: COORDINATIONNUMBER SPECIES=1-100 SWITCH={RATIONAL R_0=0.1}
d2: INCYLINDER ATOM=101 DATA=d1 DIRECTION=Z RADIUS={TANH R_0=1.5} SIGMA=0.1 LOWER=-0.1 UPPER=0.1 MEAN
PRINT ARG=d2.* FILE=colvar
\endplumedfile

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace multicolvar {

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
  void setupRegions();
  double calculateNumberInside( const Vector& cpos, Vector& derivatives, Tensor& vir, std::vector<Vector>& refders ) const ;
};

PLUMED_REGISTER_ACTION(VolumeInCylinder,"INCYLINDER")

void VolumeInCylinder::registerKeywords( Keywords& keys ) {
  ActionVolume::registerKeywords( keys );
  keys.add("atoms","ATOM","the atom whose vicinity we are interested in examining");
  keys.add("compulsory","DIRECTION","the direction of the long axis of the cylinder. Must be x, y or z");
  keys.add("compulsory","RADIUS","a switching function that gives the extent of the cylinder in the plane perpendicular to the direction");
  keys.add("compulsory","LOWER","0.0","the lower boundary on the direction parallel to the long axis of the cylinder");
  keys.add("compulsory","UPPER","0.0","the upper boundary on the direction parallel to the long axis of the cylinder");
  keys.reset_style("SIGMA","optional");
}

VolumeInCylinder::VolumeInCylinder(const ActionOptions& ao):
  Action(ao),
  ActionVolume(ao),
  docylinder(false)
{
  std::vector<AtomNumber> atom;
  parseAtomList("ATOM",atom);
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
