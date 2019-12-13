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

//+PLUMEDOC VOLUMES INSPHERE
/*
This quantity can be used to calculate functions of the distribution of collective variables for the atoms that lie in a particular, user-specified part of of the cell.

Each of the base quantities calculated by a multicolvar can can be assigned to a particular point in three
dimensional space. For example, if we have the coordination numbers for all the atoms in the
system each coordination number can be assumed to lie on the position of the central atom.
Because each base quantity can be assigned to a particular point in space we can calculate functions of the
distribution of base quantities in a particular part of the box by using:

\f[
\overline{s}_{\tau} = \frac{ \sum_i f(s_i) \sigma(r) }{ \sum_i \sigma(r) }
\f]

where the sum is over the collective variables, \f$s_i\f$, each of which can be thought to be at \f$ (x_i,y_i,z_i)\f$.
The function \f$\sigma\f$ is a \ref switchingfunction that acts on the distance between the point at which the
collective is located \f$(x_i,y_i,z_i)\f$ and the position of the atom that was specified using the ORIGIN keyword.
In other words:
\f[
r = sqrt{ ( x_i - x_0)^2 + ( y_i - y_0)^2 + ( z_i - z_0)^2}
\f]
In short this function, \f$\sigma(r_{xy})\f$, measures whether or not the CV is within a sphere that is
centered on the position of the atom specified using the keyword ORIGIN.

The function \f$(s_i)\f$ can be any of the usual LESS_THAN, MORE_THAN, WITHIN etc that are used in all other multicolvars.

When INCYLINDER is used with the \ref DENSITY action the number of atoms in the specified region is calculated

\par Examples

The input below can be use to calculate the average coordination numbers for those atoms that are within a sphere
of radius 1.5 nm that is centered on the position of atom 101.

\plumedfile
c1: COORDINATIONNUMBER SPECIES=1-100 SWITCH={RATIONAL R_0=0.1}
d2: INSPHERE ATOM=101 DATA=c1 RADIUS={TANH R_0=1.5} MEAN
PRINT ARG=d2.* FILE=colvar
\endplumedfile

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace multicolvar {

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

PLUMED_REGISTER_ACTION(VolumeInSphere,"INSPHERE")

void VolumeInSphere::registerKeywords( Keywords& keys ) {
  ActionVolume::registerKeywords( keys );
  keys.add("atoms","ATOM","the atom whose vicinity we are interested in examining");
  keys.add("compulsory","RADIUS","the switching function that tells us the extent of the spherical region of interest");
  keys.remove("SIGMA");
}

VolumeInSphere::VolumeInSphere(const ActionOptions& ao):
  Action(ao),
  ActionVolume(ao)
{
  std::vector<AtomNumber> atom;
  parseAtomList("ATOM",atom);
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
