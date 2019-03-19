/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2019 The plumed team
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
#include "ActionVolume.h"

//+PLUMEDOC VOLUMES AROUND
/*
This quantity can be used to calculate functions of the distribution of collective variables for the atoms that lie in a particular, user-specified part of of the cell.

Each of the base quantities calculated by a multicolvar can can be assigned to a particular point in three
dimensional space. For example, if we have the coordination numbers for all the atoms in the
system each coordination number can be assumed to lie on the position of the central atom.
Because each base quantity can be assigned to a particular point in space we can calculate functions of the
distribution of base quantities in a particular part of the box by using:

\f[
\overline{s}_{\tau} = \frac{ \sum_i f(s_i) w(x_i,y_i,z_i) }{ \sum_i w(x_i,y_i,z_i) }
\f]

where the sum is over the collective variables, \f$s_i\f$, each of which can be thought to be at \f$ (x_i,y_i,z_i)\f$.
The function \f$ w(x_i,y_i,z_i) \f$ measures whether or not the system is in the subregion of interest. It
is equal to:

\f[
w(x_i,y_i,z_i) = \int_{xl}^{xu} \int_{yl}^{yu} \int_{zl}^{zu} \textrm{d}x\textrm{d}y\textrm{d}z K\left( \frac{x - x_i}{\sigma} \right)K\left( \frac{y - y_i}{\sigma} \right)K\left( \frac{z - z_i}{\sigma} \right)
\f]

where \f$K\f$ is one of the kernel functions described on \ref histogrambead and \f$\sigma\f$ is a bandwidth parameter.
The function \f$(s_i)\f$ can be any of the usual LESS_THAN, MORE_THAN, WITHIN etc that are used in all other multicolvars.

When AROUND is used with the \ref DENSITY action the number of atoms in the specified region is calculated

\par Examples

The following commands tell plumed to calculate the average coordination number for the atoms
that have x (in fractional coordinates) within 2.0 nm of the com of mass c1. The final value will be labeled s.mean.
\plumedfile
COM ATOMS=1-100 LABEL=c1
COORDINATIONNUMBER SPECIES=1-100 R_0=1.0 LABEL=c
AROUND DATA=c ORIGIN=c1 XLOWER=-2.0 XUPPER=2.0 SIGMA=0.1 MEAN LABEL=s
\endplumedfile

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace multicolvar {

class VolumeAround : public ActionVolume {
private:
  Vector origin;
  bool dox, doy, doz;
  double xlow, xhigh;
  double ylow, yhigh;
  double zlow, zhigh;
public:
  static void registerKeywords( Keywords& keys );
  explicit VolumeAround(const ActionOptions& ao);
  void setupRegions();
  double calculateNumberInside( const Vector& cpos, Vector& derivatives, Tensor& vir, std::vector<Vector>& refders ) const ;
};

PLUMED_REGISTER_ACTION(VolumeAround,"AROUND")

void VolumeAround::registerKeywords( Keywords& keys ) {
  ActionVolume::registerKeywords( keys );
  keys.add("atoms","ATOM","the atom whose vicinity we are interested in examining");
  keys.add("compulsory","XLOWER","0.0","the lower boundary in x relative to the x coordinate of the atom (0 indicates use full extent of box).");
  keys.add("compulsory","XUPPER","0.0","the upper boundary in x relative to the x coordinate of the atom (0 indicates use full extent of box).");
  keys.add("compulsory","YLOWER","0.0","the lower boundary in y relative to the y coordinate of the atom (0 indicates use full extent of box).");
  keys.add("compulsory","YUPPER","0.0","the upper boundary in y relative to the y coordinate of the atom (0 indicates use full extent of box).");
  keys.add("compulsory","ZLOWER","0.0","the lower boundary in z relative to the z coordinate of the atom (0 indicates use full extent of box).");
  keys.add("compulsory","ZUPPER","0.0","the upper boundary in z relative to the z coordinate of the atom (0 indicates use full extent of box).");
}

VolumeAround::VolumeAround(const ActionOptions& ao):
  Action(ao),
  ActionVolume(ao)
{
  std::vector<AtomNumber> atom;
  parseAtomList("ATOM",atom);
  if( atom.size()!=1 ) error("should only be one atom specified");
  log.printf("  boundaries for region are calculated based on positions of atom : %d\n",atom[0].serial() );

  dox=true; parse("XLOWER",xlow); parse("XUPPER",xhigh);
  doy=true; parse("YLOWER",ylow); parse("YUPPER",yhigh);
  doz=true; parse("ZLOWER",zlow); parse("ZUPPER",zhigh);
  if( xlow==0.0 && xhigh==0.0 ) dox=false;
  if( ylow==0.0 && yhigh==0.0 ) doy=false;
  if( zlow==0.0 && zhigh==0.0 ) doz=false;
  if( !dox && !doy && !doz ) error("no subregion defined use XLOWER, XUPPER, YLOWER, YUPPER, ZLOWER, ZUPPER");
  log.printf("  boundaries for region (region of interest about atom) : x %f %f, y %f %f, z %f %f \n",xlow,xhigh,ylow,yhigh,zlow,zhigh);
  checkRead(); requestAtoms(atom);
}

void VolumeAround::setupRegions() { }

double VolumeAround::calculateNumberInside( const Vector& cpos, Vector& derivatives, Tensor& vir, std::vector<Vector>& refders ) const {
  // Setup the histogram bead
  HistogramBead bead; bead.isNotPeriodic(); bead.setKernelType( getKernelType() );

  // Calculate position of atom wrt to origin
  Vector fpos=pbcDistance( getPosition(0), cpos );
  double xcontr, ycontr, zcontr, xder, yder, zder;
  if( dox ) {
    bead.set( xlow, xhigh, getSigma() );
    xcontr=bead.calculate( fpos[0], xder );
  } else {
    xcontr=1.; xder=0.;
  }
  if( doy ) {
    bead.set( ylow, yhigh, getSigma() );
    ycontr=bead.calculate( fpos[1], yder );
  } else {
    ycontr=1.; yder=0.;
  }
  if( doz ) {
    bead.set( zlow, zhigh, getSigma() );
    zcontr=bead.calculate( fpos[2], zder );
  } else {
    zcontr=1.; zder=0.;
  }
  derivatives[0]=xder*ycontr*zcontr;
  derivatives[1]=xcontr*yder*zcontr;
  derivatives[2]=xcontr*ycontr*zder;
  // Add derivatives wrt to position of origin atom
  refders[0] = -derivatives;
  // Add virial contribution
  vir -= Tensor(fpos,derivatives);
  return xcontr*ycontr*zcontr;
}

}
}
