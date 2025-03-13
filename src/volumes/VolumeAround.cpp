/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2020 The plumed team
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
#include "tools/HistogramBead.h"
#include "ActionVolume.h"
#include "VolumeShortcut.h"

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
AROUND DATA=c ATOM=c1 XLOWER=-2.0 XUPPER=2.0 SIGMA=0.1 MEAN LABEL=s
\endplumedfile

*/
//+ENDPLUMEDOC

//+PLUMEDOC MCOLVAR AROUND_CALC
/*
Calculate a vector from the input positions with elements equal to one when the positions are in a particular part of the cell and elements equal to zero otherwise

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace volumes {

class VolumeAround {
public:
  bool dox, doy, doz;
  double sigma;
  double xlow, xhigh;
  double ylow, yhigh;
  double zlow, zhigh;
  std::string kerneltype;
  static void registerKeywords( Keywords& keys );
  void parseInput( ActionVolume<VolumeAround>* action );
  void setupRegions( ActionVolume<VolumeAround>* action, const Pbc& pbc, const std::vector<Vector>& positions ) {}
  static void parseAtoms( ActionVolume<VolumeAround>* action, std::vector<AtomNumber>& atom );
  VolumeAround& operator=( const VolumeAround& m ) {
    dox=m.dox;
    doy=m.doy;
    doz=m.doz;
    sigma=m.sigma;
    xlow=m.xlow;
    xhigh=m.xhigh;
    ylow=m.ylow;
    yhigh=m.yhigh;
    zlow=m.zlow;
    zhigh=m.zhigh;
    kerneltype=m.kerneltype;
    return *this;
  }
  static void calculateNumberInside( const VolumeInput& input, const VolumeAround& actioninput, VolumeOutput& output );
};

typedef ActionVolume<VolumeAround> Vola;
PLUMED_REGISTER_ACTION(Vola,"AROUND_CALC")
char glob_around[] = "AROUND";
typedef VolumeShortcut<glob_around> VolumeAroundShortcut;
PLUMED_REGISTER_ACTION(VolumeAroundShortcut,"AROUND")

void VolumeAround::registerKeywords( Keywords& keys ) {
  keys.setDisplayName("AROUND");
  keys.add("atoms","ORIGIN","the atom whose vicinity we are interested in examining");
  keys.add("atoms-2","ATOM","an alternative to ORIGIN");
  keys.add("compulsory","SIGMA","the width of the function to be used for kernel density estimation");
  keys.add("compulsory","KERNEL","gaussian","the type of kernel function to be used");
  keys.add("compulsory","XLOWER","0.0","the lower boundary in x relative to the x coordinate of the atom (0 indicates use full extent of box).");
  keys.add("compulsory","XUPPER","0.0","the upper boundary in x relative to the x coordinate of the atom (0 indicates use full extent of box).");
  keys.add("compulsory","YLOWER","0.0","the lower boundary in y relative to the y coordinate of the atom (0 indicates use full extent of box).");
  keys.add("compulsory","YUPPER","0.0","the upper boundary in y relative to the y coordinate of the atom (0 indicates use full extent of box).");
  keys.add("compulsory","ZLOWER","0.0","the lower boundary in z relative to the z coordinate of the atom (0 indicates use full extent of box).");
  keys.add("compulsory","ZUPPER","0.0","the upper boundary in z relative to the z coordinate of the atom (0 indicates use full extent of box).");
}

void VolumeAround::parseAtoms( ActionVolume<VolumeAround>* action, std::vector<AtomNumber>& atom ) {
  action->parseAtomList("ORIGIN",atom);
  if( atom.size()==0 ) {
    action->parseAtomList("ATOM",atom);
  }
  if( atom.size()!=1 ) {
    action->error("should only be one atom specified");
  }
  action->log.printf("  boundaries for region are calculated based on positions of atom : %d\n",atom[0].serial() );
}

void VolumeAround::parseInput( ActionVolume<VolumeAround>* action ) {
  action->parse("SIGMA",sigma);
  action->parse("KERNEL",kerneltype);
  dox=true;
  action->parse("XLOWER",xlow);
  action->parse("XUPPER",xhigh);
  doy=true;
  action->parse("YLOWER",ylow);
  action->parse("YUPPER",yhigh);
  doz=true;
  action->parse("ZLOWER",zlow);
  action->parse("ZUPPER",zhigh);
  if( xlow==0.0 && xhigh==0.0 ) {
    dox=false;
  }
  if( ylow==0.0 && yhigh==0.0 ) {
    doy=false;
  }
  if( zlow==0.0 && zhigh==0.0 ) {
    doz=false;
  }
  if( !dox && !doy && !doz ) {
    action->error("no subregion defined use XLOWER, XUPPER, YLOWER, YUPPER, ZLOWER, ZUPPER");
  }
  action->log.printf("  boundaries for region (region of interest about atom) : x %f %f, y %f %f, z %f %f \n",xlow,xhigh,ylow,yhigh,zlow,zhigh);
}

void VolumeAround::calculateNumberInside( const VolumeInput& input, const VolumeAround& actioninput, VolumeOutput& output ) {
  // Setup the histogram bead
  HistogramBead bead;
  bead.isNotPeriodic();
  bead.setKernelType( actioninput.kerneltype );

  // Calculate position of atom wrt to origin
  Vector fpos=input.pbc.distance( Vector(input.refpos[0][0],input.refpos[0][1],input.refpos[0][2]), Vector(input.cpos[0],input.cpos[1],input.cpos[2]) );
  double xcontr, ycontr, zcontr, xder, yder, zder;
  if( actioninput.dox ) {
    bead.set( actioninput.xlow, actioninput.xhigh, actioninput.sigma );
    xcontr=bead.calculate( fpos[0], xder );
  } else {
    xcontr=1.;
    xder=0.;
  }
  if( actioninput.doy ) {
    bead.set( actioninput.ylow, actioninput.yhigh, actioninput.sigma );
    ycontr=bead.calculate( fpos[1], yder );
  } else {
    ycontr=1.;
    yder=0.;
  }
  if( actioninput.doz ) {
    bead.set( actioninput.zlow, actioninput.zhigh, actioninput.sigma );
    zcontr=bead.calculate( fpos[2], zder );
  } else {
    zcontr=1.;
    zder=0.;
  }
  output.derivatives[0]=xder*ycontr*zcontr;
  output.derivatives[1]=xcontr*yder*zcontr;
  output.derivatives[2]=xcontr*ycontr*zder;
  // Add derivatives wrt to position of origin atom
  output.refders[0][0] = -output.derivatives[0];
  output.refders[0][1] = -output.derivatives[1];
  output.refders[0][2] = -output.derivatives[2];
  // Add virial contribution
  output.virial.set( 0, -Tensor(fpos,Vector(output.derivatives[0], output.derivatives[1], output.derivatives[2])) );
  output.values[0] = xcontr*ycontr*zcontr;
}

}
}
