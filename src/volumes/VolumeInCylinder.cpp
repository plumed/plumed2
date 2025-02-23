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
#include "tools/HistogramBead.h"
#include "VolumeShortcut.h"

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
d2: INCYLINDER ATOM=101 DATA=c1 DIRECTION=Z RADIUS={TANH R_0=1.5} SIGMA=0.1 LOWER=-0.1 UPPER=0.1 MEAN
PRINT ARG=d2.* FILE=colvar
\endplumedfile

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

class VolumeInCylinder {
public:
  bool docylinder;
  double min, max, sigma;
  std::string kerneltype;
  std::string swinput;
  std::vector<unsigned> dir;
  SwitchingFunction switchingFunction;
  static void registerKeywords( Keywords& keys );
  void parseInput( ActionVolume<VolumeInCylinder>* action );
  void setupRegions( ActionVolume<VolumeInCylinder>* action, const Pbc& pbc, const std::vector<Vector>& positions ) {}
  static void parseAtoms( ActionVolume<VolumeInCylinder>* action, std::vector<AtomNumber>& atom );
  VolumeInCylinder& operator=( const VolumeInCylinder& m ) {
    docylinder=m.docylinder;
    min=m.min;
    max=m.max;
    sigma=m.sigma;
    kerneltype=m.kerneltype;
    dir=m.dir;
    swinput=m.swinput;
    std::string errors;
    switchingFunction.set(swinput,errors);
    return *this;
  }
  static void calculateNumberInside( const VolumeInput& input, const VolumeInCylinder& actioninput, VolumeOutput& output );
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
  action->parse("KERNEL",kerneltype);
  std::string sdir;
  action->parse("DIRECTION",sdir);
  if( sdir=="X") {
    dir.push_back(1);
    dir.push_back(2);
    dir.push_back(0);
  } else if( sdir=="Y") {
    dir.push_back(0);
    dir.push_back(2);
    dir.push_back(1);
  } else if( sdir=="Z") {
    dir.push_back(0);
    dir.push_back(1);
    dir.push_back(2);
  } else {
    action->error(sdir + "is not a valid direction.  Should be X, Y or Z");
  }
  action->log.printf("  cylinder's long axis is along %s axis\n",sdir.c_str() );

  std::string errors;
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
    HistogramBead bead;
    bead.isNotPeriodic();
    bead.setKernelType( actioninput.kerneltype );
    bead.set( actioninput.min, actioninput.max, actioninput.sigma );
    vcylinder=bead.calculate( fpos[actioninput.dir[2]], dcylinder );
  } else {
    vcylinder=1.0;
    dcylinder=0.0;
  }

  const double dd = fpos[actioninput.dir[0]]*fpos[actioninput.dir[0]] + fpos[actioninput.dir[1]]*fpos[actioninput.dir[1]];
  double dfunc, vswitch = actioninput.switchingFunction.calculateSqr( dd, dfunc );
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
