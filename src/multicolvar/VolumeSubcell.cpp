/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

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
#include "ActionVolume.h"

//+PLUMEDOC MCOLVARF SUBCELL
/*
This quantity can be used to calculate functions of the distribution of collective 
variables for the atoms that lie in a particular, user-specified part of of the cell.

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

When REGION is used with the \ref DENSITY action the number of atoms in the specified region is calculated  

\par Examples

The following commands tell plumed to calculate the average coordination number for the atoms
that have x (in fractional coordinates) less than 0.5. The final value will be labeled s.mean.
\verbatim
COORDINATIONNUMBER SPECIES=1-100 R_0=1.0 LABEL=c
SUBCELL ARG=c XLOWER=0.0 XUPPER=0.5 SIGMA=0.1 MEAN LABEL=s
\endverbatim

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace multicolvar {

class VolumeSubcell : public ActionVolume {
private:
  Tensor cellbox;
  bool dox, doy, doz;
  double xlow, xhigh, xsmearp;
  double ylow, yhigh, ysmearp;
  double zlow, zhigh, zsmearp;
public:
  static void registerKeywords( Keywords& keys );
  VolumeSubcell(const ActionOptions& ao);
  void setupRegion();
  bool derivativesOfFractionalCoordinates(){ return true; }
  double calculateNumberInside( const Vector& cpos, HistogramBead& bead, Vector& derivatives );
}; 

PLUMED_REGISTER_ACTION(VolumeSubcell,"SUBCELL")

void VolumeSubcell::registerKeywords( Keywords& keys ){
  ActionVolume::registerKeywords( keys );
  keys.add("compulsory","XLOWER","0.0","the lower boundary in x of the subcell in fractional coordinates.");
  keys.add("compulsory","XUPPER","1.0","the upper boundary in x of the subcell in fractional coordinates.");
  keys.add("compulsory","YLOWER","0.0","the lower boundary in y of the subcell in fractional coordinates.");
  keys.add("compulsory","YUPPER","1.0","the upper boundary in y of the subcell in fractional coordinates.");
  keys.add("compulsory","ZLOWER","0.0","the lower boundary in z of the subcell in fractional coordinates.");
  keys.add("compulsory","ZUPPER","1.0","the upper boundary in z of the subcell in fractional coordinates.");
}

VolumeSubcell::VolumeSubcell(const ActionOptions& ao):
Action(ao),
ActionVolume(ao),
xsmearp(0.0),
ysmearp(0.0),
zsmearp(0.0)
{

  dox=true; parse("XLOWER",xlow); parse("XUPPER",xhigh);
  doy=true; parse("YLOWER",ylow); parse("YUPPER",yhigh);
  doz=true; parse("ZLOWER",zlow); parse("ZUPPER",zhigh);
  if( xlow==0.0 && xhigh==1.0 ) dox=false;
  if( ylow==0.0 && yhigh==1.0 ) doy=false;
  if( zlow==0.0 && zhigh==1.0 ) doz=false; 
  if( !dox && !doy && !doz ) error("no subregion defined use XLOWER, XUPPER, YLOWER, YUPPER, ZLOWER, ZUPPER");
  log.printf("  boundaries for region (in fractional coordinates) : x %f %f, y %f %f, z %f %f \n",xlow,xhigh,ylow,yhigh,zlow,zhigh);
  checkRead();
}

void VolumeSubcell::setupRegion(){
  cellbox=getBox(); 
  if( dox ) xsmearp = getSigma() / ( sqrt( cellbox(0,0)*cellbox(0,0) + cellbox(1,0)*cellbox(1,0) + cellbox(2,0)*cellbox(2,0) ) );
  if( doy ) ysmearp = getSigma() / ( sqrt( cellbox(0,1)*cellbox(0,1) + cellbox(1,1)*cellbox(1,1) + cellbox(2,1)*cellbox(2,1) ) );
  if( doz ) zsmearp = getSigma() / ( sqrt( cellbox(0,2)*cellbox(0,2) + cellbox(1,2)*cellbox(1,2) + cellbox(2,2)*cellbox(2,2) ) );
}

double VolumeSubcell::calculateNumberInside( const Vector& cpos, HistogramBead& bead, Vector& derivatives ){
  double xcontr, ycontr, zcontr, xder, yder, zder; 
  if( dox ){
     bead.set( xlow, xhigh, xsmearp );
     xcontr=bead.calculate( cpos[0], xder ); 
  } else {
     xcontr=1.; xder=0.;
  }
  if( doy ){
     bead.set( ylow, yhigh, ysmearp );
     ycontr=bead.calculate( cpos[1], yder );
  } else {
     ycontr=1.; yder=0.;
  }
  if( doz ){
     bead.set( zlow, zhigh, zsmearp );
     zcontr=bead.calculate( cpos[2], zder );
  } else {
     zcontr=1.; zder=0.;
  }
  derivatives[0]=xder*ycontr*zcontr; 
  derivatives[1]=xcontr*yder*zcontr; 
  derivatives[2]=xcontr*ycontr*zder;
  return xcontr*ycontr*zcontr;
}

}
}
