/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The plumed team
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

//+PLUMEDOC GENERIC SUBCELL
/*
This action is used in conjunction with \ref region keyword.  It is used to specify that 
a particular quantity is to be averaged for the atoms in a particular part of the cell.   
 
\par Examples

The following commands tell plumed to calculate the average coordination number for the atoms
that have x (in fractional coordinates) less than 0.5. The final value will be labeled r1_av.
\verbatim
COORDINATIONNUMBER SPECIES=1-100 R_0=1.0 LABEL=c
SUBCELL ARG=c XLOWER=0.0 XUPPER=0.5 SIGMA=0.1 LABEL=s
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
