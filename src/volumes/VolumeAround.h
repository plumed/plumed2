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
#ifndef __PLUMED_volumes_VolumeAround_h
#define __PLUMED_volumes_VolumeAround_h
#include "tools/Pbc.h"
#include "tools/HistogramBead.h"
#include "ActionVolume.h"


namespace PLMD {
namespace volumes {

class VolumeAround {
public:
  using precision=double;
  typedef VolumeIn<precision> VolumeInput;
  typedef VolumeOut<precision> VolumeOutput;
  bool dox{true}, doy{true}, doz{true};
  double sigma;
  double xlow{0.0}, xhigh{0.0};
  double ylow{0.0}, yhigh{0.0};
  double zlow{0.0}, zhigh{0.0};
  HistogramBead::KernelType kerneltype;
  static void registerKeywords( Keywords& keys );
  template<typename TPM>
  void parseInput( ActionVolume<VolumeAround, TPM>* action );
  template<typename TPM>
  void setupRegions( ActionVolume<VolumeAround, TPM>* action,
                     const Pbc& pbc,
                     const std::vector<Vector>& positions ) {}
  template<typename TPM>
  static void parseAtoms( ActionVolume<VolumeAround, TPM>* action,
                          std::vector<AtomNumber>& atom );
  static void calculateNumberInside( const VolumeInput& input,
                                     const VolumeAround& actioninput,
                                     VolumeOutput& output );
#ifdef __PLUMED_HAS_OPENACC
  void toACCDevice() const {
#pragma acc enter data copyin(this[0:1],dox,doy,doz,sigma,\
  xlow,xhigh,ylow,yhigh,zlow,zhigh,kerneltype)
  }
  void removeFromACCDevice() const {
#pragma acc exit data delete(kerneltype,zhigh,zlow,yhigh,ylow,xhigh,xlow,\
  sigma,doz,doy,dox,this[0:1])
  }
#endif //__PLUMED_HAS_OPENACC
};

void VolumeAround::registerKeywords( Keywords& keys ) {
  keys.setDisplayName("AROUND");
  keys.add("atoms","ORIGIN","the atom whose vicinity we are interested in examining");
  keys.addDeprecatedKeyword("ATOM","ORIGIN");
  keys.add("compulsory","SIGMA","the width of the function to be used for kernel density estimation");
  keys.add("compulsory","KERNEL","gaussian","the type of kernel function to be used");
  keys.add("compulsory","XLOWER","0.0","the lower boundary in x relative to the x coordinate of the atom (0 indicates use full extent of box).");
  keys.add("compulsory","XUPPER","0.0","the upper boundary in x relative to the x coordinate of the atom (0 indicates use full extent of box).");
  keys.add("compulsory","YLOWER","0.0","the lower boundary in y relative to the y coordinate of the atom (0 indicates use full extent of box).");
  keys.add("compulsory","YUPPER","0.0","the upper boundary in y relative to the y coordinate of the atom (0 indicates use full extent of box).");
  keys.add("compulsory","ZLOWER","0.0","the lower boundary in z relative to the z coordinate of the atom (0 indicates use full extent of box).");
  keys.add("compulsory","ZUPPER","0.0","the upper boundary in z relative to the z coordinate of the atom (0 indicates use full extent of box).");
}

template<typename TPM>
void VolumeAround::parseAtoms( ActionVolume<VolumeAround, TPM>* action, std::vector<AtomNumber>& atom ) {
  action->parseAtomList("ORIGIN",atom);
  if( atom.size()==0 ) {
    action->parseAtomList("ATOM",atom);
  }
  if( atom.size()!=1 ) {
    action->error("should only be one atom specified");
  }
  action->log.printf("  boundaries for region are calculated based on positions of atom : %d\n",atom[0].serial() );
}

template<typename TPM>
void VolumeAround::parseInput( ActionVolume<VolumeAround, TPM>* action ) {
  action->parse("SIGMA",sigma);
  std::string mykerneltype;
  action->parse("KERNEL",mykerneltype);
  kerneltype=HistogramBead::getKernelType(mykerneltype);
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

void VolumeAround::calculateNumberInside( const VolumeInput& input,
    const VolumeAround& actioninput,
    VolumeOutput& output ) {
  // Setup the histogram bead
  HistogramBead bead(actioninput.kerneltype, actioninput.xlow, actioninput.xhigh, actioninput.sigma );

  // Calculate position of atom wrt to origin
  Vector fpos=input.pbc.distance( Vector(input.refpos[0][0],input.refpos[0][1],input.refpos[0][2]),
                                  Vector(input.cpos[0],input.cpos[1],input.cpos[2]) );
  double xcontr=1.0;
  double xder=0.0;
  if( actioninput.dox ) {
    //bead parameters set in the constructor
    xcontr=bead.calculate( fpos[0], xder );
  }
  double ycontr=1.0;
  double yder=0.0;
  if( actioninput.doy ) {
    bead.set( actioninput.ylow, actioninput.yhigh, actioninput.sigma );
    ycontr=bead.calculate( fpos[1], yder );
  }
  double zcontr=1.0;
  double zder=0.0;
  if( actioninput.doz ) {
    bead.set( actioninput.zlow, actioninput.zhigh, actioninput.sigma );
    zcontr=bead.calculate( fpos[2], zder );
  }

  output.derivatives[0]=xder*ycontr*zcontr;
  output.derivatives[1]=xcontr*yder*zcontr;
  output.derivatives[2]=xcontr*ycontr*zder;
  // Add derivatives wrt to position of origin atom
  output.refders[0][0] = -output.derivatives[0];
  output.refders[0][1] = -output.derivatives[1];
  output.refders[0][2] = -output.derivatives[2];
  // Add virial contribution
  output.virial.set( 0, -Tensor(fpos,
                                Vector(output.derivatives[0], output.derivatives[1], output.derivatives[2])) );
  output.values[0] = xcontr*ycontr*zcontr;
}

}
}
#endif //__PLUMED_volumes_VolumeAround_h 
