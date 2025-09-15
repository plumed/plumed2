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
#ifndef __PLUMED_volumes_VolumeInCylinder_h
#define __PLUMED_volumes_VolumeInCylinder_h
#include "tools/Pbc.h"
#include "tools/SwitchingFunction.h"
#include "ActionVolume.h"
#include "tools/HistogramBead.h"


namespace PLMD {
namespace volumes {

class VolumeInCylinder {
public:
  using precision=double;
  typedef VolumeIn<precision> VolumeInput;
  typedef VolumeOut<precision> VolumeOutput;
  bool docylinder;
  double min, max, sigma;
  HistogramBead::KernelType kerneltype;
  std::array<unsigned,3> dir;
#ifdef __PLUMED_HAS_OPENACC
  SwitchingFunctionAccelerable switchingFunction;
#else
  SwitchingFunction switchingFunction;
#endif //__PLUMED_HAS_OPENACC
  static void registerKeywords( Keywords& keys );
  template<typename TPM>
  void parseInput( ActionVolume<VolumeInCylinder, TPM>* action );
  template<typename TPM>
  void setupRegions( ActionVolume<VolumeInCylinder, TPM>* action,
                     const Pbc& pbc,
                     const std::vector<Vector>& positions ) {}
  template<typename TPM>
  static void parseAtoms( ActionVolume<VolumeInCylinder, TPM>* action,
                          std::vector<AtomNumber>& atom );
  static void calculateNumberInside( const VolumeInput& input,
                                     const VolumeInCylinder& actioninput,
                                     VolumeOutput& output );
#ifdef __PLUMED_HAS_OPENACC
  void toACCDevice() const {
#pragma acc enter data copyin(this[0:1], \
docylinder, min, max, sigma, kerneltype, dir[0:3])
    switchingFunction.toACCDevice();
  }
  void removeFromACCDevice() const {
    switchingFunction.removeFromACCDevice();
#pragma acc exit data delete(dir[0:3], kerneltype, sigma, max, min, \
      docylinder, this[0:1])
  }
#endif //__PLUMED_HAS_OPENACC
};

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

template<typename TPM>
void VolumeInCylinder::parseInput( ActionVolume<VolumeInCylinder, TPM>* action ) {
  action->parse("SIGMA",sigma);
  std::string mykerneltype;
  action->parse("KERNEL",mykerneltype);
  kerneltype=HistogramBead::getKernelType(mykerneltype);
  std::string sdir;
  action->parse("DIRECTION",sdir);
  if( sdir=="X") {
    dir[0]=1;
    dir[1]=2;
    dir[2]=0;
  } else if( sdir=="Y") {
    dir[0]=0;
    dir[1]=2;
    dir[2]=1;
  } else if( sdir=="Z") {
    dir[0]=0;
    dir[1]=1;
    dir[2]=2;
  } else {
    action->error(sdir + "is not a valid direction.  Should be X, Y or Z");
  }
  action->log.printf("  cylinder's long axis is along %s axis\n",sdir.c_str() );

  std::string errors;
  std::string swinput;
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

template<typename TPM>
void VolumeInCylinder::parseAtoms( ActionVolume<VolumeInCylinder, TPM>* action,
                                   std::vector<AtomNumber>& atom ) {
  action->parseAtomList("CENTER",atom);
  if( atom.size()!=1 ) {
    action->error("should only be one atom specified");
  }
  action->log.printf("  center of cylinder is at position of atom : %d\n",atom[0].serial() );
}

void VolumeInCylinder::calculateNumberInside( const VolumeInput& input,
    const VolumeInCylinder& actioninput,
    VolumeOutput& output ) {
  // Calculate position of atom wrt to origin
  Vector fpos=input.pbc.distance( Vector(input.refpos[0][0],input.refpos[0][1],input.refpos[0][2]), Vector(input.cpos[0],input.cpos[1],input.cpos[2]) );

  double vcylinder, dcylinder;
  if( actioninput.docylinder ) {
    HistogramBead bead( actioninput.kerneltype,
                        actioninput.min, actioninput.max, actioninput.sigma );
    vcylinder=bead.calculate( fpos[actioninput.dir[2]], dcylinder );
  } else {
    vcylinder=1.0;
    dcylinder=0.0;
  }

  const double dd = fpos[actioninput.dir[0]]*fpos[actioninput.dir[0]] + fpos[actioninput.dir[1]]*fpos[actioninput.dir[1]];
  double dfunc;
  double vswitch = actioninput.switchingFunction.calculateSqr( dd, dfunc );
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
#endif //__PLUMED_volumes_VolumeInCylinder_h
