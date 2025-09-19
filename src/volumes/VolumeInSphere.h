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
#ifndef __PLUMED_volumes_VolumeInSphere_h
#define __PLUMED_volumes_VolumeInSphere_h
#include "tools/Pbc.h"
#include "tools/SwitchingFunction.h"
#include "ActionVolume.h"

namespace PLMD {
namespace volumes {

struct VolumeInSphere {
  using precision=double;
  typedef VolumeIn<precision> VolumeInput;
  typedef VolumeOut<precision> VolumeOutput;
#ifdef __PLUMED_HAS_OPENACC
  SwitchingFunctionAccelerable switchingFunction;
#else
  SwitchingFunction switchingFunction;
#endif //__PLUMED_HAS_OPENACC
  static void registerKeywords( Keywords& keys );
  template<typename TPM>
  void parseInput( ActionVolume<VolumeInSphere, TPM>* action );
  template<typename TPM>
  void setupRegions( ActionVolume<VolumeInSphere, TPM>* action,
                     const Pbc& pbc,
                     const std::vector<Vector>& positions ) {}
  template<typename TPM>
  static void parseAtoms( ActionVolume<VolumeInSphere, TPM>* action,
                          std::vector<AtomNumber>& atom );
  static void calculateNumberInside( const VolumeInput& input,
                                     const VolumeInSphere& actioninput,
                                     VolumeOutput& output );
#ifdef __PLUMED_HAS_OPENACC
  void toACCDevice() const {
#pragma acc enter data copyin(this[0:1])
    switchingFunction.toACCDevice();
  }
  void removeFromACCDevice() const {
    switchingFunction.removeFromACCDevice();
#pragma acc exit data delete(this[0:1])
  }
#endif //__PLUMED_HAS_OPENACC
};

void VolumeInSphere::registerKeywords( Keywords& keys ) {
  keys.setDisplayName("INSPHERE");
  keys.add("atoms","CENTER","the atom whose vicinity we are interested in examining");
  keys.addDeprecatedKeyword("ATOM","CENTER");
  keys.add("compulsory","RADIUS","the switching function that tells us the extent of the sphereical region of interest");
  keys.linkActionInDocs("RADIUS","LESS_THAN");
}

template<typename TPM>
void VolumeInSphere::parseInput( ActionVolume<VolumeInSphere, TPM>* action ) {
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

  action->log.printf("  radius of sphere is given by %s \n",
                     switchingFunction.description().c_str() );
}

template<typename TPM>
void VolumeInSphere::parseAtoms( ActionVolume<VolumeInSphere, TPM>* action,
                                 std::vector<AtomNumber>& atom ) {
  action->parseAtomList("CENTER",atom);
  if( atom.size()==0 ) {
    action->parseAtomList("ATOM",atom);
  }
  if( atom.size()!=1 ) {
    action->error("should only be one atom specified");
  }
  action->log.printf("  center of sphere is at position of atom : %d\n",atom[0].serial() );
}

void VolumeInSphere::calculateNumberInside( const VolumeInput& input,
    const VolumeInSphere& actioninput,
    VolumeOutput& output ) {
  // Calculate position of atom wrt to origin
  Vector fpos=input.pbc.distance( Vector(input.refpos[0][0],input.refpos[0][1],input.refpos[0][2]), Vector(input.cpos[0],input.cpos[1],input.cpos[2]) );
  double dfunc;
  output.values[0] = actioninput.switchingFunction.calculateSqr( fpos.modulo2(), dfunc );
  output.derivatives = dfunc*fpos;
  output.refders[0][0] = -output.derivatives[0];
  output.refders[0][1] = -output.derivatives[1];
  output.refders[0][2] = -output.derivatives[2];
  // Add a virial contribution
  output.virial.set( 0, -Tensor(fpos,Vector(output.derivatives[0], output.derivatives[1], output.derivatives[2])) );
}

}
}
#endif //__PLUMED_volumes_VolumeInSphere_h
