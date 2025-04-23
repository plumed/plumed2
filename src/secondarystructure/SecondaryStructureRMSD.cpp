/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2017-2023 The plumed team
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
// #ifdef __PLUMED_HAS_OPENACC
// #define __PLUMED_USE_OPENACC 1
// #endif //__PLUMED_HAS_OPENACC
#include "SecondaryStructureBase.h"
#include "core/ActionRegister.h"
#include "tools/RMSD.h"

//+PLUMEDOC MCOLVAR SECONDARY_STRUCTURE_RMSD
/*
Calclulate the distance between segments of a protein and a reference structure of interest

This action is used in the shortcuts [ALPHARMSD](ALPHARMSD.md), [ANTIBETARMSD](ANTIBETARMSD.md) and [PARABETARMSD](PARABETARMSD.md).  It calculates a
vector of RMSD values between a single reference multiple configurations and the instantaneous
positions of various groups of atoms.  For example, in the following input we define a single set of reference
set of coordinates for 3 atoms.

```plumed
c1: SECONDARY_STRUCTURE_RMSD BONDLENGTH=0.17 STRUCTURE1=1,0,0,0,1,0,0,0,1 SEGMENT1=1,2,3 SEGMENT2=4,5,6 SEGMENT3=7,8,9 SEGMENT4=10,11,12 TYPE=OPTIMAL
PRINT ARG=c1 FILE=colvar
```

A four dimensional vector is then returned that contains the RMSD distances between the 4 sets of atoms that were specified using the `SEGMENT` keywords
and the reference coordinates.  Notice that you can use multiple instances of the `STRUCTURE` keyword.  In general the the number of vectors output
is equal to the number of times the `STRUCTURE` keyword is used.

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace secondarystructure {

class SecondaryStructureRMSDInput {
/// The list of reference configurations
  std::vector<RMSD> myrmsd;
public:
/// The number of atoms
  std::size_t natoms;
/// The number of structures
  std::size_t nstructures;
/// The number of indices per task
  std::size_t nindices_per_task;
/// Are we operating without periodic boundary conditions
  bool nopbc;
/// Variables for strands cutoff
  bool align_strands;
/// The atoms involved in each of the secondary structure segments
  Matrix<unsigned> colvar_atoms;
//   void toACCDevice()const {
// #pragma acc enter data copyin(this[0:1], myrmsd[0:nstructures],natoms,nstructures,nopbc,align_strands)
//     colvar_atoms.toACCDevice();
//   }
//   void removeFromACCDevice() const {
//     colvar_atoms.removeFromACCDevice();
// #pragma acc exit data delete(align_strands,nopbc,nstructures,natoms,myrmsd[0:nstructures],this[0:1])

//   }
  static void calculateDistance( unsigned n, bool noderiv, const SecondaryStructureRMSDInput& actiondata, View<Vector> pos, ColvarOutput& output );
  void setReferenceStructure( const std::string& type, double bondlength, std::vector<Vector>& structure );
  SecondaryStructureRMSDInput& operator=( const SecondaryStructureRMSDInput& m ) {
    natoms = m.natoms;
    nstructures = m.nstructures;
    nindices_per_task = m.nindices_per_task;
    nopbc = m.nopbc;
    align_strands = m.align_strands;
    colvar_atoms=m.colvar_atoms;
    for(unsigned i=0; i<m.myrmsd.size(); ++i) {
      std::vector<Vector> mystruct( m.myrmsd[i].getReference() );
      setReferenceStructure( m.myrmsd[i].getMethod(), 0.0, mystruct );
    }
    return *this;
  }
};

typedef SecondaryStructureBase<SecondaryStructureRMSDInput> colv;
PLUMED_REGISTER_ACTION(colv,"SECONDARY_STRUCTURE_RMSD");

void SecondaryStructureRMSDInput::setReferenceStructure( const std::string& type, double bondlength, std::vector<Vector>& structure ) {
  Vector center;
  std::vector<double> align( structure.size(), 1.0 ), displace( structure.size(), 1.0 );
  for(unsigned i=0; i<structure.size(); ++i) {
    center+=structure[i]*align[i];
  }
  for(unsigned i=0; i<structure.size(); ++i) {
    structure[i] -= center;
  }
  RMSD newrmsd;
  newrmsd.clear();
  newrmsd.set(align,displace,structure,type,true,true);
  myrmsd.push_back( newrmsd );
  nstructures=myrmsd.size();
  natoms=structure.size();
}

void SecondaryStructureRMSDInput::calculateDistance( unsigned n,
    bool noderiv,
    const SecondaryStructureRMSDInput& actiondata,
    const View<Vector> pos,
    ColvarOutput& output ) {
  std::vector<Vector> myderivs( actiondata.natoms );
  output.values[n] = actiondata.myrmsd[n].calculate( pos, myderivs, false );

  if( noderiv ) {
    return;
  }
  Tensor v;
  v.zero();
  for(unsigned j=0; j<pos.size(); ++j) {
    output.derivs[n][j] = myderivs[j];
    v +=  (-1.0*Tensor( pos[j], myderivs[j] ));
  }
  output.virial.set( n, v );
}

}
}
