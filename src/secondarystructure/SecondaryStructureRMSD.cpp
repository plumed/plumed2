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
Calclulate the RMSD between segments of a protein and a reference structure of interest

This action is used in the shortcuts [ALPHARMSD](ALPHARMSD.md), [ANTIBETARMSD](ANTIBETARMSD.md) and [PARABETARMSD](PARABETARMSD.md).  It calculates a
vector of [RMSD](RMSD.md) values between a single reference multiple configurations and the instantaneous
positions of various groups of atoms.  For example, in the following input we define a single set of reference
set of coordinates for 3 atoms.

```plumed
c1: SECONDARY_STRUCTURE_RMSD ...
    ATOMS=13-24
    STRUCTURE1=1,0,0,0,1,0,0,0,1
    SEGMENT1=0,1,2 SEGMENT2=3,4,5
    SEGMENT3=6,7,8 SEGMENT4=9,10,11
    TYPE=OPTIMAL
...
PRINT ARG=c1 FILE=colvar
```

A four dimensional vector is then returned that contains the RMSD distances between the 4 sets of 3 atoms that were specified using the `SEGMENT` keywords
and the reference coordinates.  Notice that you can use multiple instances of the `STRUCTURE` keyword.  In general the the number of vectors output
is equal to the number of times the `STRUCTURE` keyword is used.  Further note that the indices specified to the SEGMENT keywords give the indices in the list
of atoms specified by the ATOMS keyword.  `SEGMENT1=0,1,2` in the above input is measuring the distance between the instaneous positions of atoms 13, 14 and 15
and the reference structure.

## Periodic boundary conditions

If you turn off periodic boundary conditions by using the NOPBC flag as shown below:

```plumed
c1: SECONDARY_STRUCTURE_RMSD ...
    ATOMS=13-24 NOPBC
    STRUCTURE1=1,0,0,0,1,0,0,0,1
    SEGMENT1=0,1,2 SEGMENT2=3,4,5
    SEGMENT3=6,7,8 SEGMENT4=9,10,11
    TYPE=SIMPLE
...
PRINT ARG=c1 FILE=colvar
```

the distances between in the instaneous structure are evaluated in a way that does not take the periodic boundary conditions
into account.  Whenver this flag is __not__ present the distances in the instaneous structure are evaluated in a way that takes the periodic boundary conditions
into account.

If you use the `ALIGN_STRANDS` flag to evaluate the distances in a way that takes the periodic boundary conditions into accounts means that
calculating the distance between the positions of the 7th and 22nd atoms ($\mathbf{x}_7$ and $\mathbf{x}_{22}$) of the instanenous structures
using the periodic boundary conditions.  The 22 atom is the repositioned at:

$$
\mathbf{x}_{22} = \mathbf{x}_7 + ( \mathbf{x}_{22} - \mathbf{x}_7 )
$$

where the difference in the second term on the right of the equality sign is evaluated using the periodic boundary conditions. New positions are then determined
for atoms 16 through 30 by adding the difference (evaluated without PBC) between the new position for atom 22 after the above transformation has been perfomed and the
original position of this atom.  This treatment is designed for making sure that the two strands of the [ANTIBETARMSD](ANTIBETARMSD.md) and [PARABETARMSD](PARABETARMSD.md)
sheets are not broken by the periodic boundaries.

## Working with beta sheet like structures

As discussed in the documentation for [ANTIBETARMSD](ANTIBETARMSD.md) and [PARABETARMSD](PARABETARMSD.md) using the STRANDS_KEYWORD cutoff can speed up the calculation
dramatically.  The reason this keyword gives such dramatric speed ups is illustrated by the example input:

```plumed
ab_cut_dists: DISTANCE ...
  ATOMS1=19,69 ATOMS2=19,79
  ATOMS3=19,89 ATOMS4=19,99
  ATOMS5=19,109
...
ab_cut: CUSTOM ARG=ab_cut_dists FUNC=step(1-x) PERIODIC=NO
ab_rmsd: SECONDARY_STRUCTURE_RMSD ...
   TYPE=OPTIMAL ALIGN_STRANDS MASK=ab_cut
   ATOMS=7,9,11,15,16,17,19,21,25,26,27,29,31,35,36,37,39,41,45,46,47,49,51,55,56,57,59,61,65,66,67,69,71,75,76,77,79,81,85,86,87,89,91,95,96,97,99,101,105,106,107,109,111,115,116,117,119,121,125,126
   SEGMENT1=0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39
   SEGMENT2=0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44
   SEGMENT3=0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49
   SEGMENT4=0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54
   SEGMENT5=0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59
   STRUCTURE1=2.263,-3.795,1.722,2.493,-2.426,2.263,3.847,-1.838,1.761,1.301,-1.517,1.921,0.852,-1.504,0.739,0.818,-0.738,2.917,-0.299,0.243,2.748,-1.421,-0.076,3.757,0.273,1.68,2.854,0.902,1.993,3.888,0.119,2.532,1.813,0.683,3.916,1.68,1.58,3.94,0.395,-0.394,5.011,1.63,-1.459,4.814,0.982,-2.962,3.559,-1.359,-2.439,2.526,-2.287,-1.189,3.006,-3.087,-2.081,1.231,-1.52,-1.524,1.324,-0.409,-2.326,0.037,-2.095,-1.858,-1.269,-1.554,-3.053,-2.199,-1.291,-0.869,-1.949,-2.512,-1.255,-2.07,-3.71,0.326,-2.363,-2.072,1.405,-2.992,-2.872,2.699,-2.129,-2.917,1.745,-4.399,-2.33,1.899,-4.545,-1.102
...
ab_ltu: LESS_THAN ARG=ab_rmsd SWITCH={RATIONAL R_0=0.1 D_0=0.0 NN=8 MM=12} MASK=ab_cut
ab_lt: CUSTOM ARG=ab_ltu,ab_cut FUNC=x*y PERIODIC=NO
ab: SUM ARG=ab_lt PERIODIC=NO
PRINT ARG=ab FILE=colvar
```

The input above illustrates the input that is used by the shortcut that computes the [ANTIBETARMSD](ANTIBETARMSD.md). You can see that the distances between the atoms on the two strands of the
beta sheet are computed in a [DISTANCE](DISTANCE.md) action before we computed the SECONDARY_STRUCTURE_RMSD values.  These distances are then transformed by a Heavyside function that is 1 if
the distance is less than 1 nm and zero otherwise.  This vector of ones and zeros is then fed into the SECONDARY_STRUCTURE_RMSD action as a mask.  By doing this we ensure that SECONDARY_STRUCTURE_RMSD
values are only computed for parts of the strand that are close together. This avoids performing unecessarily expensive calculations and is also a reasonable thing to do as if the two strands are more than
1 nm appart the residues in question are not at all similar to a beta sheet structure.

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
PLUMED_REGISTER_ACTION(colv,"SECONDARY_STRUCTURE_RMSD")

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
