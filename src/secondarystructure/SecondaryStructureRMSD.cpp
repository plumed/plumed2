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
#include "SecondaryStructureRMSD.h"
#include "core/ActionRegister.h"

//+PLUMEDOC MCOLVAR SECONDARY_STRUCTURE_RMSD
/*
Calculate the RMSD between segments of a protein and a reference structure of interest

This action is used in the shortcuts [ALPHARMSD](ALPHARMSD.md), [ANTIBETARMSD](ANTIBETARMSD.md) and [PARABETARMSD](PARABETARMSD.md).  It calculates a
vector of [RMSD](RMSD.md) values between one or more reference configurations and the instantaneous
positions of various groups of atoms.  For example, in the following input we define a single reference
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
and the reference coordinates.  Notice that you can use multiple instances of the `STRUCTURE` keyword.  In general the number of vectors output
is equal to the number of times the `STRUCTURE` keyword is used: if there is more than one, the vectors are the components `struct-1`, `struct-2`, ...
Further note that the indices specified to the SEGMENT keywords give the indices in the list
of atoms specified by the ATOMS keyword.  `SEGMENT1=0,1,2` in the above input is measuring the distance between the instantaneous positions of atoms 13, 14 and 15
and the reference structure. All the segments must contain the same number of atoms, which must be between 2 and 30.
The coordinates in STRUCTURE are in Angstrom. TYPE selects the alignment that is used to compute the [RMSD](RMSD.md): OPTIMAL, OPTIMAL-FAST or SIMPLE.

## Periodic boundary conditions

Unless the NOPBC flag is used, each segment is made whole before computing its distance from the reference: the atoms are taken in the order given
with SEGMENT, and each atom is moved to the periodic image that is closest to the previous one, as is done by [WHOLEMOLECULES](WHOLEMOLECULES.md).

If the `ALIGN_STRANDS` flag is used the segments must contain 30 atoms, which are the two strands of a beta sheet of 15 atoms each
(the atoms N, CA, CB, C and O of three residues).  The two strands are made whole separately, as the last atom of the first strand
is not bonded to the first atom of the second strand. The second strand is then moved to the periodic image where its 7th atom
(the 22nd of the segment, which is the CA of its central residue) is closest to the 7th atom of the first strand. In other words the 22nd atom is moved to

$$
\mathbf{x}_{22} = \mathbf{x}_7 + d_{\textrm{pbc}}( \mathbf{x}_{22} - \mathbf{x}_7 )
$$

where $d_{\textrm{pbc}}$ is the distance evaluated with the minimum image convention, and the atoms from the 16th to the 30th are translated by the same vector.
This treatment is designed for making sure that the two strands of the [ANTIBETARMSD](ANTIBETARMSD.md) and [PARABETARMSD](PARABETARMSD.md)
sheets are not broken by the periodic boundaries.

The periodic boundary conditions are not used at all if you use the NOPBC flag as shown below:

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

## Working with beta sheet like structures

As discussed in the documentation for [ANTIBETARMSD](ANTIBETARMSD.md) and [PARABETARMSD](PARABETARMSD.md) the STRANDS_CUTOFF keyword (1 nm by default) speeds up the calculation
dramatically.  The reason it gives such dramatic speed ups is illustrated by the example input, which is the one that these shortcuts used before
the number of segments was computed by a single action (the vectors below are still in the expanded input, and are computed if they are used):

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
ab_ltu: LESS_THAN ARG=ab_rmsd SWITCH={RATIONAL R_0=0.1 D_0=0.0 NN=8 MM=12 D_MAX=1.77828} MASK=ab_cut
ab_lt: CUSTOM ARG=ab_ltu,ab_cut FUNC=x*y PERIODIC=NO
ab: SUM ARG=ab_lt PERIODIC=NO
PRINT ARG=ab FILE=colvar
```

The input above illustrates the input that is used by the shortcut that computes the [ANTIBETARMSD](ANTIBETARMSD.md). You can see that the distances between the atoms on the two strands of the
beta sheet are computed in a [DISTANCE](DISTANCE.md) action before we computed the SECONDARY_STRUCTURE_RMSD values.  These distances are then transformed by a Heaviside function that is 1 if
the distance is less than 1 nm and zero otherwise.  This vector of ones and zeros is then fed into the SECONDARY_STRUCTURE_RMSD action as a mask.  By doing this we ensure that SECONDARY_STRUCTURE_RMSD
values are only computed for parts of the strand that are close together; the elements of the vector that are not computed are zero. This avoids performing unnecessarily expensive calculations and is also a reasonable thing to do as if the two strands are more than
1 nm apart the residues in question are not at all similar to a beta sheet structure.

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace secondarystructure {

typedef SecondaryStructureBase<SecondaryStructureRMSDInput> colv;
PLUMED_REGISTER_ACTION(colv,"SECONDARY_STRUCTURE_RMSD")

void SecondaryStructureRMSDInput::setReferenceStructure( const std::string& type, double bondlength, std::vector<Vector>& structure ) {
  if( !myrmsd.empty() && structure.size()!=natoms ) {
    plumed_merror("input structures have a different number of atoms");
  }
  // RMSD::set puts the center of the reference in the origin
  std::vector<double> align( structure.size(), 1.0 ), displace( structure.size(), 1.0 );
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
    ColvarOutput<double>& output ) {
  // RMSD::calculate needs a std::vector for the derivatives: one per thread is reused,
  // instead of allocating it for every segment and structure
  thread_local std::vector<Vector> myderivs;
  myderivs.resize( actiondata.natoms );
  output.values[n] = actiondata.myrmsd[n].calculate( pos, myderivs, false );

  if( noderiv ) {
    return;
  }
  Tensor v;
  v.zero();
  for(unsigned j=0; j<pos.size(); ++j) {
    output.derivs[n][j] = myderivs[j];
    v -= Tensor( pos[j], myderivs[j] );
  }
  output.virial.set( n, v );
}

}
}
