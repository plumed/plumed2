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
#ifdef __PLUMED_HAS_OPENACC
#define __PLUMED_USE_OPENACC 1
#endif //__PLUMED_HAS_OPENACC
#include "SecondaryStructureBase.h"
#include "core/ActionRegister.h"
#include "tools/Matrix.h"

#include<vector>

//+PLUMEDOC MCOLVAR SECONDARY_STRUCTURE_DRMSD
/*
Calclulate the DRMSD between segments of a protein and a reference structure of interest

This action is used in the shortcuts [ALPHARMSD](ALPHARMSD.md), [ANTIBETARMSD](ANTIBETARMSD.md) and [PARABETARMSD](PARABETARMSD.md).  It calculates a
vector of [DRMSD](DRMSD.md) values between a single reference multiple configurations and the instantaneous
positions of various groups of atoms.  For example, in the following input we define a single set of reference
set of coordinates for 3 atoms.

```plumed
c1: SECONDARY_STRUCTURE_DRMSD ...
    ATOMS=13-24
    BONDLENGTH=0.17 STRUCTURE1=1,0,0,0,1,0,0,0,1
    SEGMENT1=0,1,2 SEGMENT2=3,4,5
    SEGMENT3=6,7,8 SEGMENT4=9,10,11
...
PRINT ARG=c1 FILE=colvar
```

A four dimensional vector is then returned that contains the DRMSD distances between the 4 sets of 3 atoms that were specified using the `SEGMENT` keywords
and the reference coordinates.  Notice that you can use multiple instances of the `STRUCTURE` keyword.  In general the the number of vectors output
is equal to the number of times the `STRUCTURE` keyword is used.  Further note that the indices specified to the SEGMENT keywords give the indices in the list
of atoms specified by the ATOMS keyword.  `SEGMENT1=0,1,2` in the above input is measuring the distance between the instaneous positions of atoms 13, 14 and 15
and the reference structure.

## Periodic boundary conditions

If you turn off periodic boundary conditions by using the NOPBC flag as shown below:

```plumed
c1: SECONDARY_STRUCTURE_DRMSD ...
    ATOMS=13-24 NOPBC
    BONDLENGTH=0.17 STRUCTURE1=1,0,0,0,1,0,0,0,1
    SEGMENT1=0,1,2 SEGMENT2=3,4,5
    SEGMENT3=6,7,8 SEGMENT4=9,10,11
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
ab_rmsd: SECONDARY_STRUCTURE_DRMSD ...
   BONDLENGTH=0.17 ALIGN_STRANDS MASK=ab_cut
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
beta sheet are computed in a [DISTANCE](DISTANCE.md) action before we computed the SECONDARY_STRUCTURE_DRMSD values.  These distances are then transformed by a Heavyside function that is 1 if
the distance is less than 1 nm and zero otherwise.  This vector of ones and zeros is then fed into the SECONDARY_STRUCTURE_DRMSD action as a mask.  By doing this we ensure that SECONDARY_STRUCTURE_DRMSD
values are only computed for parts of the strand that are close together. This avoids performing unecessarily expensive calculations and is also a reasonable thing to do as if the two strands are more than
1 nm appart the residues in question are not at all similar to a beta sheet structure.

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace secondarystructure {

///This is a sort of compact std::vector<std::vector<T>>
template<typename T>
class flexibleMemory {
  ///Stores the data
  std::vector<T> data{};
  ///Stores the boundaries
  std::vector<size_t> sizes{};
  //helpers for openACC
  T* d;
  size_t* sz;
  void update() {
    d=data.data();
    sz=sizes.data();
  }
public:
  flexibleMemory() {
    update();
  };
  ~flexibleMemory() =default;
  flexibleMemory(const flexibleMemory& m) :
    data(m.data),
    sizes(m.sizes) {
    update();
  }
  flexibleMemory& operator=(const flexibleMemory& m) {
    if (this!=&m) {
      data=m.data;
      sizes=m.sizes;
      update();
    }
    return *this;
  }
  constexpr size_t size() const {
    return sizes.size();
  }
//Access the i-th view of the data, use this to access the elements
  constexpr View< const T> get(size_t i) const {
    //use this  at the start of your function, to not calculate again and again the initial index
    size_t init=0;
    for(size_t j=0; j<i; j++) {
      init+=sz[j];
    }
    return View<const T>(d+init,sz[i]);
  }
  void extend(const std::vector<T>& v) {
    data.insert(data.end(),v.begin(),v.end());
    sizes.push_back(v.size());
    update();
  }
  void toACCDevice()const {
#pragma acc enter data copyin(this[0:1], d[0:data.size()], sz[0:sizes.size()])
  }
  void removeFromACCDevice() const {
#pragma acc exit data delete(sz[0:sizes.size()], d[0:data.size()], this[0:1])
  }
};

class SecondaryStructureDRMSDInput {
private:
  struct pairDistance {
    double distance;
    unsigned i;
    unsigned j;
  };
/// The list of reference configurations
  flexibleMemory<pairDistance> drmsd_targets{};
  flexibleMemory<unsigned> drmsd_atoms{};
/// The general input for the secondary structure variable
public:
  size_t natoms{0};
  size_t nstructures{0};
  size_t nindices_per_task{0};
/// Are we operating without periodic boundary conditions
  bool nopbc{false};
/// Variables for strands cutoff
  bool align_strands{false};
/// The atoms involved in each of the secondary structure segments
  Matrix<unsigned> colvar_atoms{};
  static void calculateDistance( unsigned n,
                                 bool noderiv,
                                 const SecondaryStructureDRMSDInput& actiondata,
                                 View<Vector> pos,
                                 ColvarOutput& output );
  void setReferenceStructure( std::string type, double bondlength, std::vector<Vector>& structure );
  void toACCDevice()const {
#pragma acc enter data copyin(this[0:1], natoms,nstructures,nindices_per_task,nopbc,align_strands)
    drmsd_targets.toACCDevice();
    drmsd_atoms.toACCDevice();
    colvar_atoms.toACCDevice();
  }
  void removeFromACCDevice() const {
    colvar_atoms.removeFromACCDevice();
    drmsd_atoms.removeFromACCDevice();
    drmsd_targets.removeFromACCDevice();
#pragma acc exit data delete(align_strands,nopbc,nindices_per_task,nstructures,natoms,this[0:1])

  }
};

typedef SecondaryStructureBase<SecondaryStructureDRMSDInput> colv;
PLUMED_REGISTER_ACTION(colv,"SECONDARY_STRUCTURE_DRMSD")

void SecondaryStructureDRMSDInput::setReferenceStructure(
  std::string /*type*/,
  double bondlength,
  std::vector<Vector>& structure ) {
  std::vector<pairDistance> targets;
  //set is ordered and contains no duplicated data
  std::set<unsigned> atoms_targets;
  //targets.reserve( (structure.size()*(structure.size()-1))/2 );?
  for(unsigned i=0; i<structure.size()-1; ++i) {
    for(unsigned j=i+1; j<structure.size(); ++j) {
      double distance = delta( structure[i], structure[j] ).modulo();
      if(distance > bondlength) {
        targets.emplace_back(pairDistance{distance,i,j});
        atoms_targets.emplace(i);
        atoms_targets.emplace(j);
      }
    }
  }
  drmsd_targets.extend( targets );

  drmsd_atoms.extend(std::vector<unsigned> {atoms_targets.begin(),atoms_targets.end()});
  nstructures = drmsd_targets.size();
  if (natoms==0) {
    natoms = structure.size();
  } else if (natoms!=structure.size()) {
    plumed_merror("input structures have a different number of atoms");
  }

}

void SecondaryStructureDRMSDInput::calculateDistance(
  const unsigned n,
  const bool noderiv,
  const SecondaryStructureDRMSDInput& actiondata,
  const View<Vector> pos,
  ColvarOutput& output ) {
  {
    const auto targetList = actiondata.drmsd_atoms.get(n);
    const auto targetAtoms=targetList.size();
    if( !noderiv ) {
      output.virial.set( n, Tensor(0,0,0,0,0,0,0,0,0) );
      for(unsigned i=0; i<targetAtoms; ++i ) {
        output.derivs[n][targetList[i]] =Vector(0.0,0.0,0.0);
      }
    }
  }

  double drmsd = 0.0;
  Vector distance;
  Vector dd;

  const auto target = actiondata.drmsd_targets.get(n);
  const auto targetSize=target.size();

  for (unsigned i=0; i<targetSize; ++i) {
    const auto& [reference,k,j] = target[i];
    distance=delta( pos[k], pos[j] );
    const double len = distance.modulo();
    const double diff = len - reference;
    drmsd += diff*diff;

    if( !noderiv ) {
      dd=distance*(diff / len);
      output.derivs[n][k] -= dd;//-3 mul
      output.derivs[n][j] += dd;//-3mul
      output.virial.set( n, output.virial[n] - Tensor(dd,distance) );//-9 mul
    }
  }
  const double inpairs = 1./static_cast<double>(targetSize);
  output.values[n] = sqrt(inpairs*drmsd);

  if( !noderiv ) {
    double scalef = inpairs / output.values[n];

    PLMD::LoopUnroller<9>::_mul(output.virial.getView(n).data(),scalef);
    const auto targetList = actiondata.drmsd_atoms.get(n);
    for(unsigned i=0; i<actiondata.natoms; ++i ) {
      output.derivs[n][targetList[i]] *= scalef;
    }
  }
}

}
}
