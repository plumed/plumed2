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
Calclulate the distance between segments of a protein and a reference structure of interest

\par Examples

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
  SecondaryStructureDRMSDInput() = default;
  ~SecondaryStructureDRMSDInput() = default;
  SecondaryStructureDRMSDInput (const SecondaryStructureDRMSDInput& m ):
    drmsd_targets(m.drmsd_targets),
    drmsd_atoms(m.drmsd_atoms),
    natoms(m.natoms),
    nstructures(m.nstructures),
    nopbc(m.nopbc),
    align_strands(m.align_strands),
    colvar_atoms(m.colvar_atoms) {}

  SecondaryStructureDRMSDInput& operator=( const SecondaryStructureDRMSDInput& m ) {
    if (this!=&m) {
      natoms = m.natoms;
      nstructures = m.nstructures;
      nindices_per_task = m.nindices_per_task;
      nopbc = m.nopbc;
      align_strands = m.align_strands;
      colvar_atoms=m.colvar_atoms;
      drmsd_targets = m.drmsd_targets;
      drmsd_atoms = m.drmsd_atoms;
    }
    return *this;
  }
  static void calculateDistance( unsigned n,
                                 bool noderiv,
                                 const SecondaryStructureDRMSDInput& actiondata,
                                 View<Vector> pos,
                                 ColvarOutput& output );
  void setReferenceStructure( std::string type, double bondlength, std::vector<Vector>& structure );
  void toACCDevice()const {
#pragma acc enter data copyin(this[0:1], natoms,nstructures,nopbc,align_strands)
    drmsd_targets.toACCDevice();
    drmsd_atoms.toACCDevice();
    colvar_atoms.toACCDevice();
  }
  void removeFromACCDevice() const {
    colvar_atoms.removeFromACCDevice();
    drmsd_atoms.removeFromACCDevice();
    drmsd_targets.removeFromACCDevice();
#pragma acc exit data delete(align_strands,nopbc,nstructures,natoms,this[0:1])

  }
};

typedef SecondaryStructureBase<SecondaryStructureDRMSDInput> colv;
PLUMED_REGISTER_ACTION(colv,"SECONDARY_STRUCTURE_DRMSD");

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
  const auto targetList = actiondata.drmsd_atoms.get(n);
  const auto targetAtoms=targetList.size();
  if( !noderiv ) {
    output.virial.set( n, Tensor(0,0,0,0,0,0,0,0,0) );
    for(unsigned i=0; i<targetAtoms; ++i ) {
      output.derivs[n][targetList[i]] =Vector(0.0,0.0,0.0);
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
    const auto targetAtoms=targetList.size();
    for(unsigned i=0; i<actiondata.natoms; ++i ) {
      output.derivs[n][targetList[i]] *= scalef;
    }
  }
}

}
}
