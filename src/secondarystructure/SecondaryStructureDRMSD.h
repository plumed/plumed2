/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2017-2025 The plumed team
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
#ifndef __PLUMED_secondarystructure_SecondaryStructureDRMSD_h
#define __PLUMED_secondarystructure_SecondaryStructureDRMSD_h
#include "SecondaryStructureBase.h"
#include "tools/Matrix.h"

#include<vector>

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

template <typename T>
class SecondaryStructureDRMSDInput {
private:
  struct pairDistance {
    T distance;
    unsigned i;
    unsigned j;
  };
/// The list of reference configurations
  flexibleMemory<pairDistance> drmsd_targets{};
  flexibleMemory<unsigned> drmsd_atoms{};
/// The general input for the secondary structure variable
public:
  static void registerKeywords( Keywords& keys ) {
    keys.setDisplayName("SECONDARY_STRUCTURE_DRMSD");
  }
  static constexpr bool needsBondLength=true;
  typedef T precision;
  size_t natoms{0};
  size_t nstructures{0};
  size_t nindices_per_task{0};
/// Are we operating without periodic boundary conditions
  bool nopbc{false};
/// Variables for strands cutoff
  bool align_strands{false};
/// The atoms involved in each of the secondary structure segments
  Matrix<unsigned> colvar_atoms{};
  static void calculateDistance( const unsigned n,
                                 const bool noderiv,
                                 const SecondaryStructureDRMSDInput& actiondata,
                                 View<VectorT<precision>> pos,
                                 ColvarOutput<precision>& output ) {
    const auto targetList = actiondata.drmsd_atoms.get(n);
    const auto targetAtoms=targetList.size();
    if( !noderiv ) {
      output.virial.set( n, TensorT<precision>(0,0,0,0,0,0,0,0,0) );
      for(unsigned i=0; i<targetAtoms; ++i ) {
        output.derivs[n][targetList[i]] =VectorT<precision>(0.0,0.0,0.0);
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
        output.virial.set( n, output.virial[n] - TensorT<precision>(dd,distance) );//-9 mul
      }
    }
    const double inpairs = 1./static_cast<double>(targetSize);
    output.values[n] = sqrt(inpairs*drmsd);

    if( !noderiv ) {
      double scalef = inpairs / output.values[n];

      PLMD::LoopUnroller<9>::_mul(output.virial.getView(n).data(),scalef);
      for(unsigned i=0; i<actiondata.natoms; ++i ) {
        output.derivs[n][targetList[i]] *= scalef;
      }
    }
  }

  void setReferenceStructure( const std::string& /*type*/, double bondlength, std::vector<Vector>& structure ) {
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
}
}
#endif //__PLUMED_secondarystructure_SecondaryStructureDRMSD_h
