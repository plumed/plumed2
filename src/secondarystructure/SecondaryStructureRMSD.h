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
#ifndef __PLUMED_secondarystructure_SecondaryStructureRMSD_h
#define __PLUMED_secondarystructure_SecondaryStructureRMSD_h
#include "SecondaryStructureBase.h"
#include "tools/RMSD.h"

#include <vector>

namespace PLMD {
namespace secondarystructure {

/// The input of SECONDARY_STRUCTURE_RMSD: the RMSD of a segment from each of the reference structures
class SecondaryStructureRMSDInput {
/// The list of reference configurations
  std::vector<RMSD> myrmsd;
public:
  static void registerKeywords( Keywords& keys ) {
    keys.setDisplayName("SECONDARY_STRUCTURE_RMSD");
  }
  static constexpr bool needsType=true;
/// The number of atoms
  std::size_t natoms{0};
/// The number of structures
  std::size_t nstructures{0};
/// The number of indices per task
  std::size_t nindices_per_task{0};
/// Are we operating without periodic boundary conditions
  bool nopbc{false};
/// Variables for strands cutoff
  bool align_strands{false};
/// The atoms involved in each of the secondary structure segments
  Matrix<unsigned> colvar_atoms;
  static void calculateDistance( unsigned n, bool noderiv, const SecondaryStructureRMSDInput& actiondata, View<Vector> pos, ColvarOutput<double>& output );
  void setReferenceStructure( const std::string& type, double bondlength, std::vector<Vector>& structure );
  SecondaryStructureRMSDInput() = default;
  SecondaryStructureRMSDInput( const SecondaryStructureRMSDInput& m ) {
    *this = m;
  }
  SecondaryStructureRMSDInput& operator=( const SecondaryStructureRMSDInput& m ) {
    if( this==&m ) {
      return *this;
    }
    natoms = m.natoms;
    nstructures = m.nstructures;
    nindices_per_task = m.nindices_per_task;
    nopbc = m.nopbc;
    align_strands = m.align_strands;
    colvar_atoms=m.colvar_atoms;
    // setReferenceStructure appends to myrmsd, so the structures that are already here have to go
    myrmsd.clear();
    for(unsigned i=0; i<m.myrmsd.size(); ++i) {
      std::vector<Vector> mystruct( m.myrmsd[i].getReference() );
      setReferenceStructure( m.myrmsd[i].getMethod(), 0.0, mystruct );
    }
    return *this;
  }
};

}
}
#endif //__PLUMED_secondarystructure_SecondaryStructureRMSD_h
