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
#include "SecondaryStructureBase.h"
#include "core/ActionRegister.h"
#include "tools/RMSD.h"

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
/// Are we operating without periodic boundary conditions
  bool nopbc;
/// Variables for strands cutoff
  bool align_strands;
/// The atoms involved in each of the secondary structure segments
  Matrix<unsigned> colvar_atoms;
  static void calculateDistance( unsigned n, bool noderiv, const SecondaryStructureRMSDInput& actiondata, const std::vector<Vector>& pos, ColvarOutput& output );
  void setReferenceStructure( std::string type, double bondlength, std::vector<Vector>& structure );
  SecondaryStructureRMSDInput& operator=( const SecondaryStructureRMSDInput& m ) {
    natoms = m.natoms;
    nstructures = m.nstructures;
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

void SecondaryStructureRMSDInput::setReferenceStructure( std::string type, double bondlength, std::vector<Vector>& structure ) {
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

void SecondaryStructureRMSDInput::calculateDistance( unsigned n, bool noderiv, const SecondaryStructureRMSDInput& actiondata, const std::vector<Vector>& pos, ColvarOutput& output ) {
  std::vector<Vector> myderivs( pos.size() );
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
