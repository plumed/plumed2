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
#include "tools/Matrix.h"

namespace PLMD {
namespace secondarystructure {

class SecondaryStructureDRMSDInput {
private:
/// The list of reference configurations
  std::vector<std::map<std::pair<unsigned,unsigned>, double> > drmsd_targets;
/// The general input for the secondary structure variable
public:
  std::size_t natoms;
  std::size_t nstructures;
/// Are we operating without periodic boundary conditions
  bool nopbc;
/// Variables for strands cutoff
  bool align_strands;
/// The atoms involved in each of the secondary structure segments
  Matrix<unsigned> colvar_atoms;
  static void calculateDistance( unsigned n, bool noderiv, const SecondaryStructureDRMSDInput& actiondata, const std::vector<Vector>& pos, ColvarOutput& output );
  void setReferenceStructure( std::string type, double bondlength, std::vector<Vector>& structure );
  SecondaryStructureDRMSDInput& operator=( const SecondaryStructureDRMSDInput& m ) {
    natoms = m.natoms;
    nstructures = m.nstructures;
    nopbc = m.nopbc;
    align_strands = m.align_strands;
    colvar_atoms=m.colvar_atoms;
    for(unsigned i=0; i<m.drmsd_targets.size(); ++i) {
      std::map<std::pair<unsigned,unsigned>, double> targets;
      for(const auto & it : m.drmsd_targets[i] ) {
        targets[std::make_pair(it.first.first,it.first.second)] = it.second;
      }
      drmsd_targets.push_back( targets );
    }
    return *this;
  }
};

typedef SecondaryStructureBase<SecondaryStructureDRMSDInput> colv;
PLUMED_REGISTER_ACTION(colv,"SECONDARY_STRUCTURE_DRMSD");

void SecondaryStructureDRMSDInput::setReferenceStructure( std::string type, double bondlength, std::vector<Vector>& structure ) {
  std::map<std::pair<unsigned,unsigned>, double> targets;
  for(unsigned i=0; i<structure.size()-1; ++i) {
    for(unsigned j=i+1; j<structure.size(); ++j) {
      double distance = delta( structure[i], structure[j] ).modulo();
      if(distance > bondlength) {
        targets[std::make_pair(i,j)] = distance;
      }
    }
  }
  drmsd_targets.push_back( targets );
  nstructures = drmsd_targets.size();
  natoms=structure.size();
}

void SecondaryStructureDRMSDInput::calculateDistance( unsigned n, bool noderiv, const SecondaryStructureDRMSDInput& actiondata, const std::vector<Vector>& pos, ColvarOutput& output ) {
  output.virial.set( n, Tensor(0,0,0,0,0,0,0,0,0) );
  for(unsigned i=0; i<pos.size(); ++i) {
    output.derivs[n][i] = Vector(0,0,0);
  }

  double drmsd = 0;
  Vector distance;
  for(const auto & it : actiondata.drmsd_targets[n] ) {
    const unsigned k=it.first.first;
    const unsigned j=it.first.second;

    distance=delta( pos[k], pos[j] );
    const double len = distance.modulo();
    const double diff = len - it.second;
    const double der = diff / len;
    drmsd += diff*diff;

    if( !noderiv ) {
      output.derivs[n][k] += -der*distance;
      output.derivs[n][j] += der*distance;
      output.virial.set( n, output.virial[n] - der*Tensor(distance,distance) );
    }
  }
  const double inpairs = 1./static_cast<double>(actiondata.drmsd_targets[n].size());
  output.values[n] = sqrt(inpairs*drmsd);

  if( noderiv ) {
    return ;
  }

  double scalef = inpairs / output.values[n];
  output.virial.set( n, scalef*output.virial[n] );
  for(unsigned i=0; i<pos.size(); ++i ) {
    output.derivs[n][i] *= scalef;
  }
}

}
}
