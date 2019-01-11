/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2019 The plumed team
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
#include "DRMSD.h"
#include "MetricRegister.h"

namespace PLMD {

class IntermolecularDRMSD : public DRMSD {
private:
  unsigned nblocks;
  std::vector<unsigned> blocks;
public:
  explicit IntermolecularDRMSD( const ReferenceConfigurationOptions& ro );
  void read( const PDB& pdb );
  void setup_targets();
};

PLUMED_REGISTER_METRIC(IntermolecularDRMSD,"INTER-DRMSD")

IntermolecularDRMSD::IntermolecularDRMSD( const ReferenceConfigurationOptions& ro ):
  ReferenceConfiguration( ro ),
  DRMSD( ro ),
  nblocks(0)
{
}

void IntermolecularDRMSD::read( const PDB& pdb ) {
  readAtomsFromPDB( pdb, true ); nblocks = pdb.getNumberOfAtomBlocks(); blocks.resize( nblocks+1 );
  if( nblocks==1 ) error("Trying to compute intermolecular rmsd but found no TERs in input PDB");
  blocks[0]=0; for(unsigned i=0; i<nblocks; ++i) blocks[i+1]=pdb.getAtomBlockEnds()[i];
  readBounds( pdb ); setup_targets();
}

void IntermolecularDRMSD::setup_targets() {
  plumed_massert( bounds_were_set, "I am missing a call to DRMSD::setBoundsOnDistances");

  for(unsigned i=1; i<nblocks; ++i) {
    for(unsigned j=0; j<i; ++j) {
      for(unsigned iatom=blocks[i]; iatom<blocks[i+1]; ++iatom) {
        for(unsigned jatom=blocks[j]; jatom<blocks[j+1]; ++jatom) {
          double distance = delta( getReferencePosition(iatom), getReferencePosition(jatom) ).modulo();
          if(distance < upper && distance > lower ) targets[std::make_pair(iatom,jatom)] = distance;
        }
      }
    }
  }
}

}
