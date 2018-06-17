/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2018 The plumed team
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
#include "ReadReferenceConfiguration.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "tools/PDB.h"

namespace PLMD {
namespace setup {

PLUMED_REGISTER_ACTION(ReadReferenceConfiguration,"READ_ATOMS")

void ReadReferenceConfiguration::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys ); ActionAtomistic::registerKeywords( keys );
  keys.add("compulsory","REFERENCE","a file in pdb format containing the positions of the atoms in the reference structure.");
}

ReadReferenceConfiguration::ReadReferenceConfiguration(const ActionOptions&ao):
Action(ao),
ActionSetup(ao),
ActionAtomistic(ao)
{
  std::string reference; parse("REFERENCE",reference); PDB pdb;
  if( !pdb.read(reference,atoms.usingNaturalUnits(),0.1/atoms.getUnits().getLength()) ) error("mssing input file " + reference );
  log.printf("  reading reference structure from file %s \n",reference.c_str());
  log.printf("  which contains %d atoms \n", pdb.getPositions().size() ); 
  log.printf("  indices of atoms are : ");
  for(unsigned i=0;i<pdb.getPositions().size();++i) log.printf("%d ",pdb.getAtomNumbers()[i].serial() );
  log.printf("\n"); 

  // Now make virtual atoms for all these positions and set them to the pdb positions
  unsigned natoms=pdb.getPositions().size(); 
  std::vector<AtomNumber> mygroup; myindices.resize( natoms ); 
  for(unsigned i=0;i<natoms;++i){
      myindices[i] = pdb.getAtomNumbers()[i];
      AtomNumber index = atoms.addVirtualAtom( this );
      mygroup.push_back( index ); 
      atoms.setVatomMass( index, pdb.getOccupancy()[i] );
      atoms.setVatomCharge( index, pdb.getBeta()[i] );
      atoms.setVatomPosition( index, pdb.getPositions()[i] );
  }
  atoms.insertGroup( getLabel(), mygroup );
  // Set the box size if this information was read
  if( pdb.cellWasRead() ) {
      Tensor box( pdb.getBox() );
      atoms.setBox( &box[0][0] );
  }
}

ReadReferenceConfiguration::~ReadReferenceConfiguration() {
  atoms.removeVirtualAtom( this ); atoms.removeGroup( getLabel() );
}

}
}
