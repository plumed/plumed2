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

PLUMED_REGISTER_ACTION(ReadReferenceConfiguration,"READ_CONFIG")

void ReadReferenceConfiguration::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys ); ActionAtomistic::registerKeywords( keys );
  keys.add("compulsory","REFERENCE","a file in pdb format containing the positions of the atoms in the reference structure.");
  keys.add("compulsory","NUMBER","1","if there are multiple frames in the input file which structure would you like to read in here");
}

ReadReferenceConfiguration::ReadReferenceConfiguration(const ActionOptions&ao):
Action(ao),
ActionSetup(ao),
ActionAtomistic(ao),
ActionWithValue(ao),
hasatoms(false)
{
  
  std::string reference; parse("REFERENCE",reference); 
  FILE* fp=fopen(reference.c_str(),"r");
  if(!fp) error("could not open reference file " + reference );
  unsigned number; parse("NUMBER",number);
  for(unsigned i=0;i<number;++i) {
      PDB pdb; bool do_read=pdb.readFromFilepointer(fp,atoms.usingNaturalUnits(),0.1/atoms.getUnits().getLength()); 
      if(i==number-1) {
         std::vector<std::string> remark( pdb.getRemark() ); std::vector<std::string> argnames; Tools::parseVector( remark, "ARG", argnames );
         log.printf("  reading %dth reference structure from file %s \n", number, reference.c_str());
         log.printf("  which contains %d atoms and %d arguments \n", pdb.getPositions().size(), argnames.size() ); 
         if( pdb.getPositions().size()>0 ) {
             log.printf("  indices of atoms are : ");
             for(unsigned i=0;i<pdb.getPositions().size();++i) log.printf("%d ",pdb.getAtomNumbers()[i].serial() );
             log.printf("\n"); 
         }
         if( argnames.size()>0 ) {
             log.printf("  labels of arguments are : ");
             for(unsigned i=0;i<argnames.size();++i) log.printf("%s ", argnames[i].c_str() );
             log.printf("\n");
         }
         fclose(fp);
         
         // Now make virtual atoms for all these positions and set them to the pdb positions
         unsigned natoms=pdb.getPositions().size(); 
         if( natoms>0 ) {
             hasatoms=true; myindices.resize( natoms ); 
             for(unsigned i=0;i<natoms;++i){
                 myindices[i] = pdb.getAtomNumbers()[i];
                 AtomNumber index = atoms.addVirtualAtom( this );
                 mygroup.push_back( index ); 
                 atoms.setVatomMass( index, pdb.getOccupancy()[i] );
                 atoms.setVatomCharge( index, pdb.getBeta()[i] );
                 atoms.setVatomPosition( index, pdb.getPositions()[i] );
             }
             atoms.insertGroup( getLabel(), mygroup );
             // Get the block extents 
             nblocks = pdb.getNumberOfAtomBlocks(); blocks.resize( nblocks+1 );
             if( nblocks==1 ) { 
                blocks[0]=0; blocks[1]=natoms;
             } else {
                blocks[0]=0; for(unsigned i=0; i<nblocks; ++i) blocks[i+1]=pdb.getAtomBlockEnds()[i];
             } 
         }
         if( argnames.size()>0 ) {
             std::vector<unsigned> shape( 1 ); shape[0] = argnames.size();
             addValue( shape ); setNotPeriodic(); getPntrToComponent(0)->buildDataStore( getLabel() );
             for(unsigned i=0;i<argnames.size();++i) {
                 double val; Tools::parse( remark, argnames[i], val ); getPntrToComponent(0)->set( i, val );
             }
         }
         // Set the box size if this information was read
         if( pdb.cellWasRead() ) {
             Tensor box( pdb.getBox() );
             atoms.setBox( &box[0][0] );
         }
         break;
      }
      if( !do_read ) error("not enough frames input input file " + reference );
  }
}

ReadReferenceConfiguration::~ReadReferenceConfiguration() {
  if( hasatoms ) { atoms.removeVirtualAtom( this ); atoms.removeGroup( getLabel() ); }
}

void ReadReferenceConfiguration::getNatomsAndNargs( unsigned& natoms, unsigned& nargs ) const {
  // Get the number of atoms
  natoms = mygroup.size();
  // Get the number of arguments
  nargs=0; if( getNumberOfComponents()>0 ) nargs = getPntrToOutput(0)->getNumberOfValues( getLabel() );
}

void ReadReferenceConfiguration::transferDataToPlumed( const unsigned& npos, std::vector<double>& masses, std::vector<double>& charges, 
                                                       std::vector<Vector>& positions, const std::string& argname, PlumedMain& plmd ) const {
  for(unsigned i=0;i<myindices.size();++i) {
      masses[npos + i] = atoms.getVatomMass(mygroup[i]);
      charges[npos + i] = atoms.getVatomCharge(mygroup[i]);
      positions[npos + i] = atoms.getVatomPosition(mygroup[i]);
  }
  if( getNumberOfComponents()>0 ) {
      unsigned nvals = getPntrToOutput(0)->getSize(); 
      std::vector<double> valdata( nvals );
      for(unsigned i=0;i<nvals;++i) valdata[i] = getPntrToOutput(0)->get(i);
      plmd.cmd("setValue " + argname, &valdata[0] );
  }
}

}
}
