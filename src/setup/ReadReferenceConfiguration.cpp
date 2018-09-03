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
#include "SetupReferenceBase.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "tools/PDB.h"

namespace PLMD {
namespace setup {

class ReadReferenceConfiguration : public SetupReferenceBase {
public: 
  static void registerKeywords( Keywords& keys );
  explicit ReadReferenceConfiguration(const ActionOptions&ao);
};

PLUMED_REGISTER_ACTION(ReadReferenceConfiguration,"READ_CONFIG")

void ReadReferenceConfiguration::registerKeywords( Keywords& keys ) {
  SetupReferenceBase::registerKeywords( keys ); 
  keys.add("compulsory","REFERENCE","a file in pdb format containing the positions of the atoms in the reference structure.");
  keys.add("compulsory","NUMBER","1","if there are multiple frames in the input file which structure would you like to read in here");
}

ReadReferenceConfiguration::ReadReferenceConfiguration(const ActionOptions&ao):
Action(ao),
SetupReferenceBase(ao)
{
  
  std::string reference; parse("REFERENCE",reference); 
  FILE* fp=fopen(reference.c_str(),"r");
  if(!fp) error("could not open reference file " + reference );
  unsigned number; parse("NUMBER",number);
  for(unsigned i=0;i<number;++i) {
      PDB pdb; bool do_read=pdb.readFromFilepointer(fp,atoms.usingNaturalUnits(),0.1/atoms.getUnits().getLength()); 
      if(i==number-1) {
         if( pdb.getPositions().size()==0 && getNumberOfArguments()==0 ) error("found no atoms in input and names of arguments to read in were not specified in input.  Use ARG");
         log.printf("  reading %dth reference structure from file %s \n", number, reference.c_str());
         log.printf("  which contains %d atoms and %d arguments \n", pdb.getPositions().size(), getNumberOfArguments() ); 
         if( pdb.getPositions().size()>0 ) {
             log.printf("  indices of atoms are : ");
             for(unsigned i=0;i<pdb.getPositions().size();++i) log.printf("%d ",pdb.getAtomNumbers()[i].serial() );
             log.printf("\n"); 
         }
         std::vector<std::string> remark( pdb.getRemark() );
         if( getNumberOfArguments()>0 ) {
             log.printf("  labels of arguments are : ");
             for(unsigned i=0;i<getNumberOfArguments();++i) log.printf("%s ", getPntrToArgument(i)->getName().c_str() );
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
         if( getNumberOfArguments()>0 ) {
             std::vector<unsigned> shape( 1 ); shape[0] = 0; unsigned n=0;
             for(unsigned i=0;i<getNumberOfArguments();++i) shape[0] += getPntrToArgument(i)->getNumberOfValues( getLabel() );
             addValue( shape ); setNotPeriodic(); getPntrToComponent(0)->buildDataStore( getLabel() );
             for(unsigned i=0;i<getNumberOfArguments();++i) {
                 if( getPntrToArgument(i)->getRank()==0 ) {
                     double val; Tools::parse( remark, getPntrToArgument(i)->getName(), val );
                     getPntrToComponent(0)->set( n, val ); n++;
                 } else if( getPntrToArgument(i)->getRank()==1 ) {
                     for(unsigned j=0;j<getPntrToArgument(i)->getShape()[0];++j) {
                         double val; std::string num; Tools::convert( j+1, num );
                         Tools::parse( remark, getPntrToArgument(i)->getName() + "." + num, val ); 
                         getPntrToComponent(0)->set( n, val ); n++;
                     }
                 } else if( getPntrToArgument(i)->getRank()==2 ) { 
                     for(unsigned j=0;j<getPntrToArgument(i)->getShape()[0];++j) {
                         std::string jnum; Tools::convert( j+1, jnum );
                         for(unsigned k=0;k<getPntrToArgument(k)->getShape()[2];++k) {
                             double val; std::string knum; Tools::convert( k+1, knum );
                             Tools::parse( remark, getPntrToArgument(i)->getName() + "." + jnum + "." + knum, val );
                             getPntrToComponent(0)->set( n, val ); n++;
                         }
                     }
                 } else error("cannot deal with objects with ranks greater than 2");
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

}
}
