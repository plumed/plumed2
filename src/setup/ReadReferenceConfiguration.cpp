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
private:
  std::vector<std::string> read_args;
public: 
  static void registerKeywords( Keywords& keys );
  explicit ReadReferenceConfiguration(const ActionOptions&ao);
  std::string getArgName( const unsigned& k ) const ;
};

PLUMED_REGISTER_ACTION(ReadReferenceConfiguration,"READ_CONFIG")

void ReadReferenceConfiguration::registerKeywords( Keywords& keys ) {
  SetupReferenceBase::registerKeywords( keys ); 
  keys.add("compulsory","REFERENCE","a file in pdb format containing the positions of the atoms in the reference structure.");
  keys.add("compulsory","NUMBER","1","if there are multiple frames in the input file which structure would you like to read in here");
  keys.addFlag("NOALIGN",false,"when this flag is NOT present the geometric center of the molecule is calculated using the align column as weights.  This geometric center is then placed at the origin.");
  keys.add("hidden","READ_ARG","this is used by pathtool to get the arguments that must be read in");
}

ReadReferenceConfiguration::ReadReferenceConfiguration(const ActionOptions&ao):
Action(ao),
SetupReferenceBase(ao)
{
  if( getNumberOfArguments()==0 ) parseVector("READ_ARG",read_args);
  std::string reference; parse("REFERENCE",reference); 
  FILE* fp=fopen(reference.c_str(),"r");
  if(!fp) error("could not open reference file " + reference );
  unsigned number; parse("NUMBER",number); Vector center; center.zero();
  bool noalign=false; parseFlag("NOALIGN",noalign);
  for(unsigned i=0;i<number;++i) {
      PDB pdb; bool do_read=pdb.readFromFilepointer(fp,atoms.usingNaturalUnits(),0.1/atoms.getUnits().getLength()); 
      if(i==number-1) {
         if( pdb.getPositions().size()==0 && getNumberOfArguments()==0 && read_args.size()==0 ) { 
             error("found no atoms in input and names of arguments to read in were not specified in input.  Use ARG");
         }
         log.printf("  reading %dth reference structure from file %s \n", number, reference.c_str());
         log.printf("  which contains %d atoms and", pdb.getPositions().size() );
         if( read_args.size()>0 ) log.printf(" %d arguments \n", read_args.size() );
         else log.printf(" %d arguments \n", getNumberOfArguments() );
         if( pdb.getPositions().size()>0 ) {
             log.printf("  indices of atoms are : ");
             for(unsigned i=0;i<pdb.getPositions().size();++i) log.printf("%d ",pdb.getAtomNumbers()[i].serial() );
             log.printf("\n"); 

             // Compute position of center of molecule
             if( !noalign ) {
                 std::vector<double> align( pdb.getOccupancy() ); double asum=0;
                 for(unsigned i=0;i<align.size();++i) asum += align[i];
                 if( asum>epsilon ) {
                     double iasum = 1 / asum;
                     for(unsigned i=0;i<align.size();++i) align[i] *= iasum; 
                 } else {
                     double iasum = 1 / pdb.size();
                     for(unsigned i=0;i<align.size();++i) align[i] = iasum;
                 }
                 for(unsigned i=0;i<pdb.getPositions().size();++i) center += align[i]*pdb.getPositions()[i];
             }
         }
         if( getNumberOfArguments()>0 ) {
             log.printf("  labels of arguments are : ");
             for(unsigned i=0;i<getNumberOfArguments();++i) log.printf("%s ", getPntrToArgument(i)->getName().c_str() );
             log.printf("\n");
         } 
         if( read_args.size()>0 ) {
             log.printf("  labels of arguments are : ");
             for(unsigned i=0;i<read_args.size();++i) log.printf("%s ", read_args[i].c_str() );
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
                 atoms.setVatomPosition( index, pdb.getPositions()[i] - center );
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
             std::vector<unsigned> shape( 1 ); shape[0] = 0; unsigned n=0; double val;
             for(unsigned i=0;i<getNumberOfArguments();++i) shape[0] += getPntrToArgument(i)->getNumberOfValues( getLabel() );
             addValue( shape ); setNotPeriodic(); getPntrToComponent(0)->buildDataStore( getLabel() );
             for(unsigned i=0;i<getNumberOfArguments();++i) {
                 if( getPntrToArgument(i)->getRank()==0 ) {
                     if( !pdb.getArgumentValue(getPntrToArgument(i)->getName(), val)  ) error("did not find argument " + getPntrToArgument(i)->getName() + " in pdb file");
                     getPntrToComponent(0)->set( n, val ); n++;
                 } else if( getPntrToArgument(i)->getRank()==1 ) {
                     for(unsigned j=0;j<getPntrToArgument(i)->getShape()[0];++j) {
                         std::string num; Tools::convert( j+1, num );
                         if( !pdb.getArgumentValue(getPntrToArgument(i)->getName() + "." + num, val)  ) 
                             error("did not find argument " + getPntrToArgument(i)->getName() + "." + num + " in pdb file");
                         getPntrToComponent(0)->set( n, val ); n++;
                     }
                 } else if( getPntrToArgument(i)->getRank()==2 ) { 
                     for(unsigned j=0;j<getPntrToArgument(i)->getShape()[0];++j) {
                         std::string jnum; Tools::convert( j+1, jnum );
                         for(unsigned k=0;k<getPntrToArgument(k)->getShape()[2];++k) {
                             std::string knum; Tools::convert( k+1, knum );
                             if( !pdb.getArgumentValue(getPntrToArgument(i)->getName() + "." + jnum + "." + knum, val)  ) 
                                 error("did not find argument " + getPntrToArgument(i)->getName() + "." + jnum  + "." + knum + " in pdb file");
                             getPntrToComponent(0)->set( n, val ); n++;
                         }
                     }
                 } else error("cannot deal with objects with ranks greater than 2");
             }
         }
         if( read_args.size()>0 ) {
             std::vector<unsigned> shape( 1 ); shape[0] = read_args.size(); double val;
             addValue( shape ); setNotPeriodic(); getPntrToComponent(0)->buildDataStore( getLabel() );
             for(unsigned i=0;i<read_args.size();++i) {
                 if( !pdb.getArgumentValue( read_args[i], val ) ) error("did not find argument " + read_args[i] + " in pdb file");
                 getPntrToComponent(0)->set( i, val );
             }
         }
         // Set the box size if this information was read
         // Tensor box( pdb.getBoxVec() );
         // atoms.setBox( &box[0][0] );
         break;
      }
      if( !do_read ) error("not enough frames input input file " + reference );
  }
}

std::string ReadReferenceConfiguration::getArgName( const unsigned& k ) const {
  if( read_args.size()>0 ) return read_args[k];
  return SetupReferenceBase::getArgName(k);
}

}
}
