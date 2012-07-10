#include "MolInfo.h"
#include "Atoms.h"
#include "ActionRegister.h"
#include "ActionSet.h"

namespace PLMD {

//+PLUMEDOC TOPOLOGY MOLINFO
/**

This command is used to provide information on the molecules that are present in your system.

The information on the molecules in your system can either be provided in the form of a pdb file
or as a set of lists of atoms that describe the various chains in your system. If a pdb file
is used plumed the MOLINFO command will endeavor to recognize the various chains and residues that
make up the molecules in your system using the chainIDs and resnumbers from the pdb file. You can
then use this information in later commands to specify atom lists in terms residues.  For example
using this command you can find the backbone atoms in your structure automatically. 

Please be aware that the pdb parser in plumed is far from perfect. You should thus check the log file
and examine what plumed is actually doing whenenver you use the MOLINFO action.

\bug At the moment all atoms named HA1 are treated as if they are CB atoms.  This makes it possible to deal
with GLY residues in colvars like \ref ALPHARMSD. 

\par Examples

In the following example the MOLINFO command is used to provide the information on which atoms
are in the backbone of a protein to the ALPHARMSD CV.

\verbatim
MOLINFO STRUCTURE=reference.pdb
ALPHARMSD BACKBONE=all TYPE=DRMSD LESS_THAN=(SPLINE R_0=0.08 NN=8 MM=12) LABEL=a 
\endverbatim
(see also \ref ALPHARMSD)

*/
//+ENDPLUMEDOC


PLUMED_REGISTER_ACTION(MolInfo,"MOLINFO")

void MolInfo::registerKeywords( Keywords& keys ){
  ActionSetup::registerKeywords(keys);
  keys.add("compulsory","STRUCTURE","a file in pdb format containing a reference structure. "
                                    "This is used to defines the atoms in the various residues, chains, etc . " + PDB::documentation() );
  keys.add("numbered","CHAIN","(for masochists without pdb files) The atoms involved in each of the chains of interest in the structure.");
}

MolInfo::MolInfo( const ActionOptions&ao ):
Action(ao),
ActionSetup(ao)
{
  std::vector<MolInfo*> moldat=plumed.getActionSet().select<MolInfo*>();
  if( moldat.size()!=0 ) error("cannot use more than one MOLINFO action in input");

  std::vector<AtomNumber> backbone;
  parseVector("CHAIN",backbone);
  if( read_backbone.size()==0 ){
      for(unsigned i=1;;++i){
          if( !parseNumberedVector("CHAIN",i,backbone) ) break;
          read_backbone.push_back(backbone);
          backbone.resize(0);
      }
  } else {
      read_backbone.push_back(backbone);
  }
  if( read_backbone.size()==0 ){
    std::string reference; parse("STRUCTURE",reference);
    pdb.read(reference,plumed.getAtoms().usingNaturalUnits(),0.1/plumed.getAtoms().getUnits().length);
    std::vector<std::string> chains; pdb.getChainNames( chains );
    log.printf("  pdb file named %s contains %d chains \n",reference.c_str(), chains.size() );
    for(unsigned i=0;i<chains.size();++i){
       unsigned start,end; std::string errmsg;
       pdb.getResidueRange( chains[i], start, end, errmsg );
       if( errmsg.length()!=0 ) error( errmsg );
       AtomNumber astart,aend; 
       pdb.getAtomRange( chains[i], astart, aend, errmsg );
       if( errmsg.length()!=0 ) error( errmsg );
       log.printf("  chain named %s contains residues %d to %d and atoms %d to %d \n",chains[i].c_str(),start,end,astart.serial(),aend.serial());
    }
  }
  pdb.renameAtoms("HA1","CB");  // This is a hack to make this work with GLY residues 
}

void MolInfo::getBackbone( std::vector<std::string>& restrings, const std::vector<std::string>& atnames, std::vector< std::vector<AtomNumber> >& backbone ){
  if( read_backbone.size()!=0 ){
      if( restrings.size()!=1 ) error("cannot interpret anything other than all for residues when using CHAIN keywords");
      if( restrings[0]!="all" ) error("cannot interpret anything other than all for residues when using CHAIN keywords");  
      backbone.resize( read_backbone.size() );
      for(unsigned i=0;i<read_backbone.size();++i){
          backbone[i].resize( read_backbone[i].size() );
          for(unsigned j=0;j<read_backbone[i].size();++j) backbone[i][j]=read_backbone[i][j];
      }
  } else {
      if( restrings.size()==1 ){
          if( restrings[0]=="all" ){
              std::vector<std::string> chains; pdb.getChainNames( chains );
              for(unsigned i=0;i<chains.size();++i){
                  unsigned r_start, r_end; std::string errmsg, mm, nn;  
                  pdb.getResidueRange( chains[i], r_start, r_end, errmsg );
                  Tools::convert(r_start,mm); Tools::convert(r_end,nn);
                  if(i==0) restrings[0] = mm + "-" + nn;
                  else restrings.push_back(  mm + "-" + nn );
              }
          }
      }
      Tools::interpretRanges(restrings);

      // Convert the list of involved residues into a list of segments of chains
      int nk, nj; std::vector< std::vector<unsigned> > segments;
      std::vector<unsigned> thissegment;
      Tools::convert(restrings[0],nk); thissegment.push_back(nk);
      for(unsigned i=1;i<restrings.size();++i){
          Tools::convert(restrings[i-1],nk);
          Tools::convert(restrings[i],nj);
          if( (nk+1)!=nj || pdb.getChainID(nk)!=pdb.getChainID(nj) ){
             segments.push_back(thissegment);
             thissegment.resize(0); 
          } 
          thissegment.push_back(nj);
      }
      segments.push_back( thissegment );

      // And now get the backbone atoms from each segment
      backbone.resize( segments.size() ); 
      std::vector<AtomNumber> atomnumbers( atnames.size() );
      for(unsigned i=0;i<segments.size();++i){
          for(unsigned j=0;j<segments[i].size();++j){
              bool terminus=( j==0 || j==segments[i].size()-1 );
              bool foundatoms=pdb.getBackbone( segments[i][j], atnames, atomnumbers );
              if( terminus && !foundatoms ){
                  std::string num; Tools::convert( segments[i][j], num );
                  warning("Assuming residue " + num + " is a terminal group");
              } else if( !foundatoms ){
                  std::string num; Tools::convert( segments[i][j], num );
                  error("Could not find required backbone atom in residue number " + num );
              } else if( foundatoms ){
                  for(unsigned k=0;k<atomnumbers.size();++k) backbone[i].push_back( atomnumbers[k] );
              }
          }
      }
  }
  for(unsigned i=0;i<backbone.size();++i){
     if( backbone[i].size()%atnames.size()!=0 ) error("number of atoms in one of the segments of backbone is not a multiple of the number of atoms in each residue");
  }
}

}
