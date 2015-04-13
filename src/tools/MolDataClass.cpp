/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2015 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

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
#include "MolDataClass.h"
#include "Exception.h"
#include "Tools.h"
#include "PDB.h"

namespace PLMD {

unsigned MolDataClass::numberOfAtomsPerResidueInBackbone( const std::string& type ){
  if( type=="protein" ) return 5;
  else if( type=="dna" ) return 6;
  else if( type=="rna" ) return 6;   
  else return 0;
}

bool MolDataClass::allowedResidue( const std::string& type, const std::string& residuename ){
  if( type=="protein" ){
      if(residuename=="ALA") return true;
      else if(residuename=="ARG") return true;
      else if(residuename=="ASN") return true;
      else if(residuename=="ASP") return true;
      else if(residuename=="CYS") return true;
      else if(residuename=="GLN") return true;
      else if(residuename=="GLU") return true;
      else if(residuename=="GLY") return true;
      else if(residuename=="HIS") return true;
      else if(residuename=="ILE") return true;
      else if(residuename=="LEU") return true;
      else if(residuename=="LYS") return true;
      else if(residuename=="MET") return true;
      else if(residuename=="PHE") return true;
      else if(residuename=="PRO") return true;
      else if(residuename=="SER") return true;
      else if(residuename=="THR") return true;
      else if(residuename=="TRP") return true;
      else if(residuename=="TYR") return true;
      else if(residuename=="VAL") return true;     
// Terminal groups
      else if(residuename=="ACE") return true;
      else if(residuename=="NME") return true;
// Alternative residue names in common force fiels   
      else if(residuename=="GLH") return true; // neutral GLU
      else if(residuename=="ASH") return true; // neutral ASP
      else if(residuename=="HID") return true; // neutral HIS-D amber
      else if(residuename=="HSD") return true; // neutral HIS-D charmm
      else if(residuename=="HIE") return true; // neutral HIS-E amber
      else if(residuename=="HSE") return true; // neutral HIS-E charmm
      else if(residuename=="HIP") return true; // neutral HIS-P amber
      else if(residuename=="HSP") return true; // neutral HIS-P charmm
      else return false; 
  } else if( type=="dna" ){
      if(residuename=="DA") return true;
      else if(residuename=="DG") return true;
      else if(residuename=="DT") return true;
      else if(residuename=="DC") return true;
      else if(residuename=="DA5") return true;
      else if(residuename=="DA3") return true;
      else if(residuename=="DAN") return true;
      else if(residuename=="DG5") return true;
      else if(residuename=="DG3") return true;
      else if(residuename=="DGN") return true;
      else if(residuename=="DT5") return true;
      else if(residuename=="DT3") return true;
      else if(residuename=="DTN") return true;
      else if(residuename=="DC5") return true;
      else if(residuename=="DC3") return true;
      else if(residuename=="DCN") return true;
      else return false;
  } else if( type=="rna" ){
      if(residuename=="A") return true;
      else if(residuename=="G") return true;
      else if(residuename=="U") return true;
      else if(residuename=="C") return true;
      else if(residuename=="RA") return true;
      else if(residuename=="RA5") return true;
      else if(residuename=="RA3") return true;
      else if(residuename=="RAN") return true;
      else if(residuename=="RG") return true;
      else if(residuename=="RG5") return true;
      else if(residuename=="RG3") return true;
      else if(residuename=="RGN") return true;
      else if(residuename=="RU") return true;
      else if(residuename=="RU5") return true;
      else if(residuename=="RU3") return true;
      else if(residuename=="RUN") return true;
      else if(residuename=="RC") return true;
      else if(residuename=="RC5") return true;
      else if(residuename=="RC3") return true;
      else if(residuename=="RCN") return true;
      else return false;
  } 
  return false;
}

void MolDataClass::getBackboneForResidue( const std::string& type, const unsigned& residuenum, const PDB& mypdb, std::vector<AtomNumber>& atoms ){
  std::string residuename=mypdb.getResidueName( residuenum );
  plumed_massert( MolDataClass::allowedResidue( type, residuename ), "residue " + residuename + " unrecognized for molecule type " + type );
  if( type=="protein" ){
     if( residuename=="GLY"){
         atoms.resize(5);
         atoms[0]=mypdb.getNamedAtomFromResidue("N",residuenum);
         atoms[1]=mypdb.getNamedAtomFromResidue("CA",residuenum);
         atoms[2]=mypdb.getNamedAtomFromResidue("HA1",residuenum);
         atoms[3]=mypdb.getNamedAtomFromResidue("C",residuenum);
         atoms[4]=mypdb.getNamedAtomFromResidue("O",residuenum);
     } else if( residuename=="ACE"){
         atoms.resize(1); 
         atoms[0]=mypdb.getNamedAtomFromResidue("C",residuenum);
     } else if( residuename=="NME"){
         atoms.resize(1); 
         atoms[0]=mypdb.getNamedAtomFromResidue("N",residuenum);
     } else {
         atoms.resize(5);
         atoms[0]=mypdb.getNamedAtomFromResidue("N",residuenum);
         atoms[1]=mypdb.getNamedAtomFromResidue("CA",residuenum); 
         atoms[2]=mypdb.getNamedAtomFromResidue("CB",residuenum);
         atoms[3]=mypdb.getNamedAtomFromResidue("C",residuenum);
         atoms[4]=mypdb.getNamedAtomFromResidue("O",residuenum);
     }
  } else if( type=="dna" || type=="rna" ){
      atoms.resize(6);
      atoms[0]=mypdb.getNamedAtomFromResidue("P",residuenum);
      atoms[1]=mypdb.getNamedAtomFromResidue("O5\'",residuenum);
      atoms[2]=mypdb.getNamedAtomFromResidue("C5\'",residuenum);
      atoms[3]=mypdb.getNamedAtomFromResidue("C4\'",residuenum);
      atoms[4]=mypdb.getNamedAtomFromResidue("C3\'",residuenum);
      atoms[5]=mypdb.getNamedAtomFromResidue("O3\'",residuenum);
  } 
  else {
     plumed_merror(type + " is not a valid molecule type");
  }
}

bool MolDataClass::isTerminalGroup( const std::string& type, const std::string& residuename ){
  if( type=="protein" ){
      if( residuename=="ACE" ) return true;
      else if( residuename=="NME" ) return true;
      else return false;
  } else {
      plumed_merror(type + " is not a valid molecule type");
  }  
  return false;
}

void MolDataClass::specialSymbol( const std::string& type, const std::string& symbol, const PDB& mypdb, std::vector<AtomNumber>& numbers ){
  std::string name=symbol.substr(0,symbol.find_first_of('-'));
  if( type=="protein" ){
      if( name=="phi" ){
         std::size_t dash=symbol.find_first_of('-');
         unsigned resnum; Tools::convert( symbol.substr(dash+1), resnum );
         std::string resname = mypdb.getResidueName(resnum);
         if( !allowedResidue( type, resname ) || isTerminalGroup( type, resname ) ) return ;
         numbers.resize(4); 
         numbers[0]=mypdb.getNamedAtomFromResidue("C",resnum-1); 
         numbers[1]=mypdb.getNamedAtomFromResidue("N",resnum);
         numbers[2]=mypdb.getNamedAtomFromResidue("CA",resnum);
         numbers[3]=mypdb.getNamedAtomFromResidue("C",resnum);
      } else if( name=="psi" ){
         std::size_t dash=symbol.find_first_of('-');
         unsigned resnum; Tools::convert( symbol.substr(dash+1), resnum );
         std::string resname = mypdb.getResidueName(resnum);
         if( !allowedResidue( type, resname ) || isTerminalGroup( type, resname ) ) return ;
         numbers.resize(4); 
         numbers[0]=mypdb.getNamedAtomFromResidue("N",resnum); 
         numbers[1]=mypdb.getNamedAtomFromResidue("CA",resnum);
         numbers[2]=mypdb.getNamedAtomFromResidue("C",resnum);
         numbers[3]=mypdb.getNamedAtomFromResidue("N",resnum+1);
      } else if( name=="omega" ){
         std::size_t dash=symbol.find_first_of('-');
         unsigned resnum; Tools::convert( symbol.substr(dash+1), resnum );
         std::string resname = mypdb.getResidueName(resnum);
         if( !allowedResidue( type, resname ) || isTerminalGroup( type, resname ) ) return ;
         numbers.resize(4); 
         numbers[0]=mypdb.getNamedAtomFromResidue("CA",resnum); 
         numbers[1]=mypdb.getNamedAtomFromResidue("C",resnum);
         numbers[2]=mypdb.getNamedAtomFromResidue("N",resnum+1);
         numbers[3]=mypdb.getNamedAtomFromResidue("CA",resnum+1);
      } else if( name=="chi1" ){
         std::size_t dash=symbol.find_first_of('-');
         unsigned resnum; Tools::convert( symbol.substr(dash+1), resnum );
         std::string resname = mypdb.getResidueName(resnum);
         if( !allowedResidue( type, resname ) || isTerminalGroup( type, resname ) ) return ;
         if ( resname=="GLY" || resname=="ALA" ) plumed_merror("chi-1 is not defined for Alanine and Glycine");
         numbers.resize(4); 
         numbers[0]=mypdb.getNamedAtomFromResidue("N",resnum); 
         numbers[1]=mypdb.getNamedAtomFromResidue("CA",resnum);
         numbers[2]=mypdb.getNamedAtomFromResidue("CB",resnum);
         if(resname=="ILE"||resname=="VAL") 
           numbers[3]=mypdb.getNamedAtomFromResidue("CG1",resnum);
         else if(resname=="CYS") 
           numbers[3]=mypdb.getNamedAtomFromResidue("SG",resnum);
         else if(resname=="THR") 
           numbers[3]=mypdb.getNamedAtomFromResidue("OG1",resnum);
         else if(resname=="SER") 
           numbers[3]=mypdb.getNamedAtomFromResidue("OG",resnum);
         else  numbers[3]=mypdb.getNamedAtomFromResidue("CG",resnum);
      }
  } else if( type=="dna" || type=="rna" ){
      if( name=="chi" ){
          std::size_t dash=symbol.find_first_of('-');
          unsigned resnum; Tools::convert( symbol.substr(dash+1), resnum );
          std::string resname = mypdb.getResidueName(resnum); Tools::stripLeadingAndTrailingBlanks(resname);
          if( !allowedResidue( type, resname ) ) return ;
          if( resname=="DT" || resname=="DC" || resname=="U" || resname=="C" ){
              numbers.resize(4);
              numbers[0]=mypdb.getNamedAtomFromResidue("O4\'",resnum);
              numbers[1]=mypdb.getNamedAtomFromResidue("C1\'",resnum);
              numbers[2]=mypdb.getNamedAtomFromResidue("N1",resnum);
              numbers[3]=mypdb.getNamedAtomFromResidue("C2",resnum);
          } else if( resname=="DG" || resname=="DA" || resname=="G" || resname=="A" ){
              numbers.resize(4);
              numbers[0]=mypdb.getNamedAtomFromResidue("O4\'",resnum);
              numbers[1]=mypdb.getNamedAtomFromResidue("C1\'",resnum);
              numbers[2]=mypdb.getNamedAtomFromResidue("N1",resnum);
              numbers[3]=mypdb.getNamedAtomFromResidue("C4",resnum);
          }
      } else if( name == "alpha" ){
          std::size_t dash=symbol.find_first_of('-');
          unsigned resnum; Tools::convert( symbol.substr(dash+1), resnum );
          std::string resname = mypdb.getResidueName(resnum); Tools::stripLeadingAndTrailingBlanks(resname);
          if( !allowedResidue( type, resname ) ) return ;
          numbers.resize(4);
          numbers[0]=mypdb.getNamedAtomFromResidue("O3\'",resnum-1);
          numbers[1]=mypdb.getNamedAtomFromResidue("P",resnum);
          numbers[2]=mypdb.getNamedAtomFromResidue("O5\'",resnum);
          numbers[3]=mypdb.getNamedAtomFromResidue("C5\'",resnum);
      } else if( name == "beta" ){
          std::size_t dash=symbol.find_first_of('-');
          unsigned resnum; Tools::convert( symbol.substr(dash+1), resnum );
          std::string resname = mypdb.getResidueName(resnum); Tools::stripLeadingAndTrailingBlanks(resname);
          if( !allowedResidue( type, resname ) ) return ;
          numbers.resize(4);
          numbers[0]=mypdb.getNamedAtomFromResidue("P",resnum);
          numbers[1]=mypdb.getNamedAtomFromResidue("O5\'",resnum);
          numbers[2]=mypdb.getNamedAtomFromResidue("C5\'",resnum);
          numbers[3]=mypdb.getNamedAtomFromResidue("C4\'",resnum);
      } else if( name == "gamma" ){
          std::size_t dash=symbol.find_first_of('-');
          unsigned resnum; Tools::convert( symbol.substr(dash+1), resnum );
          std::string resname = mypdb.getResidueName(resnum); Tools::stripLeadingAndTrailingBlanks(resname);
          if( !allowedResidue( type, resname ) ) return ;
          numbers.resize(4);
          numbers[0]=mypdb.getNamedAtomFromResidue("O5\'",resnum);
          numbers[1]=mypdb.getNamedAtomFromResidue("C5\'",resnum);
          numbers[2]=mypdb.getNamedAtomFromResidue("C4\'",resnum);
          numbers[3]=mypdb.getNamedAtomFromResidue("C3\'",resnum);
      } else if( name == "delta"){
          std::size_t dash=symbol.find_first_of('-');
          unsigned resnum; Tools::convert( symbol.substr(dash+1), resnum );
          std::string resname = mypdb.getResidueName(resnum); Tools::stripLeadingAndTrailingBlanks(resname);
          if( !allowedResidue( type, resname ) ) return ;
          numbers.resize(4);
          numbers[0]=mypdb.getNamedAtomFromResidue("C5\'",resnum);
          numbers[1]=mypdb.getNamedAtomFromResidue("C4\'",resnum);
          numbers[2]=mypdb.getNamedAtomFromResidue("C3\'",resnum);
          numbers[3]=mypdb.getNamedAtomFromResidue("O3\'",resnum);
      } else if( name == "epsilon" ){
          std::size_t dash=symbol.find_first_of('-');
          unsigned resnum; Tools::convert( symbol.substr(dash+1), resnum );
          std::string resname = mypdb.getResidueName(resnum); Tools::stripLeadingAndTrailingBlanks(resname);
          if( !allowedResidue( type, resname ) ) return ;
          numbers.resize(4);
          numbers[0]=mypdb.getNamedAtomFromResidue("C4\'",resnum);
          numbers[1]=mypdb.getNamedAtomFromResidue("C3\'",resnum);
          numbers[2]=mypdb.getNamedAtomFromResidue("O3\'",resnum);
          numbers[3]=mypdb.getNamedAtomFromResidue("P",resnum+1);
      } else if( name == "zeta" ){
          std::size_t dash=symbol.find_first_of('-');
          unsigned resnum; Tools::convert( symbol.substr(dash+1), resnum );
          std::string resname = mypdb.getResidueName(resnum); Tools::stripLeadingAndTrailingBlanks(resname);
          if( !allowedResidue( type, resname ) ) return ;
          numbers.resize(4);
          numbers[0]=mypdb.getNamedAtomFromResidue("C3\'",resnum);
          numbers[1]=mypdb.getNamedAtomFromResidue("O3\'",resnum);
          numbers[2]=mypdb.getNamedAtomFromResidue("P",resnum+1);
          numbers[3]=mypdb.getNamedAtomFromResidue("O5\'",resnum+1);
      } else if( name == "v0" ){
          std::size_t dash=symbol.find_first_of('-');
          unsigned resnum; Tools::convert( symbol.substr(dash+1), resnum );
          std::string resname = mypdb.getResidueName(resnum); Tools::stripLeadingAndTrailingBlanks(resname);
          if( !allowedResidue( type, resname ) ) return ;
          numbers.resize(4);
          numbers[0]=mypdb.getNamedAtomFromResidue("C4\'",resnum);
          numbers[1]=mypdb.getNamedAtomFromResidue("O4\'",resnum);
          numbers[2]=mypdb.getNamedAtomFromResidue("C1\'",resnum);
          numbers[3]=mypdb.getNamedAtomFromResidue("C2\'",resnum);
      } else if( name == "v1" ){
          std::size_t dash=symbol.find_first_of('-');
          unsigned resnum; Tools::convert( symbol.substr(dash+1), resnum );
          std::string resname = mypdb.getResidueName(resnum); Tools::stripLeadingAndTrailingBlanks(resname);
          if( !allowedResidue( type, resname ) ) return ;
          numbers.resize(4);
          numbers[0]=mypdb.getNamedAtomFromResidue("O4\'",resnum);
          numbers[1]=mypdb.getNamedAtomFromResidue("C1\'",resnum);
          numbers[2]=mypdb.getNamedAtomFromResidue("C2\'",resnum);
          numbers[3]=mypdb.getNamedAtomFromResidue("C3\'",resnum);
      } else if( name == "v2" ){
          std::size_t dash=symbol.find_first_of('-');
          unsigned resnum; Tools::convert( symbol.substr(dash+1), resnum );
          std::string resname = mypdb.getResidueName(resnum); Tools::stripLeadingAndTrailingBlanks(resname);
          if( !allowedResidue( type, resname ) ) return ;
          numbers.resize(4);
          numbers[0]=mypdb.getNamedAtomFromResidue("C1\'",resnum);
          numbers[1]=mypdb.getNamedAtomFromResidue("C2\'",resnum);
          numbers[2]=mypdb.getNamedAtomFromResidue("C3\'",resnum);
          numbers[3]=mypdb.getNamedAtomFromResidue("C4\'",resnum);
      } else if( name == "v3" ){
          std::size_t dash=symbol.find_first_of('-');
          unsigned resnum; Tools::convert( symbol.substr(dash+1), resnum );
          std::string resname = mypdb.getResidueName(resnum); Tools::stripLeadingAndTrailingBlanks(resname);
          if( !allowedResidue( type, resname ) ) return ;
          numbers.resize(4);
          numbers[0]=mypdb.getNamedAtomFromResidue("C2\'",resnum);
          numbers[1]=mypdb.getNamedAtomFromResidue("C3\'",resnum);
          numbers[2]=mypdb.getNamedAtomFromResidue("C4\'",resnum);
          numbers[3]=mypdb.getNamedAtomFromResidue("O4\'",resnum);
      } else if( name == "v4" ){
          std::size_t dash=symbol.find_first_of('-');
          unsigned resnum; Tools::convert( symbol.substr(dash+1), resnum );
          std::string resname = mypdb.getResidueName(resnum); Tools::stripLeadingAndTrailingBlanks(resname);
          if( !allowedResidue( type, resname ) ) return ;
          numbers.resize(4);
          numbers[0]=mypdb.getNamedAtomFromResidue("C3\'",resnum);
          numbers[1]=mypdb.getNamedAtomFromResidue("C4\'",resnum);
          numbers[2]=mypdb.getNamedAtomFromResidue("O4\'",resnum);
          numbers[3]=mypdb.getNamedAtomFromResidue("C1\'",resnum);
      } else if( name == "back" ){
          std::size_t dash=symbol.find_first_of('-');
          unsigned resnum; Tools::convert( symbol.substr(dash+1), resnum );
          std::string resname = mypdb.getResidueName(resnum); Tools::stripLeadingAndTrailingBlanks(resname);
          if( !allowedResidue( type, resname ) ) return ;
          numbers.resize(6);
          numbers[0]=mypdb.getNamedAtomFromResidue("P",resnum);
          numbers[1]=mypdb.getNamedAtomFromResidue("O5\'",resnum);
          numbers[2]=mypdb.getNamedAtomFromResidue("C5\'",resnum);
          numbers[3]=mypdb.getNamedAtomFromResidue("C4\'",resnum);
          numbers[4]=mypdb.getNamedAtomFromResidue("C3\'",resnum);
          numbers[5]=mypdb.getNamedAtomFromResidue("O3\'",resnum);
      } else if( name == "sugar" ){
          std::size_t dash=symbol.find_first_of('-');
          unsigned resnum; Tools::convert( symbol.substr(dash+1), resnum );
          std::string resname = mypdb.getResidueName(resnum); Tools::stripLeadingAndTrailingBlanks(resname);
          if( !allowedResidue( type, resname ) ) return ;
          numbers.resize(5);
          numbers[0]=mypdb.getNamedAtomFromResidue("C4\'",resnum);
          numbers[1]=mypdb.getNamedAtomFromResidue("O4\'",resnum);
          numbers[2]=mypdb.getNamedAtomFromResidue("C1\'",resnum);
          numbers[3]=mypdb.getNamedAtomFromResidue("C2\'",resnum);
          numbers[4]=mypdb.getNamedAtomFromResidue("C3\'",resnum);
      } else if( name == "base" ){
          std::size_t dash=symbol.find_first_of('-');
          unsigned resnum; Tools::convert( symbol.substr(dash+1), resnum );
          std::string resname = mypdb.getResidueName(resnum); Tools::stripLeadingAndTrailingBlanks(resname);
          if( !allowedResidue( type, resname ) ) return ;
          if( resname=="DC" || resname=="C" || resname=="DCN" || resname=="DC5" || resname=="DC3" || resname=="RC5" || resname=="RC3" || resname=="RCN" ){
              numbers.resize(8);
              numbers[0]=mypdb.getNamedAtomFromResidue("N1",resnum);
              numbers[1]=mypdb.getNamedAtomFromResidue("C2",resnum);
              numbers[2]=mypdb.getNamedAtomFromResidue("O2",resnum);
              numbers[3]=mypdb.getNamedAtomFromResidue("N3",resnum);
              numbers[4]=mypdb.getNamedAtomFromResidue("C4",resnum);
              numbers[5]=mypdb.getNamedAtomFromResidue("N4",resnum);
              numbers[6]=mypdb.getNamedAtomFromResidue("C5",resnum);
              numbers[7]=mypdb.getNamedAtomFromResidue("C6",resnum);
          } else if( resname=="U" || resname=="RU" || resname=="RU5" || resname=="RU3" || resname=="RUN" ){
              numbers.resize(8);
              numbers[0]=mypdb.getNamedAtomFromResidue("N1",resnum);
              numbers[1]=mypdb.getNamedAtomFromResidue("C2",resnum);
              numbers[2]=mypdb.getNamedAtomFromResidue("O2",resnum);
              numbers[3]=mypdb.getNamedAtomFromResidue("N3",resnum);
              numbers[4]=mypdb.getNamedAtomFromResidue("C4",resnum);
              numbers[5]=mypdb.getNamedAtomFromResidue("O4",resnum);
              numbers[6]=mypdb.getNamedAtomFromResidue("C5",resnum);
              numbers[7]=mypdb.getNamedAtomFromResidue("C6",resnum);
          } else if( resname=="DT" || resname=="DT5" || resname=="DT3" || resname=="DTN" ){
              numbers.resize(9);
              numbers[0]=mypdb.getNamedAtomFromResidue("N1",resnum);
              numbers[1]=mypdb.getNamedAtomFromResidue("C2",resnum);
              numbers[2]=mypdb.getNamedAtomFromResidue("O2",resnum);
              numbers[3]=mypdb.getNamedAtomFromResidue("N3",resnum);
              numbers[4]=mypdb.getNamedAtomFromResidue("C4",resnum);
              numbers[5]=mypdb.getNamedAtomFromResidue("O4",resnum);
              numbers[6]=mypdb.getNamedAtomFromResidue("C5",resnum);
              numbers[7]=mypdb.getNamedAtomFromResidue("C7",resnum);
              numbers[8]=mypdb.getNamedAtomFromResidue("C6",resnum);
          } else if( resname=="DG" || resname=="G" || resname=="DGN" || resname=="DG5" || resname=="DG3" || resname=="RG5" || resname=="RG3" || resname=="RGN" ){
              numbers.resize(11);
              numbers[0]=mypdb.getNamedAtomFromResidue("N9",resnum);
              numbers[1]=mypdb.getNamedAtomFromResidue("C4",resnum);
              numbers[2]=mypdb.getNamedAtomFromResidue("N3",resnum);
              numbers[3]=mypdb.getNamedAtomFromResidue("C2",resnum);
              numbers[4]=mypdb.getNamedAtomFromResidue("N2",resnum);
              numbers[5]=mypdb.getNamedAtomFromResidue("N1",resnum);
              numbers[6]=mypdb.getNamedAtomFromResidue("C6",resnum);
              numbers[7]=mypdb.getNamedAtomFromResidue("O6",resnum);
              numbers[8]=mypdb.getNamedAtomFromResidue("C5",resnum);
              numbers[9]=mypdb.getNamedAtomFromResidue("N7",resnum);
              numbers[10]=mypdb.getNamedAtomFromResidue("C8",resnum);
          }  else if( resname=="DA" || resname=="A" || resname=="DAN" || resname=="DA5" || resname=="DA3" || resname=="RA5" || resname=="RA3" || resname=="RAN" ){
              numbers.resize(10);
              numbers[0]=mypdb.getNamedAtomFromResidue("N9",resnum);
              numbers[1]=mypdb.getNamedAtomFromResidue("C4",resnum);
              numbers[2]=mypdb.getNamedAtomFromResidue("N1",resnum);
              numbers[3]=mypdb.getNamedAtomFromResidue("C2",resnum);
              numbers[4]=mypdb.getNamedAtomFromResidue("N3",resnum);
              numbers[5]=mypdb.getNamedAtomFromResidue("C6",resnum);
              numbers[6]=mypdb.getNamedAtomFromResidue("N6",resnum);
              numbers[7]=mypdb.getNamedAtomFromResidue("C5",resnum);
              numbers[8]=mypdb.getNamedAtomFromResidue("N7",resnum);
              numbers[9]=mypdb.getNamedAtomFromResidue("C8",resnum);
          }
      }
  }
  else {
      plumed_merror(type + " is not a valid molecule type"); 
  }
}  

}
