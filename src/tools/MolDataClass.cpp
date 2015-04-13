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
  if(type=="protein" || type=="rna" || type=="dna"){
    numbers.resize(0);
    std::size_t dash=symbol.find_first_of('-');
    std::string name=symbol.substr(0,dash);
    unsigned resnum; Tools::convert( symbol.substr(dash+1), resnum );
    std::string resname = mypdb.getResidueName(resnum);
    Tools::stripLeadingAndTrailingBlanks(resname);
    if(allowedResidue("protein",resname)){
      if( name=="phi" && !isTerminalGroup("protein",resname) ){
        numbers.push_back(mypdb.getNamedAtomFromResidue("C",resnum-1));
        numbers.push_back(mypdb.getNamedAtomFromResidue("N",resnum));
        numbers.push_back(mypdb.getNamedAtomFromResidue("CA",resnum));
        numbers.push_back(mypdb.getNamedAtomFromResidue("C",resnum));
      } else if( name=="psi" && !isTerminalGroup("protein",resname) ){
        numbers.push_back(mypdb.getNamedAtomFromResidue("N",resnum));
        numbers.push_back(mypdb.getNamedAtomFromResidue("CA",resnum));
        numbers.push_back(mypdb.getNamedAtomFromResidue("C",resnum));
        numbers.push_back(mypdb.getNamedAtomFromResidue("N",resnum+1));
      } else if( name=="omega" && !isTerminalGroup("protein",resname) ){
        numbers.push_back(mypdb.getNamedAtomFromResidue("CA",resnum));
        numbers.push_back(mypdb.getNamedAtomFromResidue("C",resnum));
        numbers.push_back(mypdb.getNamedAtomFromResidue("N",resnum+1));
        numbers.push_back(mypdb.getNamedAtomFromResidue("CA",resnum+1));
      } else if( name=="chi1" && !isTerminalGroup("protein",resname) ){
        if ( resname=="GLY" || resname=="ALA" ) plumed_merror("chi-1 is not defined for Alanine and Glycine");
        numbers.push_back(mypdb.getNamedAtomFromResidue("N",resnum));
        numbers.push_back(mypdb.getNamedAtomFromResidue("CA",resnum));
        numbers.push_back(mypdb.getNamedAtomFromResidue("CB",resnum));
        if(resname=="ILE"||resname=="VAL")
          numbers.push_back(mypdb.getNamedAtomFromResidue("CG1",resnum));
        else if(resname=="CYS")
          numbers.push_back(mypdb.getNamedAtomFromResidue("SG",resnum));
        else if(resname=="THR")
          numbers.push_back(mypdb.getNamedAtomFromResidue("OG1",resnum));
        else if(resname=="SER")
          numbers.push_back(mypdb.getNamedAtomFromResidue("OG",resnum));
        else
          numbers.push_back(mypdb.getNamedAtomFromResidue("CG",resnum));
      } else plumed_merror("protein name not recognized "+name);
    } else if( allowedResidue("rna",resname) || allowedResidue("dna",resname)){
      std::string basetype;
      if(resname.find_first_of("A")!=std::string::npos) basetype+="A";
      if(resname.find_first_of("U")!=std::string::npos) basetype+="U";
      if(resname.find_first_of("T")!=std::string::npos) basetype+="T";
      if(resname.find_first_of("C")!=std::string::npos) basetype+="C";
      if(resname.find_first_of("G")!=std::string::npos) basetype+="G";
      plumed_massert(basetype.length()==1,"cannot find type of rna/dna residue "+resname+" "+basetype);
      if( name=="chi" ) {
        numbers.push_back(mypdb.getNamedAtomFromResidue("O4\'",resnum));
        numbers.push_back(mypdb.getNamedAtomFromResidue("C1\'",resnum));
        numbers.push_back(mypdb.getNamedAtomFromResidue("N1",resnum));
        if(basetype=="T" || basetype=="U" || basetype=="C"){
          numbers.push_back(mypdb.getNamedAtomFromResidue("C2",resnum));
        } else if(basetype=="G" || basetype=="A"){
          numbers.push_back(mypdb.getNamedAtomFromResidue("C4",resnum));
        } else plumed_error();
      } else if( name=="alpha" ) {
        numbers.push_back(mypdb.getNamedAtomFromResidue("O3\'",resnum-1));
        numbers.push_back(mypdb.getNamedAtomFromResidue("P",resnum));
        numbers.push_back(mypdb.getNamedAtomFromResidue("O5\'",resnum));
        numbers.push_back(mypdb.getNamedAtomFromResidue("C5\'",resnum));
      } else if( name=="beta" ) {
        numbers.push_back(mypdb.getNamedAtomFromResidue("P",resnum));
        numbers.push_back(mypdb.getNamedAtomFromResidue("O5\'",resnum));
        numbers.push_back(mypdb.getNamedAtomFromResidue("C5\'",resnum));
        numbers.push_back(mypdb.getNamedAtomFromResidue("C4\'",resnum));
      } else if( name=="gamma" ) {
        numbers.push_back(mypdb.getNamedAtomFromResidue("O5\'",resnum));
        numbers.push_back(mypdb.getNamedAtomFromResidue("C5\'",resnum));
        numbers.push_back(mypdb.getNamedAtomFromResidue("C4\'",resnum));
        numbers.push_back(mypdb.getNamedAtomFromResidue("C3\'",resnum));
      } else if( name=="delta" ) {
        numbers.push_back(mypdb.getNamedAtomFromResidue("C5\'",resnum));
        numbers.push_back(mypdb.getNamedAtomFromResidue("C4\'",resnum));
        numbers.push_back(mypdb.getNamedAtomFromResidue("C3\'",resnum));
        numbers.push_back(mypdb.getNamedAtomFromResidue("O3\'",resnum));
      } else if( name=="epsilon" ) {
        numbers.push_back(mypdb.getNamedAtomFromResidue("C4\'",resnum));
        numbers.push_back(mypdb.getNamedAtomFromResidue("C3\'",resnum));
        numbers.push_back(mypdb.getNamedAtomFromResidue("O3\'",resnum));
        numbers.push_back(mypdb.getNamedAtomFromResidue("P",resnum+1));
      } else if( name=="zeta" ) {
        numbers.push_back(mypdb.getNamedAtomFromResidue("C3\'",resnum));
        numbers.push_back(mypdb.getNamedAtomFromResidue("O3\'",resnum));
        numbers.push_back(mypdb.getNamedAtomFromResidue("P",resnum+1));
        numbers.push_back(mypdb.getNamedAtomFromResidue("O5\'",resnum+1));
      } else if( name=="v0" ) {
        numbers.push_back(mypdb.getNamedAtomFromResidue("C4\'",resnum));
        numbers.push_back(mypdb.getNamedAtomFromResidue("O4\'",resnum));
        numbers.push_back(mypdb.getNamedAtomFromResidue("C1\'",resnum));
        numbers.push_back(mypdb.getNamedAtomFromResidue("C2\'",resnum));
      } else if( name=="v1" ) {
        numbers.push_back(mypdb.getNamedAtomFromResidue("O4\'",resnum));
        numbers.push_back(mypdb.getNamedAtomFromResidue("C1\'",resnum));
        numbers.push_back(mypdb.getNamedAtomFromResidue("C2\'",resnum));
        numbers.push_back(mypdb.getNamedAtomFromResidue("C3\'",resnum));
      } else if( name=="v2" ) {
        numbers.push_back(mypdb.getNamedAtomFromResidue("C1\'",resnum));
        numbers.push_back(mypdb.getNamedAtomFromResidue("C2\'",resnum));
        numbers.push_back(mypdb.getNamedAtomFromResidue("C3\'",resnum));
        numbers.push_back(mypdb.getNamedAtomFromResidue("C4\'",resnum));
      } else if( name=="v3" ) {
        numbers.push_back(mypdb.getNamedAtomFromResidue("C2\'",resnum));
        numbers.push_back(mypdb.getNamedAtomFromResidue("C3\'",resnum));
        numbers.push_back(mypdb.getNamedAtomFromResidue("C4\'",resnum));
        numbers.push_back(mypdb.getNamedAtomFromResidue("O4\'",resnum));
      } else if( name=="v4" ) {
        numbers.push_back(mypdb.getNamedAtomFromResidue("C3\'",resnum));
        numbers.push_back(mypdb.getNamedAtomFromResidue("C4\'",resnum));
        numbers.push_back(mypdb.getNamedAtomFromResidue("O4\'",resnum));
        numbers.push_back(mypdb.getNamedAtomFromResidue("C1\'",resnum));
      } else if( name=="back" ) {
        numbers.push_back(mypdb.getNamedAtomFromResidue("P",resnum));
        numbers.push_back(mypdb.getNamedAtomFromResidue("O5\'",resnum));
        numbers.push_back(mypdb.getNamedAtomFromResidue("C5\'",resnum));
        numbers.push_back(mypdb.getNamedAtomFromResidue("C4\'",resnum));
        numbers.push_back(mypdb.getNamedAtomFromResidue("C3\'",resnum));
        numbers.push_back(mypdb.getNamedAtomFromResidue("O3\'",resnum));
      } else if( name=="sugar" ) {
        numbers.push_back(mypdb.getNamedAtomFromResidue("C4\'",resnum));
        numbers.push_back(mypdb.getNamedAtomFromResidue("O4\'",resnum));
        numbers.push_back(mypdb.getNamedAtomFromResidue("C1\'",resnum));
        numbers.push_back(mypdb.getNamedAtomFromResidue("C2\'",resnum));
        numbers.push_back(mypdb.getNamedAtomFromResidue("C3\'",resnum));
      } else if( name=="base" ) {
        if(basetype=="C"){
          numbers.push_back(mypdb.getNamedAtomFromResidue("N1",resnum));
          numbers.push_back(mypdb.getNamedAtomFromResidue("C2",resnum));
          numbers.push_back(mypdb.getNamedAtomFromResidue("O2",resnum));
          numbers.push_back(mypdb.getNamedAtomFromResidue("N3",resnum));
          numbers.push_back(mypdb.getNamedAtomFromResidue("C4",resnum));
          numbers.push_back(mypdb.getNamedAtomFromResidue("N4",resnum));
          numbers.push_back(mypdb.getNamedAtomFromResidue("C5",resnum));
          numbers.push_back(mypdb.getNamedAtomFromResidue("C6",resnum));
        } else if(basetype=="U"){
          numbers.push_back(mypdb.getNamedAtomFromResidue("N1",resnum));
          numbers.push_back(mypdb.getNamedAtomFromResidue("C2",resnum));
          numbers.push_back(mypdb.getNamedAtomFromResidue("O2",resnum));
          numbers.push_back(mypdb.getNamedAtomFromResidue("N3",resnum));
          numbers.push_back(mypdb.getNamedAtomFromResidue("C4",resnum));
          numbers.push_back(mypdb.getNamedAtomFromResidue("O4",resnum));
          numbers.push_back(mypdb.getNamedAtomFromResidue("C5",resnum));
          numbers.push_back(mypdb.getNamedAtomFromResidue("C6",resnum));
        } else if(basetype=="T"){
          numbers.push_back(mypdb.getNamedAtomFromResidue("N1",resnum));
          numbers.push_back(mypdb.getNamedAtomFromResidue("C2",resnum));
          numbers.push_back(mypdb.getNamedAtomFromResidue("O2",resnum));
          numbers.push_back(mypdb.getNamedAtomFromResidue("N3",resnum));
          numbers.push_back(mypdb.getNamedAtomFromResidue("C4",resnum));
          numbers.push_back(mypdb.getNamedAtomFromResidue("O4",resnum));
          numbers.push_back(mypdb.getNamedAtomFromResidue("C5",resnum));
          numbers.push_back(mypdb.getNamedAtomFromResidue("C7",resnum));
          numbers.push_back(mypdb.getNamedAtomFromResidue("C6",resnum));
        } else if(basetype=="G"){
          numbers.push_back(mypdb.getNamedAtomFromResidue("N9",resnum));
          numbers.push_back(mypdb.getNamedAtomFromResidue("C4",resnum));
          numbers.push_back(mypdb.getNamedAtomFromResidue("N3",resnum));
          numbers.push_back(mypdb.getNamedAtomFromResidue("C2",resnum));
          numbers.push_back(mypdb.getNamedAtomFromResidue("N2",resnum));
          numbers.push_back(mypdb.getNamedAtomFromResidue("N1",resnum));
          numbers.push_back(mypdb.getNamedAtomFromResidue("C6",resnum));
          numbers.push_back(mypdb.getNamedAtomFromResidue("O6",resnum));
          numbers.push_back(mypdb.getNamedAtomFromResidue("C5",resnum));
          numbers.push_back(mypdb.getNamedAtomFromResidue("N7",resnum));
          numbers.push_back(mypdb.getNamedAtomFromResidue("C8",resnum));
        } else if(basetype=="A"){
          numbers.push_back(mypdb.getNamedAtomFromResidue("N9",resnum));
          numbers.push_back(mypdb.getNamedAtomFromResidue("C4",resnum));
          numbers.push_back(mypdb.getNamedAtomFromResidue("N1",resnum));
          numbers.push_back(mypdb.getNamedAtomFromResidue("C2",resnum));
          numbers.push_back(mypdb.getNamedAtomFromResidue("N3",resnum));
          numbers.push_back(mypdb.getNamedAtomFromResidue("C6",resnum));
          numbers.push_back(mypdb.getNamedAtomFromResidue("N6",resnum));
          numbers.push_back(mypdb.getNamedAtomFromResidue("C5",resnum));
          numbers.push_back(mypdb.getNamedAtomFromResidue("N7",resnum));
          numbers.push_back(mypdb.getNamedAtomFromResidue("C8",resnum));
        } else plumed_error();
      } else plumed_merror("RNA/DNA name not recognized "+name);
    }
  }
  else {
      plumed_merror(type + " is not a valid molecule type"); 
  }
}  

}
