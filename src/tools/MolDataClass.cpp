/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

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
  } else {
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
  if( type=="protein" ){
      if( symbol.find("phi")!=std::string::npos ){
         std::size_t dash=symbol.find_first_of('-');
         unsigned resnum; Tools::convert( symbol.substr(dash+1), resnum );
         std::string resname = mypdb.getResidueName(resnum);
         if( !allowedResidue( type, resname ) || isTerminalGroup( type, resname ) ) return ;
         numbers.resize(4); 
         numbers[0]=mypdb.getNamedAtomFromResidue("C",resnum-1); 
         numbers[1]=mypdb.getNamedAtomFromResidue("N",resnum);
         numbers[2]=mypdb.getNamedAtomFromResidue("CA",resnum);
         numbers[3]=mypdb.getNamedAtomFromResidue("C",resnum);
      } else if( symbol.find("psi")!=std::string::npos ){
         std::size_t dash=symbol.find_first_of('-');
         unsigned resnum; Tools::convert( symbol.substr(dash+1), resnum );
         std::string resname = mypdb.getResidueName(resnum);
         if( !allowedResidue( type, resname ) || isTerminalGroup( type, resname ) ) return ;
         numbers.resize(4); 
         numbers[0]=mypdb.getNamedAtomFromResidue("N",resnum); 
         numbers[1]=mypdb.getNamedAtomFromResidue("CA",resnum);
         numbers[2]=mypdb.getNamedAtomFromResidue("C",resnum);
         numbers[3]=mypdb.getNamedAtomFromResidue("N",resnum+1);
      }
  } else {
      plumed_merror(type + " is not a valid molecule type"); 
  }
}  

}
