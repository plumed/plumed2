/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2023 The plumed team
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
#include "MolDataClass.h"
#include "Exception.h"
#include "Tools.h"
#include "PDB.h"
#include <cctype>

namespace PLMD {

unsigned MolDataClass::numberOfAtomsPerResidueInBackbone( const std::string& type ) {
  if( type=="protein" ) {
    return 5;
  } else if( type=="dna" ) {
    return 6;
  } else if( type=="rna" ) {
    return 6;
  } else {
    return 0;
  }
}

std::string MolDataClass::baseResidueName( const std::string& residuename ) {
  if( residuename.length()==4 && residuename.find(' ')==std::string::npos ) {
    const std::string base=residuename.substr(0,3);
    if( allowedResidue("protein",base) ) {
      return base;
    }
  }
  return residuename;
}

bool MolDataClass::allowedResidue( const std::string& type, const std::string& residuename ) {
  // the four-character variants of the amino acids written by GROMACS (HISE, LYSH, ASPH, CYS2...)
  if( type=="protein" && residuename.length()==4 && baseResidueName(residuename)!=residuename ) {
    return true;
  }
  if( type=="protein" ) {
    if(residuename=="ALA") {
      return true;
    } else if(residuename=="ARG") {
      return true;
    } else if(residuename=="ASN") {
      return true;
    } else if(residuename=="ASP") {
      return true;
    } else if(residuename=="CYS") {
      return true;
    } else if(residuename=="GLN") {
      return true;
    } else if(residuename=="GLU") {
      return true;
    } else if(residuename=="GLY") {
      return true;
    } else if(residuename=="HIS") {
      return true;
    } else if(residuename=="ILE") {
      return true;
    } else if(residuename=="LEU") {
      return true;
    } else if(residuename=="LYS") {
      return true;
    } else if(residuename=="MET") {
      return true;
    } else if(residuename=="PHE") {
      return true;
    } else if(residuename=="PRO") {
      return true;
    } else if(residuename=="SER") {
      return true;
    } else if(residuename=="THR") {
      return true;
    } else if(residuename=="TRP") {
      return true;
    } else if(residuename=="TYR") {
      return true;
    } else if(residuename=="VAL") {
      return true;
    }
// Terminal groups
    else if(residuename=="ACE") {
      return true;
    } else if(residuename=="NME") {
      return true;
    } else if(residuename=="NH2") {
      return true;
    }
// Alternative residue names in common force fields
    else if(residuename=="GLH") {
      return true;  // neutral GLU
    } else if(residuename=="ASH") {
      return true;  // neutral ASP
    } else if(residuename=="HID") {
      return true;  // HIS-D amber
    } else if(residuename=="HSD") {
      return true;  // HIS-D charmm
    } else if(residuename=="HIE") {
      return true;  // HIS-E amber
    } else if(residuename=="HSE") {
      return true;  // HIS-E charmm
    } else if(residuename=="HIP") {
      return true;  // HIS-P amber
    } else if(residuename=="HSP") {
      return true;  // HIS-P charmm
    }
// Weird amino acids
    else if(residuename=="NLE") {
      return true;
    } else if(residuename=="SFO") {
      return true;
    } else {
      return false;
    }
  } else if( type=="dna" ) {
    if(residuename=="A") {
      return true;
    } else if(residuename=="A5") {
      return true;
    } else if(residuename=="A3") {
      return true;
    } else if(residuename=="AN") {
      return true;
    } else if(residuename=="G") {
      return true;
    } else if(residuename=="G5") {
      return true;
    } else if(residuename=="G3") {
      return true;
    } else if(residuename=="GN") {
      return true;
    } else if(residuename=="T") {
      return true;
    } else if(residuename=="T5") {
      return true;
    } else if(residuename=="T3") {
      return true;
    } else if(residuename=="TN") {
      return true;
    } else if(residuename=="C") {
      return true;
    } else if(residuename=="C5") {
      return true;
    } else if(residuename=="C3") {
      return true;
    } else if(residuename=="CN") {
      return true;
    } else if(residuename=="DA") {
      return true;
    } else if(residuename=="DA5") {
      return true;
    } else if(residuename=="DA3") {
      return true;
    } else if(residuename=="DAN") {
      return true;
    } else if(residuename=="DG") {
      return true;
    } else if(residuename=="DG5") {
      return true;
    } else if(residuename=="DG3") {
      return true;
    } else if(residuename=="DGN") {
      return true;
    } else if(residuename=="DT") {
      return true;
    } else if(residuename=="DT5") {
      return true;
    } else if(residuename=="DT3") {
      return true;
    } else if(residuename=="DTN") {
      return true;
    } else if(residuename=="DC") {
      return true;
    } else if(residuename=="DC5") {
      return true;
    } else if(residuename=="DC3") {
      return true;
    } else if(residuename=="DCN") {
      return true;
    } else {
      return false;
    }
  } else if( type=="rna" ) {
    if(residuename=="A") {
      return true;
    } else if(residuename=="A5") {
      return true;
    } else if(residuename=="A3") {
      return true;
    } else if(residuename=="AN") {
      return true;
    } else if(residuename=="G") {
      return true;
    } else if(residuename=="G5") {
      return true;
    } else if(residuename=="G3") {
      return true;
    } else if(residuename=="GN") {
      return true;
    } else if(residuename=="U") {
      return true;
    } else if(residuename=="U5") {
      return true;
    } else if(residuename=="U3") {
      return true;
    } else if(residuename=="UN") {
      return true;
    } else if(residuename=="C") {
      return true;
    } else if(residuename=="C5") {
      return true;
    } else if(residuename=="C3") {
      return true;
    } else if(residuename=="CN") {
      return true;
    } else if(residuename=="RA") {
      return true;
    } else if(residuename=="RA5") {
      return true;
    } else if(residuename=="RA3") {
      return true;
    } else if(residuename=="RAN") {
      return true;
    } else if(residuename=="RG") {
      return true;
    } else if(residuename=="RG5") {
      return true;
    } else if(residuename=="RG3") {
      return true;
    } else if(residuename=="RGN") {
      return true;
    } else if(residuename=="RU") {
      return true;
    } else if(residuename=="RU5") {
      return true;
    } else if(residuename=="RU3") {
      return true;
    } else if(residuename=="RUN") {
      return true;
    } else if(residuename=="RC") {
      return true;
    } else if(residuename=="RC5") {
      return true;
    } else if(residuename=="RC3") {
      return true;
    } else if(residuename=="RCN") {
      return true;
    } else {
      return false;
    }
  } else if( type=="water" ) {
    if(residuename=="SOL") {
      return true;
    }
    if(residuename=="WAT") {
      return true;
    }
    if(residuename=="HOH") {
      return true;  // wwPDB
    }
    return false;
  } else if( type=="ion" ) {
    if(residuename=="IB+") {
      return true;
    }
    if(residuename=="CA") {
      return true;
    }
    if(residuename=="CL") {
      return true;
    }
    if(residuename=="NA") {
      return true;
    }
    if(residuename=="MG") {
      return true;
    }
    if(residuename=="K") {
      return true;
    }
    if(residuename=="RB") {
      return true;
    }
    if(residuename=="CS") {
      return true;
    }
    if(residuename=="LI") {
      return true;
    }
    if(residuename=="ZN") {
      return true;
    }
    return false;
  }
  return false;
}

void MolDataClass::getBackboneForResidue( const std::string& type, int residuenum, const std::string& chain, const PDB& mypdb, std::vector<AtomNumber>& atoms ) {
  getBackboneForResidue( type, residuenum, chain, "*", mypdb, atoms );
}

void MolDataClass::getBackboneForResidue( const std::string& type, int residuenum, const std::string& chain, const std::string& icode, const PDB& mypdb, std::vector<AtomNumber>& atoms ) {
  std::string residuename=baseResidueName( mypdb.getResidueName( residuenum, chain, icode ) );
  plumed_massert( MolDataClass::allowedResidue( type, residuename ), "residue " + residuename + " unrecognized for molecule type " + type );
  if( type=="protein" ) {
    if( residuename=="GLY") {
      atoms.resize(5);
      atoms[0]=mypdb.getNamedAtomFromResidueAndChain("N",residuenum,chain,icode);
      atoms[1]=mypdb.getNamedAtomFromResidueAndChain("CA",residuenum,chain,icode);
      atoms[2]=mypdb.getNamedAtomFromResidueAndChain("HA2",residuenum,chain,icode);
      atoms[3]=mypdb.getNamedAtomFromResidueAndChain("C",residuenum,chain,icode);
      atoms[4]=mypdb.getNamedAtomFromResidueAndChain("O",residuenum,chain,icode);
    } else if( residuename=="ACE") {
      atoms.resize(1);
      atoms[0]=mypdb.getNamedAtomFromResidueAndChain("C",residuenum,chain,icode);
    } else if( residuename=="NME"||residuename=="NH2") {
      atoms.resize(1);
      atoms[0]=mypdb.getNamedAtomFromResidueAndChain("N",residuenum,chain,icode);
    } else {
      atoms.resize(5);
      atoms[0]=mypdb.getNamedAtomFromResidueAndChain("N",residuenum,chain,icode);
      atoms[1]=mypdb.getNamedAtomFromResidueAndChain("CA",residuenum,chain,icode);
      atoms[2]=mypdb.getNamedAtomFromResidueAndChain("CB",residuenum,chain,icode);
      atoms[3]=mypdb.getNamedAtomFromResidueAndChain("C",residuenum,chain,icode);
      atoms[4]=mypdb.getNamedAtomFromResidueAndChain("O",residuenum,chain,icode);
    }
  } else if( type=="dna" || type=="rna" ) {
    atoms.resize(6);
    atoms[0]=mypdb.getNamedAtomFromResidueAndChain("P",residuenum,chain,icode);
    atoms[1]=mypdb.getNamedAtomFromResidueAndChain("O5\'",residuenum,chain,icode);
    atoms[2]=mypdb.getNamedAtomFromResidueAndChain("C5\'",residuenum,chain,icode);
    atoms[3]=mypdb.getNamedAtomFromResidueAndChain("C4\'",residuenum,chain,icode);
    atoms[4]=mypdb.getNamedAtomFromResidueAndChain("C3\'",residuenum,chain,icode);
    atoms[5]=mypdb.getNamedAtomFromResidueAndChain("O3\'",residuenum,chain,icode);
  } else {
    plumed_merror(type + " is not a valid molecule type");
  }
}

bool MolDataClass::isTerminalGroup( const std::string& type, const std::string& residuename ) {
  const std::string name=baseResidueName(residuename);
  if( type=="protein" ) {
    if( name=="ACE" ) {
      return true;
    } else if( name=="NME" ) {
      return true;
    } else if( name=="NH2" ) {
      return true;
    } else {
      return false;
    }
  } else {
    plumed_merror(type + " is not a valid molecule type");
  }
  return false;
}

void MolDataClass::specialSymbol( const std::string& type, const std::string& symbol, const PDB& mypdb, std::vector<AtomNumber>& numbers ) {
  if(type=="protein" || type=="rna" || type=="dna") {
// symbol should be something like
// phi-123 i.e. phi torsion of residue 123 of first chain
// psi-A321 i.e. psi torsion of residue 321 of chain A
// psi-4_321 is psi torsion of residue 321 of chain 4
// psi-A_321 is equivalent to psi-A321
    numbers.resize(0);

// special cases:
    if(symbol=="protein") {
      const auto & all = mypdb.getAtomNumbers();
      for(auto a : all) {
        auto resname=mypdb.getResidueName(a);
        Tools::stripLeadingAndTrailingBlanks(resname);
        if(allowedResidue("protein",resname)) {
          numbers.push_back(a);
        }
      }
      return;
    }

    if(symbol=="nucleic") {
      const auto & all = mypdb.getAtomNumbers();
      for(auto a : all) {
        auto resname=mypdb.getResidueName(a);
        Tools::stripLeadingAndTrailingBlanks(resname);
        if(allowedResidue("dna",resname) || allowedResidue("rna",resname)) {
          numbers.push_back(a);
        }
      }
      return;
    }

    if(symbol=="ions") {
      const auto & all = mypdb.getAtomNumbers();
      for(auto a : all) {
        auto resname=mypdb.getResidueName(a);
        Tools::stripLeadingAndTrailingBlanks(resname);
        if(allowedResidue("ion",resname)) {
          numbers.push_back(a);
        }
      }
      return;
    }

    if(symbol=="water") {
      const auto & all = mypdb.getAtomNumbers();
      for(auto a : all) {
        auto resname=mypdb.getResidueName(a);
        Tools::stripLeadingAndTrailingBlanks(resname);
        if(allowedResidue("water",resname)) {
          numbers.push_back(a);
        }
      }
      return;
    }

    if(symbol=="hydrogens") {
      const auto & all = mypdb.getAtomNumbers();
      for(auto a : all) {
        auto atomname=mypdb.getAtomName(a);
        Tools::stripLeadingAndTrailingBlanks(atomname);
        auto notnumber=atomname.find_first_not_of("0123456789");
        if(notnumber!=std::string::npos && atomname[notnumber]=='H') {
          numbers.push_back(a);
        }
      }
      return;
    }

    if(symbol=="nonhydrogens") {
      const auto & all = mypdb.getAtomNumbers();
      for(auto a : all) {
        auto atomname=mypdb.getAtomName(a);
        Tools::stripLeadingAndTrailingBlanks(atomname);
        auto notnumber=atomname.find_first_not_of("0123456789");
        if(!(notnumber!=std::string::npos && atomname[notnumber]=='H')) {
          numbers.push_back(a);
        }
      }
      return;
    }


    std::size_t dash=symbol.find_first_of('-');
    if(dash==std::string::npos) {
      plumed_error() << "Unrecognized symbol "<<symbol;
    }

    std::string name=symbol.substr(0,dash);
    // what follows the dash is [chain][_]number[insertion code], e.g. @phi-12, @phi-A12, @phi-A1_12, @phi-A_52A.
    // The chain is what precedes the last underscore or, without underscores, what precedes the first digit;
    // without a chain, the residue is taken from the first chain that contains it.
    // A minus sign right before the residue number is part of the number (negative residue number),
    // e.g. @phi--3, @phi-A-3 (chain A, residue -3) and @phi-A_-3
    const std::string rest=symbol.substr(dash+1);
    std::string chainid;
    std::string resid;
    const std::size_t underscore=rest.find_last_of('_');
    if(underscore != std::string::npos) {
      chainid=rest.substr(0,underscore);
      resid=rest.substr(underscore+1);
    } else {
      std::size_t firstnum=rest.find_first_of("0123456789");
      if(firstnum==std::string::npos) {
        plumed_error() << "cannot find a residue number in "<<symbol;
      }
      if(firstnum>0 && rest[firstnum-1]=='-') {
        firstnum--;
      }
      chainid=(firstnum==0 ? std::string("*") : rest.substr(0,firstnum));
      resid=rest.substr(firstnum);
    }
    // the residue number can be followed by an insertion code (@phi-A_52A)
    std::size_t endnum=(resid.length()>0 && resid[0]=='-' ? 1 : 0);
    const std::size_t startdigits=endnum;
    while(endnum<resid.length() && std::isdigit(static_cast<unsigned char>(resid[endnum]))) {
      endnum++;
    }
    int resnum=0;
    if(endnum==startdigits || !Tools::convertNoexcept(resid.substr(0,endnum),resnum)) {
      plumed_error() << "cannot find a residue number in "<<symbol;
    }
    const std::string icode=resid.substr(endnum);
    // the residues before and after, for the dihedrals that involve them
    int prevnum=0, nextnum=0;
    std::string previcode, nexticode;
    mypdb.getAdjacentResidue(resnum,chainid,icode,-1,prevnum,previcode);
    mypdb.getAdjacentResidue(resnum,chainid,icode,+1,nextnum,nexticode);
    std::string resname;
    resname= mypdb.getResidueName(resnum,chainid,icode);
    Tools::stripLeadingAndTrailingBlanks(resname);
    resname=baseResidueName(resname);
    if(allowedResidue("protein",resname)) {
      if( name=="phi" && !isTerminalGroup("protein",resname) ) {
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C",prevnum,chainid,previcode));
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("N",resnum,chainid,icode));
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("CA",resnum,chainid,icode));
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C",resnum,chainid,icode));
      } else if( name=="psi" && !isTerminalGroup("protein",resname) ) {
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("N",resnum,chainid,icode));
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("CA",resnum,chainid,icode));
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C",resnum,chainid,icode));
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("N",nextnum,chainid,nexticode));
      } else if( name=="omega" && !isTerminalGroup("protein",resname) ) {
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("CA",prevnum,chainid,previcode));
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C",prevnum,chainid,previcode));
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("N",resnum,chainid,icode));
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("CA",resnum,chainid,icode));
      } else if( name=="chi1" && !isTerminalGroup("protein",resname) ) {
        if ( resname=="GLY" || resname=="ALA" || resname=="SFO" ) {
          plumed_merror("chi-1 is not defined for ALA, GLY and SFO");
        }
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("N",resnum,chainid,icode));
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("CA",resnum,chainid,icode));
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("CB",resnum,chainid,icode));
        if(resname=="ILE"||resname=="VAL") {
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("CG1",resnum,chainid,icode));
        } else if(resname=="CYS") {
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("SG",resnum,chainid,icode));
        } else if(resname=="THR") {
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("OG1",resnum,chainid,icode));
        } else if(resname=="SER") {
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("OG",resnum,chainid,icode));
        } else {
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("CG",resnum,chainid,icode));
        }
      } else if( name=="chi2" && !isTerminalGroup("protein",resname) ) {
        if ( resname=="GLY" || resname=="ALA" || resname=="SFO" || resname=="CYS" || resname=="SER" ||
             resname=="THR" || resname=="VAL" ) {
          plumed_merror("chi-2 is not defined for ALA, GLY, CYS, SER, THR, VAL and SFO");
        }
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("CA",resnum,chainid,icode));
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("CB",resnum,chainid,icode));

        if(resname=="ILE") {
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("CG1",resnum,chainid,icode));
        } else {
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("CG",resnum,chainid,icode));
        }
        if(resname=="ASN" || resname=="ASP") {
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("OD1",resnum,chainid,icode));
        } else if(resname=="HIS") {
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("ND1",resnum,chainid,icode));
        } else if(resname=="LEU" || resname=="PHE" || resname=="TRP" || resname=="TYR") {
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("CD1",resnum,chainid,icode));
        } else if(resname=="MET") {
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("SD",resnum,chainid,icode));
        } else {
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("CD",resnum,chainid,icode));
        }
      } else if( name=="chi3" && !isTerminalGroup("protein",resname) ) {
        if (!( resname=="ARG" || resname=="GLN" || resname=="GLU" || resname=="LYS" ||
               resname=="MET" )) {
          plumed_merror("chi-3 is defined only for ARG, GLN, GLU, LYS and MET");
        }
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("CB",resnum,chainid,icode));
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("CG",resnum,chainid,icode));

        if(resname=="MET") {
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("SD",resnum,chainid,icode));
        } else {
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("CD",resnum,chainid,icode));
        }
        if(resname=="GLN" || resname=="GLU") {
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("OE1",resnum,chainid,icode));
        } else if(resname=="LYS" || resname=="MET") {
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("CE",resnum,chainid,icode));
        } else {
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("NE",resnum,chainid,icode));
        }
      } else if( name=="chi4" && !isTerminalGroup("protein",resname) ) {
        if (!( resname=="ARG" || resname=="LYS" )) {
          plumed_merror("chi-4 is defined only for ARG and LYS");
        }
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("CG",resnum,chainid,icode));
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("CD",resnum,chainid,icode));

        if(resname=="ARG") {
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("NE",resnum,chainid,icode));
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("CZ",resnum,chainid,icode));
        } else {
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("CE",resnum,chainid,icode));
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("NZ",resnum,chainid,icode));
        }
      } else if( name=="chi5" && !isTerminalGroup("protein",resname) ) {
        if (!( resname=="ARG" )) {
          plumed_merror("chi-5 is defined only for ARG");
        }
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("CD",resnum,chainid,icode));
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("NE",resnum,chainid,icode));
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("CZ",resnum,chainid,icode));
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("NH1",resnum,chainid,icode));
      } else if( name=="sidechain" && !isTerminalGroup("protein",resname) ) {
        if(resname=="GLY") {
          plumed_merror("sidechain is not defined for GLY");
        }
        std::vector<AtomNumber> tmpnum1(mypdb.getAtomsInResidue(resnum,chainid,icode));
        for(unsigned i=0; i<tmpnum1.size(); i++) {
          if(mypdb.getAtomName(tmpnum1[i])!="N"&&mypdb.getAtomName(tmpnum1[i])!="CA"&&
              mypdb.getAtomName(tmpnum1[i])!="C"&&mypdb.getAtomName(tmpnum1[i])!="O"&&
              mypdb.getAtomName(tmpnum1[i])!="HA"&&mypdb.getAtomName(tmpnum1[i])!="H"&&
              mypdb.getAtomName(tmpnum1[i])!="HN"&&mypdb.getAtomName(tmpnum1[i])!="H1"&&
              mypdb.getAtomName(tmpnum1[i])!="H2"&&mypdb.getAtomName(tmpnum1[i])!="H3"&&
              mypdb.getAtomName(tmpnum1[i])!="OXT"&&mypdb.getAtomName(tmpnum1[i])!="OT1"&&
              mypdb.getAtomName(tmpnum1[i])!="OT2"&&mypdb.getAtomName(tmpnum1[i])!="OC1"&&
              mypdb.getAtomName(tmpnum1[i])!="OC2") {
            numbers.push_back( tmpnum1[i] );
          }
        }
      } else if( name=="back" && !isTerminalGroup("protein",resname) ) {
        std::vector<AtomNumber> tmpnum1(mypdb.getAtomsInResidue(resnum,chainid,icode));
        for(unsigned i=0; i<tmpnum1.size(); i++) {
          if(mypdb.getAtomName(tmpnum1[i])=="N"||mypdb.getAtomName(tmpnum1[i])=="CA"||
              mypdb.getAtomName(tmpnum1[i])=="C"||mypdb.getAtomName(tmpnum1[i])=="O"||
              mypdb.getAtomName(tmpnum1[i])=="HA"||mypdb.getAtomName(tmpnum1[i])=="H"||
              mypdb.getAtomName(tmpnum1[i])=="HN"||mypdb.getAtomName(tmpnum1[i])=="H1"||
              mypdb.getAtomName(tmpnum1[i])=="H2"||mypdb.getAtomName(tmpnum1[i])=="H3"||
              mypdb.getAtomName(tmpnum1[i])=="OXT"||mypdb.getAtomName(tmpnum1[i])=="OT1"||
              mypdb.getAtomName(tmpnum1[i])=="OT2"||mypdb.getAtomName(tmpnum1[i])=="OC1"||
              mypdb.getAtomName(tmpnum1[i])=="OC2"||mypdb.getAtomName(tmpnum1[i])=="HA1"||
              mypdb.getAtomName(tmpnum1[i])=="HA2"||mypdb.getAtomName(tmpnum1[i])=="HA3") {
            numbers.push_back( tmpnum1[i] );
          }
        }
      } else {
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain(name,resnum,chainid,icode));
      }

    } else if( allowedResidue("rna",resname) || allowedResidue("dna",resname)) {
      std::string basetype;
      if(resname.find_first_of("A")!=std::string::npos) {
        basetype+="A";
      }
      if(resname.find_first_of("U")!=std::string::npos) {
        basetype+="U";
      }
      if(resname.find_first_of("T")!=std::string::npos) {
        basetype+="T";
      }
      if(resname.find_first_of("C")!=std::string::npos) {
        basetype+="C";
      }
      if(resname.find_first_of("G")!=std::string::npos) {
        basetype+="G";
      }
      plumed_massert(basetype.length()==1,"cannot find type of rna/dna residue "+resname+" "+basetype);
      if( name=="chi" ) {
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("O4\'",resnum,chainid,icode));
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C1\'",resnum,chainid,icode));
        if(basetype=="T" || basetype=="U" || basetype=="C") {
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("N1",resnum,chainid,icode));
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C2",resnum,chainid,icode));
        } else if(basetype=="G" || basetype=="A") {
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("N9",resnum,chainid,icode));
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C4",resnum,chainid,icode));
        } else {
          plumed_error();
        }
      } else if( name=="alpha" ) {
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("O3\'",prevnum,chainid,previcode));
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("P",resnum,chainid,icode));
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("O5\'",resnum,chainid,icode));
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C5\'",resnum,chainid,icode));
      } else if( name=="beta" ) {
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("P",resnum,chainid,icode));
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("O5\'",resnum,chainid,icode));
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C5\'",resnum,chainid,icode));
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C4\'",resnum,chainid,icode));
      } else if( name=="gamma" ) {
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("O5\'",resnum,chainid,icode));
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C5\'",resnum,chainid,icode));
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C4\'",resnum,chainid,icode));
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C3\'",resnum,chainid,icode));
      } else if( name=="delta" ) {
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C5\'",resnum,chainid,icode));
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C4\'",resnum,chainid,icode));
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C3\'",resnum,chainid,icode));
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("O3\'",resnum,chainid,icode));
      } else if( name=="epsilon" ) {
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C4\'",resnum,chainid,icode));
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C3\'",resnum,chainid,icode));
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("O3\'",resnum,chainid,icode));
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("P",nextnum,chainid,nexticode));
      } else if( name=="zeta" ) {
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C3\'",resnum,chainid,icode));
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("O3\'",resnum,chainid,icode));
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("P",nextnum,chainid,nexticode));
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("O5\'",nextnum,chainid,nexticode));
      } else if( name=="v0" ) {
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C4\'",resnum,chainid,icode));
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("O4\'",resnum,chainid,icode));
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C1\'",resnum,chainid,icode));
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C2\'",resnum,chainid,icode));
      } else if( name=="v1" ) {
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("O4\'",resnum,chainid,icode));
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C1\'",resnum,chainid,icode));
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C2\'",resnum,chainid,icode));
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C3\'",resnum,chainid,icode));
      } else if( name=="v2" ) {
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C1\'",resnum,chainid,icode));
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C2\'",resnum,chainid,icode));
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C3\'",resnum,chainid,icode));
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C4\'",resnum,chainid,icode));
      } else if( name=="v3" ) {
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C2\'",resnum,chainid,icode));
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C3\'",resnum,chainid,icode));
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C4\'",resnum,chainid,icode));
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("O4\'",resnum,chainid,icode));
      } else if( name=="v4" ) {
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C3\'",resnum,chainid,icode));
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C4\'",resnum,chainid,icode));
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("O4\'",resnum,chainid,icode));
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C1\'",resnum,chainid,icode));
      } else if( name=="back" ) {
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("P",resnum,chainid,icode));
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("O5\'",resnum,chainid,icode));
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C5\'",resnum,chainid,icode));
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C4\'",resnum,chainid,icode));
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C3\'",resnum,chainid,icode));
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("O3\'",resnum,chainid,icode));
      } else if( name=="sugar" ) {
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C4\'",resnum,chainid,icode));
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("O4\'",resnum,chainid,icode));
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C1\'",resnum,chainid,icode));
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C2\'",resnum,chainid,icode));
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C3\'",resnum,chainid,icode));
      } else if( name=="base" ) {
        if(basetype=="C") {
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("N1",resnum,chainid,icode));
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C2",resnum,chainid,icode));
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("O2",resnum,chainid,icode));
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("N3",resnum,chainid,icode));
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C4",resnum,chainid,icode));
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("N4",resnum,chainid,icode));
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C5",resnum,chainid,icode));
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C6",resnum,chainid,icode));
        } else if(basetype=="U") {
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("N1",resnum,chainid,icode));
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C2",resnum,chainid,icode));
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("O2",resnum,chainid,icode));
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("N3",resnum,chainid,icode));
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C4",resnum,chainid,icode));
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("O4",resnum,chainid,icode));
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C5",resnum,chainid,icode));
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C6",resnum,chainid,icode));
        } else if(basetype=="T") {
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("N1",resnum,chainid,icode));
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C2",resnum,chainid,icode));
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("O2",resnum,chainid,icode));
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("N3",resnum,chainid,icode));
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C4",resnum,chainid,icode));
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("O4",resnum,chainid,icode));
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C5",resnum,chainid,icode));
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C7",resnum,chainid,icode));
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C6",resnum,chainid,icode));
        } else if(basetype=="G") {
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("N9",resnum,chainid,icode));
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C4",resnum,chainid,icode));
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("N3",resnum,chainid,icode));
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C2",resnum,chainid,icode));
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("N2",resnum,chainid,icode));
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("N1",resnum,chainid,icode));
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C6",resnum,chainid,icode));
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("O6",resnum,chainid,icode));
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C5",resnum,chainid,icode));
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("N7",resnum,chainid,icode));
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C8",resnum,chainid,icode));
        } else if(basetype=="A") {
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("N9",resnum,chainid,icode));
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C4",resnum,chainid,icode));
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("N1",resnum,chainid,icode));
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C2",resnum,chainid,icode));
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("N3",resnum,chainid,icode));
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C6",resnum,chainid,icode));
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("N6",resnum,chainid,icode));
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C5",resnum,chainid,icode));
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("N7",resnum,chainid,icode));
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C8",resnum,chainid,icode));
        } else {
          plumed_error();
        }
      } else if( name=="lcs" ) {
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C2",resnum,chainid,icode));
        if(basetype=="T" || basetype=="U" || basetype=="C") {
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C4",resnum,chainid,icode));
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C6",resnum,chainid,icode));
        } else if(basetype=="G" || basetype=="A") {
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C6",resnum,chainid,icode));
          numbers.push_back(mypdb.getNamedAtomFromResidueAndChain("C4",resnum,chainid,icode));
        } else {
          plumed_error();
        }
      } else {
        numbers.push_back(mypdb.getNamedAtomFromResidueAndChain(name,resnum,chainid,icode));
      }
    } else {
      numbers.push_back(mypdb.getNamedAtomFromResidueAndChain(name,resnum,chainid,icode));
    }
  } else {
    plumed_merror(type + " is not a valid molecule type");
  }
}

}
