/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2017-2023 The plumed team
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
#include "OnebeadMapping.h"
#include "tools/Exception.h"
#include "tools/Log.h"
#include "tools/OpenMP.h"
#include "tools/PDB.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <map>

namespace PLMD {
namespace isdb {

namespace {
// reports an invalid input: the handler is expected to throw
[[noreturn]] void fail(const OnebeadMapping::ErrorHandler &error, const std::string &msg) {
  error(msg);
  plumed_merror(msg);
}
}

// Detect a monoatomic ion from its residue name. Detection is by RESIDUE name only.
// The protein/nucleic-acid branches are matched before this is reached.
// Accepted aliases follow the RCSB PDB het codes plus the AMBER and CHARMM ion residue names.
int OnebeadTypes::ionTypeFromResname(const std::string &resname) {
  std::string r = resname;
  r.erase(std::remove_if(r.begin(), r.end(), ::isspace), r.end());
  std::transform(r.begin(), r.end(), r.begin(), ::toupper);
  if(r=="NA"  || r=="SOD" || r=="NA+") {
    return ION_NA;
  }
  if(r=="K"   || r=="POT" || r=="K+") {
    return ION_K;
  }
  if(r=="CL"  || r=="CLA" || r=="CL-") {
    return ION_CL;
  }
  if(r=="CA"  || r=="CAL" || r=="CA2+") {
    return ION_CA;
  }
  if(r=="MG"  || r=="MG2+") {
    return ION_MG;
  }
  if(r=="ZN"  || r=="ZN2" || r=="ZN2+") {
    return ION_ZN;
  }
  if(r=="FE2") {
    return ION_FE2;
  }
  if(r=="FE3") {
    return ION_FE3;
  }
  if(r=="MN"  || r=="MN2+") {
    return ION_MN;
  }
  return -1;
}


double OnebeadTypes::ionMass(unsigned ion) {
  switch(ion) {
  case ION_NA:
    return 22.990;
  case ION_K:
    return 39.098;
  case ION_CL:
    return 35.453;
  case ION_CA:
    return 40.078;
  case ION_MG:
    return 24.305;
  case ION_ZN:
    return 65.38;
  case ION_FE2:
    return 55.845;
  case ION_FE3:
    return 55.845;
  case ION_MN:
    return 54.938;
  default:
    plumed_merror("ionMass called with an unknown ion id");
  }
}

// Maps a monosaccharide residue name onto one of the eight ONEBEAD glycan bead
// types. Both the PDB chemical-component codes and the usual force-field spellings
// are accepted. hasN tells whether the residue carries an amide nitrogen in the
// TEMPLATE: it resolves the N-acetylhexosamines when the input file names them
// after their parent hexose, which is the convention used by several builders.
// Some builders rename the amino acid that carries a glycan (GLYCAM writes NLN for
// an N-linked Asn, OLS/OLT for an O-linked Ser/Thr). Only a hydrogen is lost to the
// glycosidic bond, so the bead is the parent amino acid and is treated as such. The
// alias is applied to the name used everywhere downstream, including the LCPO key,
// so no duplicate tables are needed.
std::string OnebeadTypes::onebeadResidueAlias(const std::string &resname) {
  if(resname=="NLN") {
    return "ASN";
  }
  if(resname=="OLS") {
    return "SER";
  }
  if(resname=="OLT") {
    return "THR";
  }
  return resname;
}

// Residues that the ONEBEAD model represents with a single bead per residue.
// PDB::allowedResidue already covers all of them, CYX included; the explicit test is
// kept only so that this stays correct if the module is ever pointed at
// MolDataClass::allowedResidue, whose protein list does not contain CYX.
bool OnebeadTypes::isOnebeadProteinResidue(const PDB &pdb, const std::string &resname) {
  return pdb.allowedResidue("protein",resname) || resname=="CYX";
}

int OnebeadTypes::glycanTypeFromResname(const std::string &resname, bool hasN) {
  const std::string &r = resname;

  // explicitly N-acetylated hexosamines
  if(r=="NAG" || r=="NDG" || r=="AGLCNA" || r=="BGLCNA") {
    return GLYC_NAG;
  }
  if(r=="NGA" || r=="A2G" || r=="AGALNA" || r=="BGALNA") {
    return GLYC_NGA;
  }
  // glucose / galactose: N-acetylated or not
  if(r=="GLC" || r=="BGC" || r=="AGLC" || r=="BGLC") {
    return hasN ? GLYC_NAG : GLYC_GLC;
  }
  if(r=="GAL" || r=="GLA" || r=="AGAL" || r=="BGAL") {
    return hasN ? GLYC_NGA : GLYC_GAL;
  }
  // mannose: the two anomers have their own parameter sets
  if(r=="MAN" || r=="AMAN") {
    return GLYC_MAN;
  }
  if(r=="BMA" || r=="BMAN") {
    return GLYC_BMA;
  }
  // 6-deoxyhexoses. FUC/FUL are alpha/beta-L-fucose and FCA/FCB their D enantiomers:
  // all four have exactly the same composition and differ only in configuration, which
  // no scattering experiment resolves. Rhamnose (RAM, RM4) is deliberately NOT here:
  // it shares the formula but is a different sugar, so it belongs to the unsupported
  // list rather than being silently given fucose parameters.
  if(r=="FUC" || r=="FUL" || r=="FCA" || r=="FCB" || r=="AFUC" || r=="BFUC") {
    return GLYC_FUC;
  }
  // N-acetylneuraminic acid
  if(r=="SIA" || r=="SLB" || r=="ANE5AC" || r=="BNE5AC") {
    return GLYC_SIA;
  }
  return -1;
}

// Monosaccharides that the ONEBEAD model recognises but has no parameters for. They
// are listed separately so that the user gets a specific message instead of the
// generic "Residue not known".
bool OnebeadTypes::isUnsupportedGlycan(const std::string &r) {
  return (r=="RAM" || r=="RM4" || r=="XXR" ||               // rhamnose
          r=="XYL" || r=="XYS" || r=="XYP" || r=="LXZ" ||   // xylose
          r=="GCU" || r=="BDP" || r=="GCV" ||               // glucuronic acid
          r=="IDS" || r=="IDR" ||                           // iduronic acid
          r=="SGN" || r=="SUS" || r=="UAP" ||               // sulfated sugars
          r=="NGC" || r=="NGE" ||                           // Neu5Gc
          r=="RIB" || r=="ARA" || r=="ARB" || r=="AHR" ||   // pentoses
          r=="GLP" || r=="PA1" || r=="GCS");                // glucosamine
}

OnebeadMapping::OnebeadMapping(const PDB &pdb, const std::vector<AtomNumber> &atoms, Log &log, const ErrorHandler &error) {
  // bead of each atom (in the order of the TEMPLATE, checked below to be that of atoms)
  // and number of atoms of each bead
  std::vector<unsigned> atoms_per_bead;
  std::vector<std::string> chains;
  pdb.getChainNames( chains );

  const unsigned n_map = std::max(atoms.size(), pdb.getPositions().size());
  glyc_type_.assign(n_map, -1);
  glyc_bridge_.assign(n_map, false);
  // Atom numbers in TEMPLATE traversal order, one entry per atom, ions included.
  std::vector<AtomNumber> atom_numbers;
  atom_numbers.reserve(atoms.size());
  // Audit trail. Accepting several spellings for the same residue is only safe if the
  // user can see what each one was taken to mean, so every distinct residue name that
  // reaches the glycan or alias path is recorded with its resolution and its count and
  // printed below.
  std::map<std::string, std::pair<int, unsigned> > glyc_seen;
  std::map<std::string, std::pair<std::string, unsigned> > alias_seen;
  std::vector<unsigned> glyc_count(NONEBEAD, 0);
  unsigned glyc_nbridge = 0;

  // Heavy-atom positions of the TEMPLATE. They are used below to tell a glycosidic
  // (bridging) oxygen from a free hydroxyl, which selects the LCPO atom class of
  // every glycan oxygen. Note that the TEMPLATE is read with scale=1, so these
  // coordinates are in Angstrom, like the LCPO radii.
  const std::vector<Vector>     &tmpl_pos = pdb.getPositions();
  const std::vector<AtomNumber> &tmpl_num = pdb.getAtomNumbers();
  std::vector<bool> tmpl_heavy(tmpl_pos.size(), false);
  for(unsigned k=0; k<tmpl_num.size(); ++k) {
    std::string an = pdb.getAtomName(tmpl_num[k]);
    if(an.empty()) {
      continue;
    }
    char e = an.at(0);
    if(isdigit(e) && an.size()>1) {
      e = an.at(1);
    }
    tmpl_heavy[k] = (e!='H');
  }
  // Number of heavy atoms covalently bound to a: C-O and C-N bonds are ~1.4 Angstrom
  // while the shortest non-bonded contacts are well above 2.2 Angstrom, so 1.8 is a
  // safe threshold. Only glycan oxygens are queried, so the cost is negligible and
  // it is paid once, at setup.
  const double bond_cut = 1.8;
  auto nHeavyBonded = [&](AtomNumber a) {
    const Vector pa = pdb.getPosition(a);
    unsigned nb = 0;
    for(unsigned k=0; k<tmpl_pos.size(); ++k) {
      if(!tmpl_heavy[k] || tmpl_num[k].index()==a.index()) {
        continue;
      }
      const Vector d = tmpl_pos[k]-pa;
      if(std::fabs(d[0])>bond_cut || std::fabs(d[1])>bond_cut || std::fabs(d[2])>bond_cut) {
        continue;
      }
      if(d.modulo2() < bond_cut*bond_cut) {
        nb++;
      }
    }
    return nb;
  };
  // Same neighbourhood, but returning the template indices instead of the count.
  auto listHeavyBonded = [&](AtomNumber a) {
    const Vector pa = pdb.getPosition(a);
    std::vector<unsigned> nb;
    for(unsigned k=0; k<tmpl_pos.size(); ++k) {
      if(!tmpl_heavy[k] || tmpl_num[k].index()==a.index()) {
        continue;
      }
      const Vector d = tmpl_pos[k]-pa;
      if(std::fabs(d[0])>bond_cut || std::fabs(d[1])>bond_cut || std::fabs(d[2])>bond_cut) {
        continue;
      }
      if(d.modulo2() < bond_cut*bond_cut) {
        nb.push_back(k);
      }
    }
    return nb;
  };

  // Every per-atom array is appended here and only here, so that they cannot drift
  // out of step. Returns the running index of the atom just added.
  auto registerAtom = [&](AtomNumber a, const std::string &rname, double mass) {
    atom_name_.push_back(pdb.getAtomName(a));
    residue_name_.push_back(rname);
    atom_numbers.push_back(a);
    mass_.push_back(mass);
    return atom_name_.size()-1;
  };

  // cycle over chains
  for(unsigned i=0; i<chains.size(); ++i) {
    // cycle over the residues present in the chain
    // (this works also with negative residue numbers and with gaps in the numbering)
    for(const int res : pdb.getResiduesInChain(chains[i])) {
      std::string Rname = pdb.getResidueName(res, chains[i]);
      Rname.erase(std::remove_if(Rname.begin(), Rname.end(), ::isspace),Rname.end());
      const std::string Rraw = Rname;
      Rname = onebeadResidueAlias(Rname);
      if(Rname!=Rraw) {
        alias_seen[Rraw].first = Rname;
        alias_seen[Rraw].second++;
      }
      std::vector<AtomNumber> res_atoms = pdb.getAtomsInResidue(res, chains[i]);
      // monoatomic ion: a single bead per ion, detected by residue name. We register the
      // atom name/residue name (so that the LCPO typing can recognise and skip it) and
      // use the tabulated atomic mass. The mass is immaterial for a one-atom bead (the COM
      // is the atom itself) but is set correctly for clarity.
      int ion = ionTypeFromResname(Rname);
      if(ion >= 0) {
        if(res_atoms.size()!=1) {
          fail(error, "Monoatomic ion residue "+Rname+" of chain "+chains[i]+" in the TEMPLATE must contain exactly one atom; found "+std::to_string(res_atoms.size())+".");
        }
        const double m = ionMass(ion);
        for(unsigned a_i=0; a_i<res_atoms.size(); a_i++) {
          registerAtom(res_atoms[a_i], Rname, m);
          bead_of_.push_back(atoms_per_bead.size());
        }
        atoms_per_bead.push_back(res_atoms.size());
        continue;
      }
      // Glycans: a monosaccharide is a single bead, like an amino acid. The presence
      // of an amide nitrogen is what separates GlcNAc from Glc and GalNAc from Gal,
      // so that both the PDB chemical-component codes (NAG/NDG, NGA/A2G) and the
      // force-field conventions that keep the parent hexose name are handled.
      bool res_hasN = false;
      for(unsigned a=0; a<res_atoms.size(); a++) {
        std::string an = pdb.getAtomName(res_atoms[a]);
        if(!an.empty() && an.at(0)=='N') {
          res_hasN = true;
          break;
        }
      }
      for(unsigned a=0; a<res_atoms.size(); a++) {
        for(unsigned b=a+1; b<res_atoms.size(); b++) {
          if(pdb.getAtomName(res_atoms[a])==pdb.getAtomName(res_atoms[b])) {
            fail(error, "Atom "+pdb.getAtomName(res_atoms[a])+" appears more than once in residue "+
                 Rname+" of chain "+chains[i]+" in the TEMPLATE. This usually means the file still "
                 "contains alternate locations (altLoc) or several models. Keep one of each before "
                 "running ONEBEAD.");
          }
        }
      }
      const int glyc = glycanTypeFromResname(Rname, res_hasN);
      const bool is_onebead_res = (glyc>=0) || isOnebeadProteinResidue(pdb,Rname);
      if(glyc>=0) {
        glyc_count[glyc]++;
        glyc_seen[Rraw].first = glyc;
        glyc_seen[Rraw].second++;
      } else if(isUnsupportedGlycan(Rname)) {
        fail(error, "Residue "+Rname+" is a monosaccharide for which no ONEBEAD parameters are "
             "available. Remove it, or exclude it from ATOMS=.");
      }
      std::vector<unsigned> tmp_residue_atom;
      tmp_residue_atom.resize(3,0);
      // cycle over atoms
      for(unsigned a=0; a<res_atoms.size(); a++) {
        // operations shared among all beads
        std::string Aname=pdb.getAtomName(res_atoms[a]);
        char type;
        char first = Aname.at(0);
        // We assume that element symbol is first letter, if not a number
        if (!isdigit(first)) {
          type = first;
          // otherwise is the second
        } else {
          type = Aname.at(1);
        }
        double amass = 0.;
        if (type == 'H') {
          amass = 1.008;
        } else if(type == 'C') {
          amass = 12.011;
        } else if(type == 'N') {
          amass = 14.007;
        } else if(type == 'O') {
          amass = 15.999;
        } else if(type == 'S') {
          amass = 32.065;
        } else if(type == 'P') {
          amass = 30.974;
        } else {
          fail(error, "Unknown element in mass extraction\n");
        }
        const unsigned idx = registerAtom(res_atoms[a], Rname, amass);
        if(glyc>=0) {
          glyc_type_[idx] = glyc;
          // a glycan oxygen bound to two heavy atoms is either the ring oxygen or
          // a glycosidic bridge; both are ethers as far as LCPO is concerned
          if(type=='O' && nHeavyBonded(res_atoms[a])>1) {
            glyc_bridge_[idx] = true;
            glyc_nbridge++;
          }
          bead_of_.push_back(atoms_per_bead.size());
        } else if(isOnebeadProteinResidue(pdb,Rname)) {
          bead_of_.push_back(atoms_per_bead.size());
        } else {
          // check for nucleic acids
          // Pentose bead
          if( Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  || Aname=="O3'"  ||
              Aname=="C3'"  || Aname=="O2'"  || Aname=="C2'"  || Aname=="C1'"  || Aname=="H5'"  ||
              Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  || Aname=="H2'"  || Aname=="H2''" ||
              Aname=="H2'2" || Aname=="H1'"  || Aname=="HO5'" || Aname=="HO3'" || Aname=="HO2'" ||
              Aname=="H5'1" || Aname=="H5'2" || Aname=="HO'2" || Aname=="H2'1" || Aname=="H5T"  ||
              Aname=="H3T" ) {
            bead_of_.push_back(atoms_per_bead.size()+0);
            tmp_residue_atom[0]++;
          }
          // Nucleobase bead
          else if(Aname=="N1"  || Aname=="N2"  || Aname=="N3"  || Aname=="N4"  || Aname=="N6"  ||
                  Aname=="N7"  || Aname=="N9"  || Aname=="C2"  || Aname=="C4"  || Aname=="C5"  ||
                  Aname=="C6"  || Aname=="C7"  || Aname=="C8"  || Aname=="O2"  || Aname=="O4"  ||
                  Aname=="O6"  || Aname=="H1"  || Aname=="H2"  || Aname=="H3"  || Aname=="H5"  ||
                  Aname=="H6"  || Aname=="H8"  || Aname=="H21" || Aname=="H22" || Aname=="H41" ||
                  Aname=="H42" || Aname=="H61" || Aname=="H62" || Aname=="H71" || Aname=="H72" ||
                  Aname=="H73" ) {
            bead_of_.push_back(atoms_per_bead.size()+1);
            tmp_residue_atom[1]++;
          }
          // PO bead
          else if(Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3" || Aname=="O1P" ||
                  Aname=="O2P" || Aname=="O3P" || Aname=="HP"  || Aname=="HOP3" ) {
            bead_of_.push_back(atoms_per_bead.size()+2);
            tmp_residue_atom[2]++;
          }
          // error
          else {
            fail(error, "Atom name "+Aname+" cannot be indexed to any bead. Check the PDB.");
          }
        }
      }
      if(is_onebead_res) {
        atoms_per_bead.push_back(res_atoms.size());
      } else {
        if(tmp_residue_atom[0]==0) {
          fail(error, "Nucleotide residue "+Rname+" of chain "+chains[i]+" in the TEMPLATE has no atoms assigned to the pentose bead.");
        }
        if(tmp_residue_atom[1]==0) {
          fail(error, "Nucleotide residue "+Rname+" of chain "+chains[i]+" in the TEMPLATE has no atoms assigned to the nucleobase bead.");
        }
        atoms_per_bead.push_back(tmp_residue_atom[0]);
        atoms_per_bead.push_back(tmp_residue_atom[1]);
        if(tmp_residue_atom[2]>0) {
          atoms_per_bead.push_back(tmp_residue_atom[2]);
        }
      }
    }
  }
  unsigned glyc_tot = 0;
  for(unsigned s=GLYC_FUC; s<=GLYC_SIA; ++s) {
    glyc_tot += glyc_count[s];
  }
  // Protein-glycan junction. The side-chain atom that receives the glycan loses its
  // hydrogen and gains a heavy neighbour: ASN ND2 becomes a secondary amide, SER OG /
  // THR OG1 / TYR OH become ethers. That is the same change of LCPO class as for a
  // glycosidic oxygen, so the acceptor is flagged here and readLCPOparam() picks up
  // the "_L" entry. The bead composition is untouched: the form factors are not
  // affected, only the solvent-accessible surface.
  unsigned nlink = 0;
  if(glyc_tot>0) {
    for(unsigned idx=0; idx<glyc_type_.size() && idx<atom_name_.size(); ++idx) {
      if(glyc_type_[idx]<0) {
        continue;
      }
      // the anomeric carbon is C1 in a pyranose and C2 in a sialic acid
      const std::string &an = atom_name_[idx];
      const bool anomeric = (glyc_type_[idx]==GLYC_SIA) ? (an=="C2") : (an=="C1");
      if(!anomeric || idx>=atom_numbers.size()) {
        continue;
      }
      const std::vector<unsigned> nb = listHeavyBonded(atom_numbers[idx]);
      for(unsigned n=0; n<nb.size(); ++n) {
        // find the running index of the bonded atom
        for(unsigned j=0; j<atom_numbers.size(); ++j) {
          if(atom_numbers[j].index()!=tmpl_num[nb[n]].index()) {
            continue;
          }
          if(glyc_type_[j]<0 && j<glyc_bridge_.size() && !glyc_bridge_[j]) {
            glyc_bridge_[j] = true;
            nlink++;
          }
          break;
        }
      }
    }
  }
  if(!alias_seen.empty()) {
    log.printf("  Residue-name aliases applied:\n");
    for(std::map<std::string, std::pair<std::string, unsigned> >::const_iterator it=alias_seen.begin();
        it!=alias_seen.end(); ++it) {
      log.printf("    %-8s -> %-8s (%u residues)\n",
                 it->first.c_str(), it->second.first.c_str(), it->second.second);
    }
  }
  if(glyc_tot>0) {
    const char* gname[8] = {"Fuc","Man(a)","Man(b)","Gal","Glc","GlcNAc","GalNAc","Neu5Ac"};
    log.printf("  Glycan beads: %u (", glyc_tot);
    for(unsigned s=GLYC_FUC; s<=GLYC_SIA; ++s) {
      if(glyc_count[s]>0) {
        log.printf(" %s:%u", gname[s-GLYC_FUC], glyc_count[s]);
      }
    }
    log.printf(" ), %u bridging/ring oxygens, %u protein/nucleic linkages\n",
               glyc_nbridge, nlink);
    log.printf("  Glycan residue names as read, and the bead type each was resolved to:\n");
    for(std::map<std::string, std::pair<int, unsigned> >::const_iterator it=glyc_seen.begin();
        it!=glyc_seen.end(); ++it) {
      log.printf("    %-8s -> %-8s (%u residues)\n", it->first.c_str(),
                 gname[it->second.first-GLYC_FUC], it->second.second);
    }
  }
  // the per-atom arrays are all appended by registerAtom(), so these can only fail if
  // a future branch bypasses it
  plumed_massert(atom_numbers.size()==atom_name_.size(),
                 "OnebeadMapping: atom_numbers and AtomResidueName out of step");
  plumed_massert(atom_numbers.size()==residue_name_.size(),
                 "OnebeadMapping: AtomResidueName rows out of step");
  plumed_massert(atom_numbers.size()==mass_.size(),
                 "OnebeadMapping: atoms_masses out of step");
  plumed_massert(atom_numbers.size()==bead_of_.size(),
                 "OnebeadMapping: residue_atom out of step");

  if(atom_numbers.size()!=atoms.size()) {
    std::string msg = "The TEMPLATE describes " + std::to_string(atom_numbers.size()) +
                      " atoms but ATOMS selects " + std::to_string(atoms.size()) + ". For ONEBEAD the "
                      "TEMPLATE must contain exactly the atoms selected by ATOMS, in the same order. ";
    if(atom_numbers.size()>atoms.size()) {
      msg += "Remove from the TEMPLATE whatever is not in ATOMS (water, ions, other chains).";
    } else {
      msg += "The TEMPLATE is missing atoms that ATOMS selects.";
    }
    fail(error, msg);
  }
  for(unsigned i=0; i<atoms.size(); ++i) {
    if(atom_numbers[i].index()!=atoms[i].index()) {
      fail(error, "TEMPLATE and ATOMS disagree at position " + std::to_string(i+1) + ": the TEMPLATE has atom " +
           std::to_string(atom_numbers[i].serial()) + " (" + atom_name_[i] + " of " +
           residue_name_[i] + ") where ATOMS has atom " + std::to_string(atoms[i].serial()) +
           ". For ONEBEAD the two must list the same atoms in the same order; note that the TEMPLATE is "
           "read chain by chain and, within a chain, by ascending residue number, so the atoms of each "
           "chain must be contiguous and its residues in ascending order.");
    }
  }
  glyc_type_.resize(atoms.size());
  glyc_bridge_.resize(atoms.size());

  // the bead types (0 if a bead is never reached, as before the split from SAXS)
  bead_type_.assign(atoms_per_bead.size(), 0);
  assignBeadTypes(pdb, atoms, error);

  // the atoms of each bead, in increasing order, with their weights in the centre
  // of mass of the bead (the masses are summed in the order of the atoms)
  const unsigned nbeads = atoms_per_bead.size();
  const unsigned natoms = atoms.size();
  bead_off_.assign(nbeads+1,0);
  for(unsigned a=0; a<natoms; ++a) {
    bead_off_[bead_of_[a]+1]++;
  }
  for(unsigned r=0; r<nbeads; ++r) {
    bead_off_[r+1] += bead_off_[r];
  }
  bead_atoms_.resize(natoms);
  std::vector<unsigned> fill(bead_off_.begin(), bead_off_.end()-1);
  std::vector<double> bead_mass(nbeads,0.);
  for(unsigned a=0; a<natoms; ++a) {
    bead_atoms_[fill[bead_of_[a]]++] = a;
    bead_mass[bead_of_[a]] += mass_[a];
  }
  bead_inv_mass_.resize(nbeads);
  for(unsigned r=0; r<nbeads; ++r) {
    bead_inv_mass_[r] = 1.0/bead_mass[r];
  }
  atom_weight_.resize(natoms);
  for(unsigned a=0; a<natoms; ++a) {
    atom_weight_[a] = mass_[a]*bead_inv_mass_[bead_of_[a]];
  }
}

void OnebeadMapping::assignBeadTypes(const PDB &pdb, const std::vector<AtomNumber> &atoms, const ErrorHandler &error) {
  for(unsigned i=0; i<atoms.size(); ++i) {
    std::string Aname = pdb.getAtomName(atoms[i]);
    std::string Rname = pdb.getResidueName(atoms[i]);
    Rname.erase(std::remove_if(Rname.begin(), Rname.end(), ::isspace),Rname.end());
    Rname = onebeadResidueAlias(Rname);
    if(Rname=="ALA") {
      bead_type_[bead_of_[i]]=ALA;
    } else if(Rname=="ARG") {
      bead_type_[bead_of_[i]]=ARG;
    } else if(Rname=="ASN") {
      bead_type_[bead_of_[i]]=ASN;
    } else if(Rname=="ASP") {
      bead_type_[bead_of_[i]]=ASP;
    } else if(Rname=="CYS") {
      bead_type_[bead_of_[i]]=CYS;
    } else if(Rname=="CYX") {
      bead_type_[bead_of_[i]]=CYX;
    } else if(Rname=="GLN") {
      bead_type_[bead_of_[i]]=GLN;
    } else if(Rname=="GLU") {
      bead_type_[bead_of_[i]]=GLU;
    } else if(Rname=="GLY") {
      bead_type_[bead_of_[i]]=GLY;
    } else if(Rname=="HIS") {
      bead_type_[bead_of_[i]]=HIS;
    } else if(Rname=="HID") {
      bead_type_[bead_of_[i]]=HIS;
    } else if(Rname=="HIE") {
      bead_type_[bead_of_[i]]=HIS;
    } else if(Rname=="HIP") {
      bead_type_[bead_of_[i]]=HIP;
      // CHARMM NAMING FOR PROTONATION STATES OF HISTIDINE
    } else if(Rname=="HSD") {
      bead_type_[bead_of_[i]]=HIS;
    } else if(Rname=="HSE") {
      bead_type_[bead_of_[i]]=HIS;
    } else if(Rname=="HSP") {
      bead_type_[bead_of_[i]]=HIP;
    } else if(Rname=="ILE") {
      bead_type_[bead_of_[i]]=ILE;
    } else if(Rname=="LEU") {
      bead_type_[bead_of_[i]]=LEU;
    } else if(Rname=="LYS") {
      bead_type_[bead_of_[i]]=LYS;
    } else if(Rname=="MET") {
      bead_type_[bead_of_[i]]=MET;
    } else if(Rname=="PHE") {
      bead_type_[bead_of_[i]]=PHE;
    } else if(Rname=="PRO") {
      bead_type_[bead_of_[i]]=PRO;
    } else if(Rname=="SER") {
      bead_type_[bead_of_[i]]=SER;
    } else if(Rname=="THR") {
      bead_type_[bead_of_[i]]=THR;
    } else if(Rname=="TRP") {
      bead_type_[bead_of_[i]]=TRP;
    } else if(Rname=="TYR") {
      bead_type_[bead_of_[i]]=TYR;
    } else if(Rname=="VAL") {
      bead_type_[bead_of_[i]]=VAL;
    }
    // NUCLEIC ACIDS
    // nucleobases are not automatically populated as an additional check on the health of the PDB.
    // RNA - G
    else if(Rname=="G") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="HP" ||
          Aname=="O1P" || Aname=="O2P"  ) {
        bead_type_[bead_of_[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="O2'"  || Aname=="C2'"  ||
                Aname=="C1'"  || Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  ||
                Aname=="H3'"  || Aname=="H2'"  || Aname=="H1'"  || Aname=="H3T"  ||
                Aname=="HO5'" || Aname=="HO3'" || Aname=="HO2'" || Aname=="H5'1" ||
                Aname=="H5'2" || Aname=="HO'2" || Aname=="H2'1" || Aname=="H5T" ) {
        bead_type_[bead_of_[i]]=BB_RNA;
      } else if( Aname=="N1" || Aname=="C2"  || Aname=="N2" || Aname=="N3" ||
                 Aname=="C4" || Aname=="C5"  || Aname=="C6" || Aname=="O6" ||
                 Aname=="N7" || Aname=="C8"  || Aname=="N9" || Aname=="H1" ||
                 Aname=="H8" || Aname=="H21" || Aname=="H22" ) {
        bead_type_[bead_of_[i]]=BASE_G;
      } else {
        fail(error, "Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // RNA - G3
    } else if(Rname=="G3") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3" ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" ) {
        bead_type_[bead_of_[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'" || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'" || Aname=="O2'"  || Aname=="C2'"  ||
                Aname=="C1'"  || Aname=="H5'" || Aname=="H5''" || Aname=="H4'"  ||
                Aname=="H3'"  || Aname=="H2'" || Aname=="H1'"  || Aname=="H3T"  ||
                Aname=="H2'1" || Aname=="HO3'"|| Aname=="HO2'" || Aname=="H5'1" ||
                Aname=="H5'2" || Aname=="HO'2" ) {
        bead_type_[bead_of_[i]]=BB_RNA_3;
      } else if(Aname=="N1" || Aname=="C2"  || Aname=="N2" || Aname=="N3" ||
                Aname=="C4" || Aname=="C5"  || Aname=="C6" || Aname=="O6" ||
                Aname=="N7" || Aname=="C8"  || Aname=="N9" || Aname=="H1" ||
                Aname=="H8" || Aname=="H21" || Aname=="H22" ) {
        bead_type_[bead_of_[i]]=BASE_G;
      } else {
        fail(error, "Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // RNA - G5
    } else if(Rname=="G5") {
      if( Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
          Aname=="O3'"  || Aname=="C3'"  || Aname=="O2'"  || Aname=="C2'"  ||
          Aname=="C1'"  || Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  ||
          Aname=="H3'"  || Aname=="H2'"  || Aname=="H1'"  || Aname=="HO5'" ||
          Aname=="HO2'" || Aname=="H5'1" || Aname=="H5'2" || Aname=="HO'2" ||
          Aname=="H2'1" || Aname=="H5T" ) {
        bead_type_[bead_of_[i]]=BB_RNA_5;
      } else if(Aname=="N1" || Aname=="C2"  || Aname=="N2" || Aname=="N3" ||
                Aname=="C4" || Aname=="C5"  || Aname=="C6" || Aname=="O6" ||
                Aname=="N7" || Aname=="C8"  || Aname=="N9" || Aname=="H1" ||
                Aname=="H8" || Aname=="H21" || Aname=="H22" ) {
        bead_type_[bead_of_[i]]=BASE_G;
      } else {
        fail(error, "Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // RNA - GT
    } else if(Rname=="GT") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3"  ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" || Aname=="HOP3" ) {
        bead_type_[bead_of_[i]]=BB_PO3;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="O2'"  || Aname=="C2'"  ||
                Aname=="C1'"  || Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  ||
                Aname=="H3'"  || Aname=="H2'"  || Aname=="H1'"  || Aname=="H3T"  ||
                Aname=="HO5'" || Aname=="HO3'" || Aname=="HO2'" || Aname=="H5'1" ||
                Aname=="H5'2" || Aname=="HO'2" || Aname=="H2'1" || Aname=="H5T" ) {
        bead_type_[bead_of_[i]]=BB_RNA;
      } else if( Aname=="N1" || Aname=="C2"  || Aname=="N2" || Aname=="N3" ||
                 Aname=="C4" || Aname=="C5"  || Aname=="C6" || Aname=="O6" ||
                 Aname=="N7" || Aname=="C8"  || Aname=="N9" || Aname=="H1" ||
                 Aname=="H8" || Aname=="H21" || Aname=="H22" ) {
        bead_type_[bead_of_[i]]=BASE_G;
      } else {
        fail(error, "Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // RNA - U
    } else if(Rname=="U") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="HP" ||
          Aname=="O1P" || Aname=="O2P"  ) {
        bead_type_[bead_of_[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="O2'"  || Aname=="C2'"  ||
                Aname=="C1'"  || Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  ||
                Aname=="H3'"  || Aname=="H2'"  || Aname=="H1'"  || Aname=="H3T"  ||
                Aname=="HO5'" || Aname=="HO3'" || Aname=="HO2'" || Aname=="H5'1" ||
                Aname=="H5'2" || Aname=="HO'2" || Aname=="H2'1" || Aname=="H5T" ) {
        bead_type_[bead_of_[i]]=BB_RNA;
      } else if( Aname=="N1" || Aname=="C2"  || Aname=="O2" || Aname=="N3" ||
                 Aname=="C4" || Aname=="O4"  || Aname=="C5" || Aname=="C6" ||
                 Aname=="H3" || Aname=="H5"  || Aname=="H6") {
        bead_type_[bead_of_[i]]=BASE_U;
      } else {
        fail(error, "Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // RNA - U3
    } else if(Rname=="U3") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3" ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" ) {
        bead_type_[bead_of_[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'" || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'" || Aname=="O2'"  || Aname=="C2'"  ||
                Aname=="C1'"  || Aname=="H5'" || Aname=="H5''" || Aname=="H4'"  ||
                Aname=="H3'"  || Aname=="H2'" || Aname=="H1'"  || Aname=="H3T"  ||
                Aname=="H2'1" || Aname=="HO3'"|| Aname=="HO2'" || Aname=="H5'1" ||
                Aname=="H5'2" || Aname=="HO'2" ) {
        bead_type_[bead_of_[i]]=BB_RNA_3;
      } else if(Aname=="N1" || Aname=="C2"  || Aname=="O2" || Aname=="N3" ||
                Aname=="C4" || Aname=="O4"  || Aname=="C5" || Aname=="C6" ||
                Aname=="H3" || Aname=="H5"  || Aname=="H6") {
        bead_type_[bead_of_[i]]=BASE_U;
      } else {
        fail(error, "Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // RNA - U5
    } else if(Rname=="U5") {
      if( Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
          Aname=="O3'"  || Aname=="C3'"  || Aname=="O2'"  || Aname=="C2'"  ||
          Aname=="C1'"  || Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  ||
          Aname=="H3'"  || Aname=="H2'"  || Aname=="H1'"  || Aname=="HO5'" ||
          Aname=="HO2'" || Aname=="H5'1" || Aname=="H5'2" || Aname=="HO'2" ||
          Aname=="H2'1" || Aname=="H5T" ) {
        bead_type_[bead_of_[i]]=BB_RNA_5;
      } else if(Aname=="N1" || Aname=="C2"  || Aname=="O2" || Aname=="N3" ||
                Aname=="C4" || Aname=="O4"  || Aname=="C5" || Aname=="C6" ||
                Aname=="H3" || Aname=="H5"  || Aname=="H6") {
        bead_type_[bead_of_[i]]=BASE_U;
      } else {
        fail(error, "Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // RNA - UT
    } else if(Rname=="UT") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3"  ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" || Aname=="HOP3" ) {
        bead_type_[bead_of_[i]]=BB_PO3;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="O2'"  || Aname=="C2'"  ||
                Aname=="C1'"  || Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  ||
                Aname=="H3'"  || Aname=="H2'"  || Aname=="H1'"  || Aname=="H3T"  ||
                Aname=="HO5'" || Aname=="HO3'" || Aname=="HO2'" || Aname=="H5'1" ||
                Aname=="H5'2" || Aname=="HO'2" || Aname=="H2'1" || Aname=="H5T" ) {
        bead_type_[bead_of_[i]]=BB_RNA;
      } else if( Aname=="N1" || Aname=="C2"  || Aname=="O2" || Aname=="N3" ||
                 Aname=="C4" || Aname=="O4"  || Aname=="C5" || Aname=="C6" ||
                 Aname=="H3" || Aname=="H5"  || Aname=="H6") {
        bead_type_[bead_of_[i]]=BASE_U;
      } else {
        fail(error, "Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // RNA - A
    } else if(Rname=="A") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="HP" ||
          Aname=="O1P" || Aname=="O2P"  ) {
        bead_type_[bead_of_[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="O2'"  || Aname=="C2'"  ||
                Aname=="C1'"  || Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  ||
                Aname=="H3'"  || Aname=="H2'"  || Aname=="H1'"  || Aname=="H3T"  ||
                Aname=="HO5'" || Aname=="HO3'" || Aname=="HO2'" || Aname=="H5'1" ||
                Aname=="H5'2" || Aname=="HO'2" || Aname=="H2'1" || Aname=="H5T" ) {
        bead_type_[bead_of_[i]]=BB_RNA;
      } else if(Aname=="N1"  || Aname=="C2" || Aname=="N3" || Aname=="C4" ||
                Aname=="C5"  || Aname=="C6" || Aname=="N6" || Aname=="N7" ||
                Aname=="C8"  || Aname=="N9" || Aname=="H2" || Aname=="H8" ||
                Aname=="H61" || Aname=="H62" ) {
        bead_type_[bead_of_[i]]=BASE_A;
      } else {
        fail(error, "Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // RNA - A3
    } else if(Rname=="A3") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3" ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" ) {
        bead_type_[bead_of_[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'" || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'" || Aname=="O2'"  || Aname=="C2'"  ||
                Aname=="C1'"  || Aname=="H5'" || Aname=="H5''" || Aname=="H4'"  ||
                Aname=="H3'"  || Aname=="H2'" || Aname=="H1'"  || Aname=="H3T"  ||
                Aname=="H2'1" || Aname=="HO3'"|| Aname=="HO2'" || Aname=="H5'1" ||
                Aname=="H5'2" || Aname=="HO'2" ) {
        bead_type_[bead_of_[i]]=BB_RNA_3;
      } else if(Aname=="N1"  || Aname=="C2" || Aname=="N3" || Aname=="C4" ||
                Aname=="C5"  || Aname=="C6" || Aname=="N6" || Aname=="N7" ||
                Aname=="C8"  || Aname=="N9" || Aname=="H2" || Aname=="H8" ||
                Aname=="H61" || Aname=="H62" ) {
        bead_type_[bead_of_[i]]=BASE_A;
      } else {
        fail(error, "Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // RNA - A5
    } else if(Rname=="A5") {
      if( Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
          Aname=="O3'"  || Aname=="C3'"  || Aname=="O2'"  || Aname=="C2'"  ||
          Aname=="C1'"  || Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  ||
          Aname=="H3'"  || Aname=="H2'"  || Aname=="H1'"  || Aname=="HO5'" ||
          Aname=="HO2'" || Aname=="H5'1" || Aname=="H5'2" || Aname=="HO'2" ||
          Aname=="H2'1" || Aname=="H5T" ) {
        bead_type_[bead_of_[i]]=BB_RNA_5;
      } else if(Aname=="N1"  || Aname=="C2" || Aname=="N3" || Aname=="C4" ||
                Aname=="C5"  || Aname=="C6" || Aname=="N6" || Aname=="N7" ||
                Aname=="C8"  || Aname=="N9" || Aname=="H2" || Aname=="H8" ||
                Aname=="H61" || Aname=="H62" ) {
        bead_type_[bead_of_[i]]=BASE_A;
      } else {
        fail(error, "Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // RNA - AT
    } else if(Rname=="AT") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3"  ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" || Aname=="HOP3" ) {
        bead_type_[bead_of_[i]]=BB_PO3;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="O2'"  || Aname=="C2'"  ||
                Aname=="C1'"  || Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  ||
                Aname=="H3'"  || Aname=="H2'"  || Aname=="H1'"  || Aname=="H3T"  ||
                Aname=="HO5'" || Aname=="HO3'" || Aname=="HO2'" || Aname=="H5'1" ||
                Aname=="H5'2" || Aname=="HO'2" || Aname=="H2'1" || Aname=="H5T" ) {
        bead_type_[bead_of_[i]]=BB_RNA;
      } else if(Aname=="N1"  || Aname=="C2" || Aname=="N3" || Aname=="C4" ||
                Aname=="C5"  || Aname=="C6" || Aname=="N6" || Aname=="N7" ||
                Aname=="C8"  || Aname=="N9" || Aname=="H2" || Aname=="H8" ||
                Aname=="H61" || Aname=="H62" ) {
        bead_type_[bead_of_[i]]=BASE_A;
      } else {
        fail(error, "Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // RNA - C
    } else if(Rname=="C") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="HP" ||
          Aname=="O1P" || Aname=="O2P"  ) {
        bead_type_[bead_of_[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="O2'"  || Aname=="C2'"  ||
                Aname=="C1'"  || Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  ||
                Aname=="H3'"  || Aname=="H2'"  || Aname=="H1'"  || Aname=="H3T"  ||
                Aname=="HO5'" || Aname=="HO3'" || Aname=="HO2'" || Aname=="H5'1" ||
                Aname=="H5'2" || Aname=="HO'2" || Aname=="H2'1" || Aname=="H5T" ) {
        bead_type_[bead_of_[i]]=BB_RNA;
      } else if(Aname=="N1" || Aname=="C2" || Aname=="O2"  || Aname=="N3" ||
                Aname=="C4" || Aname=="N4" || Aname=="C5"  || Aname=="C6" ||
                Aname=="H5" || Aname=="H6" || Aname=="H41" || Aname=="H42" ) {
        bead_type_[bead_of_[i]]=BASE_C;
      } else {
        fail(error, "Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // RNA - C3
    } else if(Rname=="C3") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3" ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" ) {
        bead_type_[bead_of_[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'" || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'" || Aname=="O2'"  || Aname=="C2'"  ||
                Aname=="C1'"  || Aname=="H5'" || Aname=="H5''" || Aname=="H4'"  ||
                Aname=="H3'"  || Aname=="H2'" || Aname=="H1'"  || Aname=="H3T"  ||
                Aname=="H2'1" || Aname=="HO3'"|| Aname=="HO2'" || Aname=="H5'1" ||
                Aname=="H5'2" || Aname=="HO'2" ) {
        bead_type_[bead_of_[i]]=BB_RNA_3;
      } else if(Aname=="N1" || Aname=="C2" || Aname=="O2"  || Aname=="N3" ||
                Aname=="C4" || Aname=="N4" || Aname=="C5"  || Aname=="C6" ||
                Aname=="H5" || Aname=="H6" || Aname=="H41" || Aname=="H42" ) {
        bead_type_[bead_of_[i]]=BASE_C;
      } else {
        fail(error, "Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // RNA - C5
    } else if(Rname=="C5") {
      if( Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
          Aname=="O3'"  || Aname=="C3'"  || Aname=="O2'"  || Aname=="C2'"  ||
          Aname=="C1'"  || Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  ||
          Aname=="H3'"  || Aname=="H2'"  || Aname=="H1'"  || Aname=="HO5'" ||
          Aname=="HO2'" || Aname=="H5'1" || Aname=="H5'2" || Aname=="HO'2" ||
          Aname=="H2'1" || Aname=="H5T" ) {
        bead_type_[bead_of_[i]]=BB_RNA_5;
      } else if(Aname=="N1" || Aname=="C2" || Aname=="O2"  || Aname=="N3" ||
                Aname=="C4" || Aname=="N4" || Aname=="C5"  || Aname=="C6" ||
                Aname=="H5" || Aname=="H6" || Aname=="H41" || Aname=="H42" ) {
        bead_type_[bead_of_[i]]=BASE_C;
      } else {
        fail(error, "Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // RNA - CT
    } else if(Rname=="CT") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3"  ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" || Aname=="HOP3" ) {
        bead_type_[bead_of_[i]]=BB_PO3;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="O2'"  || Aname=="C2'"  ||
                Aname=="C1'"  || Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  ||
                Aname=="H3'"  || Aname=="H2'"  || Aname=="H1'"  || Aname=="H3T"  ||
                Aname=="HO5'" || Aname=="HO3'" || Aname=="HO2'" || Aname=="H5'1" ||
                Aname=="H5'2" || Aname=="HO'2" || Aname=="H2'1" || Aname=="H5T" ) {
        bead_type_[bead_of_[i]]=BB_RNA;
      } else if(Aname=="N1" || Aname=="C2" || Aname=="O2"  || Aname=="N3" ||
                Aname=="C4" || Aname=="N4" || Aname=="C5"  || Aname=="C6" ||
                Aname=="H5" || Aname=="H6" || Aname=="H41" || Aname=="H42" ) {
        bead_type_[bead_of_[i]]=BASE_C;
      } else {
        fail(error, "Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - G
    } else if(Rname=="DG") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="HP" ||
          Aname=="O1P" || Aname=="O2P"  ) {
        bead_type_[bead_of_[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  || Aname=="C1'"  ||
                Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
                Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO5'" ||
                Aname=="HO3'" || Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" ||
                Aname=="H2'2" || Aname=="H5T"  || Aname=="H3T" ) {
        bead_type_[bead_of_[i]]=BB_DNA;
      } else if(Aname=="N1" || Aname=="C2"  || Aname=="N2" || Aname=="N3" ||
                Aname=="C4" || Aname=="C5"  || Aname=="C6" || Aname=="O6" ||
                Aname=="N7" || Aname=="C8"  || Aname=="N9" || Aname=="H1" ||
                Aname=="H8" || Aname=="H21" || Aname=="H22" ) {
        bead_type_[bead_of_[i]]=BASE_G;
      } else {
        fail(error, "Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - G3
    } else if(Rname=="DG3") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3" ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" ) {
        bead_type_[bead_of_[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  || Aname=="C1'"  ||
                Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
                Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO3'" ||
                Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" || Aname=="H2'2" ||
                Aname=="H3T" ) {
        bead_type_[bead_of_[i]]=BB_DNA_3;
      } else if(Aname=="N1" || Aname=="C2"  || Aname=="N2" || Aname=="N3" ||
                Aname=="C4" || Aname=="C5"  || Aname=="C6" || Aname=="O6" ||
                Aname=="N7" || Aname=="C8"  || Aname=="N9" || Aname=="H1" ||
                Aname=="H8" || Aname=="H21" || Aname=="H22" ) {
        bead_type_[bead_of_[i]]=BASE_G;
      } else {
        fail(error, "Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - G5
    } else if(Rname=="DG5") {
      if( Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
          Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  ||  Aname=="C1'" ||
          Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
          Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO5'" ||
          Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" || Aname=="H2'2" ||
          Aname=="H5T" ) {
        bead_type_[bead_of_[i]]=BB_DNA_5;
      } else if(Aname=="N1" || Aname=="C2"  || Aname=="N2" || Aname=="N3" ||
                Aname=="C4" || Aname=="C5"  || Aname=="C6" || Aname=="O6" ||
                Aname=="N7" || Aname=="C8"  || Aname=="N9" || Aname=="H1" ||
                Aname=="H8" || Aname=="H21" || Aname=="H22" ) {
        bead_type_[bead_of_[i]]=BASE_G;
      } else {
        fail(error, "Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - GT
    } else if(Rname=="DGT") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3"  ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" || Aname=="HOP3" ) {
        bead_type_[bead_of_[i]]=BB_PO3;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  || Aname=="C1'"  ||
                Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
                Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO5'" ||
                Aname=="HO3'" || Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" ||
                Aname=="H2'2" || Aname=="H5T"  || Aname=="H3T" ) {
        bead_type_[bead_of_[i]]=BB_DNA;
      } else if(Aname=="N1" || Aname=="C2"  || Aname=="N2" || Aname=="N3" ||
                Aname=="C4" || Aname=="C5"  || Aname=="C6" || Aname=="O6" ||
                Aname=="N7" || Aname=="C8"  || Aname=="N9" || Aname=="H1" ||
                Aname=="H8" || Aname=="H21" || Aname=="H22" ) {
        bead_type_[bead_of_[i]]=BASE_G;
      } else {
        fail(error, "Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - T
    } else if(Rname=="DT") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="HP" ||
          Aname=="O1P" || Aname=="O2P"  ) {
        bead_type_[bead_of_[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  || Aname=="C1'"  ||
                Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
                Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO5'" ||
                Aname=="HO3'" || Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" ||
                Aname=="H2'2" || Aname=="H5T"  || Aname=="H3T" ) {
        bead_type_[bead_of_[i]]=BB_DNA;
      } else if(Aname=="N1"  || Aname=="C2"  || Aname=="O2"  || Aname=="N3"  ||
                Aname=="C4"  || Aname=="O4"  || Aname=="C5"  || Aname=="C6"  ||
                Aname=="C7"  || Aname=="H3"  || Aname=="H6"  || Aname=="H71" ||
                Aname=="H72" || Aname=="H73" ) {
        bead_type_[bead_of_[i]]=BASE_T;
      } else {
        fail(error, "Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - T3
    } else if(Rname=="DT3") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3" ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" ) {
        bead_type_[bead_of_[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  || Aname=="C1'"  ||
                Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
                Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO3'" ||
                Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" || Aname=="H2'2" ||
                Aname=="H3T" ) {
        bead_type_[bead_of_[i]]=BB_DNA_3;
      } else if(Aname=="N1"  || Aname=="C2"  || Aname=="O2"  || Aname=="N3"  ||
                Aname=="C4"  || Aname=="O4"  || Aname=="C5"  || Aname=="C6"  ||
                Aname=="C7"  || Aname=="H3"  || Aname=="H6"  || Aname=="H71" ||
                Aname=="H72" || Aname=="H73" ) {
        bead_type_[bead_of_[i]]=BASE_T;
      } else {
        fail(error, "Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - T5
    } else if(Rname=="DT5") {
      if( Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
          Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  ||  Aname=="C1'" ||
          Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
          Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO5'" ||
          Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" || Aname=="H2'2" ||
          Aname=="H5T" ) {
        bead_type_[bead_of_[i]]=BB_DNA_5;
      } else if(Aname=="N1"  || Aname=="C2"  || Aname=="O2"  || Aname=="N3"  ||
                Aname=="C4"  || Aname=="O4"  || Aname=="C5"  || Aname=="C6"  ||
                Aname=="C7"  || Aname=="H3"  || Aname=="H6"  || Aname=="H71" ||
                Aname=="H72" || Aname=="H73" ) {
        bead_type_[bead_of_[i]]=BASE_T;
      } else {
        fail(error, "Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - TT
    } else if(Rname=="DTT") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3"  ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" || Aname=="HOP3" ) {
        bead_type_[bead_of_[i]]=BB_PO3;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  || Aname=="C1'"  ||
                Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
                Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO5'" ||
                Aname=="HO3'" || Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" ||
                Aname=="H2'2" || Aname=="H5T"  || Aname=="H3T" ) {
        bead_type_[bead_of_[i]]=BB_DNA;
      } else if(Aname=="N1"  || Aname=="C2"  || Aname=="O2"  || Aname=="N3"  ||
                Aname=="C4"  || Aname=="O4"  || Aname=="C5"  || Aname=="C6"  ||
                Aname=="C7"  || Aname=="H3"  || Aname=="H6"  || Aname=="H71" ||
                Aname=="H72" || Aname=="H73" ) {
        bead_type_[bead_of_[i]]=BASE_T;
      } else {
        fail(error, "Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - A
    } else if(Rname=="DA") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="HP" ||
          Aname=="O1P" || Aname=="O2P"  ) {
        bead_type_[bead_of_[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  || Aname=="C1'"  ||
                Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
                Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO5'" ||
                Aname=="HO3'" || Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" ||
                Aname=="H2'2" || Aname=="H5T"  || Aname=="H3T" ) {
        bead_type_[bead_of_[i]]=BB_DNA;
      } else if(Aname=="N1"  || Aname=="C2" || Aname=="N3" || Aname=="C4" ||
                Aname=="C5"  || Aname=="C6" || Aname=="N6" || Aname=="N7" ||
                Aname=="C8"  || Aname=="N9" || Aname=="H2" || Aname=="H8" ||
                Aname=="H61" || Aname=="H62" ) {
        bead_type_[bead_of_[i]]=BASE_A;
      } else {
        fail(error, "Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - A3
    } else if(Rname=="DA3") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3" ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" ) {
        bead_type_[bead_of_[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  || Aname=="C1'"  ||
                Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
                Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO3'" ||
                Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" || Aname=="H2'2" ||
                Aname=="H3T" ) {
        bead_type_[bead_of_[i]]=BB_DNA_3;
      } else if(Aname=="N1"  || Aname=="C2" || Aname=="N3" || Aname=="C4" ||
                Aname=="C5"  || Aname=="C6" || Aname=="N6" || Aname=="N7" ||
                Aname=="C8"  || Aname=="N9" || Aname=="H2" || Aname=="H8" ||
                Aname=="H61" || Aname=="H62" ) {
        bead_type_[bead_of_[i]]=BASE_A;
      } else {
        fail(error, "Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - A5
    } else if(Rname=="DA5") {
      if( Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
          Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  ||  Aname=="C1'" ||
          Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
          Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO5'" ||
          Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" || Aname=="H2'2" ||
          Aname=="H5T" ) {
        bead_type_[bead_of_[i]]=BB_DNA_5;
      } else if(Aname=="N1"  || Aname=="C2" || Aname=="N3" || Aname=="C4" ||
                Aname=="C5"  || Aname=="C6" || Aname=="N6" || Aname=="N7" ||
                Aname=="C8"  || Aname=="N9" || Aname=="H2" || Aname=="H8" ||
                Aname=="H61" || Aname=="H62" ) {
        bead_type_[bead_of_[i]]=BASE_A;
      } else {
        fail(error, "Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - AT
    } else if(Rname=="DAT") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3"  ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" || Aname=="HOP3" ) {
        bead_type_[bead_of_[i]]=BB_PO3;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  || Aname=="C1'"  ||
                Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
                Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO5'" ||
                Aname=="HO3'" || Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" ||
                Aname=="H2'2" || Aname=="H5T"  || Aname=="H3T" ) {
        bead_type_[bead_of_[i]]=BB_DNA;
      } else if(Aname=="N1"  || Aname=="C2" || Aname=="N3" || Aname=="C4" ||
                Aname=="C5"  || Aname=="C6" || Aname=="N6" || Aname=="N7" ||
                Aname=="C8"  || Aname=="N9" || Aname=="H2" || Aname=="H8" ||
                Aname=="H61" || Aname=="H62" ) {
        bead_type_[bead_of_[i]]=BASE_A;
      } else {
        fail(error, "Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - C
    } else if(Rname=="DC") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="HP" ||
          Aname=="O1P" || Aname=="O2P"  ) {
        bead_type_[bead_of_[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  || Aname=="C1'"  ||
                Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
                Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO5'" ||
                Aname=="HO3'" || Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" ||
                Aname=="H2'2" || Aname=="H5T"  || Aname=="H3T" ) {
        bead_type_[bead_of_[i]]=BB_DNA;
      } else if(Aname=="N1" || Aname=="C2" || Aname=="O2"  || Aname=="N3" ||
                Aname=="C4" || Aname=="N4" || Aname=="C5"  || Aname=="C6" ||
                Aname=="H5" || Aname=="H6" || Aname=="H41" || Aname=="H42" ) {
        bead_type_[bead_of_[i]]=BASE_C;
      } else {
        fail(error, "Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - C3
    } else if(Rname=="DC3") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3" ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" ) {
        bead_type_[bead_of_[i]]=BB_PO2;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  || Aname=="C1'"  ||
                Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
                Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO3'" ||
                Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" || Aname=="H2'2" ||
                Aname=="H3T" ) {
        bead_type_[bead_of_[i]]=BB_DNA_3;
      } else if(Aname=="N1" || Aname=="C2" || Aname=="O2"  || Aname=="N3" ||
                Aname=="C4" || Aname=="N4" || Aname=="C5"  || Aname=="C6" ||
                Aname=="H5" || Aname=="H6" || Aname=="H41" || Aname=="H42" ) {
        bead_type_[bead_of_[i]]=BASE_C;
      } else {
        fail(error, "Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - C5
    } else if(Rname=="DC5") {
      if( Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
          Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  ||  Aname=="C1'" ||
          Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
          Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO5'" ||
          Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" || Aname=="H2'2" ||
          Aname=="H5T" ) {
        bead_type_[bead_of_[i]]=BB_DNA_5;
      } else if(Aname=="N1" || Aname=="C2" || Aname=="O2"  || Aname=="N3" ||
                Aname=="C4" || Aname=="N4" || Aname=="C5"  || Aname=="C6" ||
                Aname=="H5" || Aname=="H6" || Aname=="H41" || Aname=="H42" ) {
        bead_type_[bead_of_[i]]=BASE_C;
      } else {
        fail(error, "Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
      // DNA - CT
    } else if(Rname=="DCT") {
      if( Aname=="P"   || Aname=="OP1" || Aname=="OP2" || Aname=="OP3"  ||
          Aname=="O1P" || Aname=="O2P" || Aname=="O3P" || Aname=="HOP3" ) {
        bead_type_[bead_of_[i]]=BB_PO3;
      } else if(Aname=="O5'"  || Aname=="C5'"  || Aname=="O4'"  || Aname=="C4'"  ||
                Aname=="O3'"  || Aname=="C3'"  || Aname=="C2'"  || Aname=="C1'"  ||
                Aname=="H5'"  || Aname=="H5''" || Aname=="H4'"  || Aname=="H3'"  ||
                Aname=="H2'"  || Aname=="H2''" || Aname=="H1'"  || Aname=="HO5'" ||
                Aname=="HO3'" || Aname=="H5'1" || Aname=="H5'2" || Aname=="H2'1" ||
                Aname=="H2'2" || Aname=="H5T"  || Aname=="H3T" ) {
        bead_type_[bead_of_[i]]=BB_DNA;
      } else if(Aname=="N1" || Aname=="C2" || Aname=="O2"  || Aname=="N3" ||
                Aname=="C4" || Aname=="N4" || Aname=="C5"  || Aname=="C6" ||
                Aname=="H5" || Aname=="H6" || Aname=="H41" || Aname=="H42" ) {
        bead_type_[bead_of_[i]]=BASE_C;
      } else {
        fail(error, "Atom name "+Aname+" is not defined for residue "+Rname+". Check the PDB.");
      }
    } else if(glyc_type_[i] >= 0) {
      // monosaccharide: the bead type was resolved once, from the residue name and
      // from the presence of an N-acetyl group, in the constructor
      bead_type_[bead_of_[i]] = static_cast<unsigned>(glyc_type_[i]);
    } else if(isUnsupportedGlycan(Rname)) {
      fail(error, "Residue "+Rname+" is a monosaccharide for which no ONEBEAD parameters are available.");
    } else if(ionTypeFromResname(Rname) >= 0) {
      // monoatomic ion: assign its onebead bead type. The SAXS/SANS terms for this
      // slot are filled analytically (see SAXS::getOnebeadIonparam_saxs/_sans).
      bead_type_[bead_of_[i]] = obSlot(ionTypeFromResname(Rname));
    } else {
      fail(error, "Residue not known: "+Rname);
    }
  }
}

void OnebeadMapping::getBeadPositions(const std::vector<Vector> &pos, std::vector<Vector> &beads) const {
  const unsigned nbeads = getNumberOfBeads();
  beads.resize(nbeads);
  #pragma omp parallel for num_threads(OpenMP::getNumThreads()) schedule(static)
  for(unsigned r=0; r<nbeads; ++r) {
    Vector sum_pos(0,0,0);
    for(unsigned e=bead_off_[r]; e<bead_off_[r+1]; ++e) {
      const unsigned a = bead_atoms_[e];
      sum_pos += mass_[a] * pos[a];
    }
    beads[r] = sum_pos*bead_inv_mass_[r];
  }
}

}
}
