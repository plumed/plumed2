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
#ifndef __PLUMED_isdb_OnebeadMapping_h
#define __PLUMED_isdb_OnebeadMapping_h

#include "tools/AtomNumber.h"
#include "tools/Vector.h"

#include <functional>
#include <string>
#include <vector>

namespace PLMD {

class Log;
class PDB;

namespace isdb {

/// Identifiers of the ONEBEAD bead types and of the monoatomic ions, with the
/// residue-name helpers that classify the residues of an all-atom structure.
/// The enumerators are used as indices of parameter tables (e.g. the SAXS/SANS
/// form factors), so their order must not change.
struct OnebeadTypes {
  /// Canonical ion identifiers. The ATOMISTIC (AFF_*) and ONEBEAD (OB_*) ion
  /// slots of SAXS are listed in this same order.
  enum { ION_NA, ION_K, ION_CL, ION_CA, ION_MG, ION_ZN, ION_FE2, ION_FE3, ION_MN, NIONS };
  /// Bead types of the ONEBEAD model.
  enum { TRP,
         TYR,
         PHE,
         HIS,
         HIP,
         ARG,
         LYS,
         CYS,
         CYX,
         ASP,
         GLU,
         ILE,
         LEU,
         MET,
         ASN,
         PRO,
         GLN,
         SER,
         THR,
         VAL,
         ALA,
         GLY,
         BB_PO2,
         BB_PO3,
         BB_DNA,
         BB_DNA_5,
         BB_DNA_3,
         BB_RNA,
         BB_RNA_5,
         BB_RNA_3,
         BASE_A,
         BASE_C,
         BASE_T,
         BASE_G,
         BASE_U,
         // Monosaccharide bead types (one bead per residue, exactly as for the
         // amino acids). GLYC_NAG/GLYC_NGA are the N-acetylated hexosamines
         // (GlcNAc/GalNAc).
         GLYC_FUC, GLYC_MAN, GLYC_BMA, GLYC_GAL, GLYC_GLC, GLYC_NAG, GLYC_NGA, GLYC_SIA,
         // Monoatomic-ion bead types (one atom per bead), in the same order as ION_*.
         OB_NA, OB_K, OB_CL, OB_CA, OB_MG, OB_ZN, OB_FE2, OB_FE3, OB_MN,
         NONEBEAD
       };
  /// The ONEBEAD bead type of an ion
  static unsigned obSlot(unsigned ion) {
    return OB_NA + ion;
  }
  /// The canonical ION_* id for a residue name, or -1 if it is not a known ion.
  static int ionTypeFromResname(const std::string &resname);
  /// The atomic mass (g/mol) of an ion.
  static double ionMass(unsigned ion);
  /// The canonical name of a residue whose name encodes a glycosylation site
  /// (e.g. NLN -> ASN); other names are returned unchanged.
  static std::string onebeadResidueAlias(const std::string &resname);
  /// True when the residue is an amino acid (a single bead), CYX included.
  static bool isOnebeadProteinResidue(const PDB &pdb, const std::string &resname);
  /// The GLYC_* bead type of a monosaccharide residue name, or -1 if the residue
  /// is not a supported glycan. hasN must be true when the residue carries an
  /// amide nitrogen: that is what distinguishes GlcNAc from Glc and GalNAc from
  /// Gal when the input file names an N-acetylhexosamine after its parent hexose.
  static int glycanTypeFromResname(const std::string &resname, bool hasN);
  /// True for monosaccharides that are recognised but not parameterised.
  static bool isUnsupportedGlycan(const std::string &resname);
};

/// The ONEBEAD coarse-grained mapping of an all-atom structure.
///
/// Every amino acid is one bead, every nucleotide is two or three beads (sugar,
/// base and, when present, phosphate), every monosaccharide is one bead and every
/// monoatomic ion is one bead. A bead is the centre of mass of its atoms.
///
/// The mapping is built from an all-atom PDB (the template), read chain by chain
/// and, within a chain, residue by residue. The template must contain exactly the
/// atoms of the list passed to the constructor, in the same order. Atom indices in
/// this class are positions in that list.
class OnebeadMapping :
  public OnebeadTypes {
public:
  /// Called with a message when the input is not valid. It must not return
  /// (e.g. Action::error); if it does, the constructor throws.
  typedef std::function<void(const std::string &)> ErrorHandler;
  /// Build the mapping of atoms. The template must be read with scale 1
  /// (Angstrom). A summary of the residue-name aliases and of the glycans is
  /// written to log.
  OnebeadMapping(const PDB &tmpl, const std::vector<AtomNumber> &atoms, Log &log, const ErrorHandler &error);
  unsigned getNumberOfAtoms() const {
    return bead_of_.size();
  }
  unsigned getNumberOfBeads() const {
    return bead_type_.size();
  }
  /// The bead of atom a.
  unsigned getBead(unsigned a) const {
    return bead_of_[a];
  }
  /// The bead type (TRP, ..., OB_MN) of bead b.
  unsigned getBeadType(unsigned b) const {
    return bead_type_[b];
  }
  /// The atoms of bead b are getBeadAtoms()[k] for getBeadOffsets()[b] <= k < getBeadOffsets()[b+1],
  /// in increasing order.
  const std::vector<unsigned> & getBeadOffsets() const {
    return bead_off_;
  }
  const std::vector<unsigned> & getBeadAtoms() const {
    return bead_atoms_;
  }
  /// The mass of atom a (g/mol), from its element.
  double getAtomMass(unsigned a) const {
    return mass_[a];
  }
  /// The mass of atom a divided by the mass of its bead: the derivative of the
  /// bead position with respect to the position of the atom.
  double getAtomWeight(unsigned a) const {
    return atom_weight_[a];
  }
  /// The atom name of atom a, as in the template.
  const std::string & getAtomName(unsigned a) const {
    return atom_name_[a];
  }
  /// The residue name of atom a, after onebeadResidueAlias().
  const std::string & getResidueName(unsigned a) const {
    return residue_name_[a];
  }
  /// The GLYC_* bead type of atom a, or -1 if the atom is not in a monosaccharide.
  int getGlycanType(unsigned a) const {
    return glyc_type_[a];
  }
  /// True for the oxygens bound to two heavy atoms in a monosaccharide (ring and
  /// glycosidic oxygens), and for the protein or nucleic-acid atoms that carry a
  /// glycan.
  bool isGlycanBridge(unsigned a) const {
    return glyc_bridge_[a];
  }
  /// The positions of the beads (centres of mass) from those of the atoms. The
  /// sums are done in the order of the atoms, so the result does not depend on
  /// the number of threads.
  void getBeadPositions(const std::vector<Vector> &pos, std::vector<Vector> &beads) const;
private:
  // per atom
  std::vector<unsigned> bead_of_;
  std::vector<double> mass_;
  std::vector<double> atom_weight_;
  std::vector<std::string> atom_name_;
  std::vector<std::string> residue_name_;
  std::vector<int> glyc_type_;
  std::vector<bool> glyc_bridge_;
  // per bead
  std::vector<unsigned> bead_type_;
  std::vector<unsigned> bead_off_;
  std::vector<unsigned> bead_atoms_;
  std::vector<double> bead_inv_mass_;
  // the bead types, from the residue and atom names
  void assignBeadTypes(const PDB &tmpl, const std::vector<AtomNumber> &atoms, const ErrorHandler &error);
};

}
}

#endif
