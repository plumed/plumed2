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
#ifndef __PLUMED_isdb_LCPO_h
#define __PLUMED_isdb_LCPO_h

#include "tools/Vector.h"

#include <cstddef>
#include <functional>
#include <map>
#include <string>
#include <vector>

namespace PLMD {
namespace isdb {

/// Solvent-accessible surface area of the atoms with the LCPO method
/// (Weiser, Shenkin, Still, J. Comput. Chem. 20, 217 (1999)), for proteins,
/// nucleic acids and glycans.
///
/// The parameters of each atom are assigned from its atom and residue names.
/// Atoms of monoatomic ions and hydrogens have no parameters and are skipped.
/// The neighbours are found with a cell list and kept in a Verlet list, which
/// is rebuilt only when an atom has moved by more than about half of its skin;
/// the result is the same as with an all-pairs search.
class LCPO {
public:
  /// Called with a message when the input is not valid. It must not return
  /// (e.g. Action::error); if it does, the constructor throws.
  typedef std::function<void(const std::string &)> ErrorHandler;
  /// The parameters of the atoms (radius in Angstrom, including the probe
  /// radius of 1.4 Angstrom, and the four LCPO coefficients), keyed by
  /// "RESIDUE_ATOM". Glycan oxygens have a second entry with the suffix "_L",
  /// used when they carry a glycosidic bond.
  static std::map<std::string, std::vector<double> > getParameterTable();
  /// Assign the parameters of the atoms. For atom i: atomNames[i] and
  /// residueNames[i] are its names; glycan[i] tells whether it belongs to a
  /// monosaccharide and bridge[i] whether it is an oxygen (or a protein atom)
  /// that carries a glycosidic bond. The N-terminal nitrogen and the C-terminal
  /// oxygen are the first and the last atoms that are neither in a glycan nor
  /// in an ion.
  LCPO(const std::vector<std::string> &atomNames, const std::vector<std::string> &residueNames,
       const std::vector<bool> &glycan, const std::vector<bool> &bridge, const ErrorHandler &error);
  unsigned getNumberOfAtoms() const {
    return param_.size();
  }
  /// True if atom i can have a nonzero area.
  bool hasArea(unsigned i) const {
    return param_[i].size()>1 && param_[i][1]>0.0;
  }
  /// The areas (Angstrom^2) of the atoms with an area, from the positions (nm)
  /// of all the atoms. Only the atoms with an area of index first, first+stride,
  /// ... in the list of the atoms with an area are computed (to split the work
  /// over MPI ranks); area is resized to the number of atoms and is zero for
  /// the others. The neighbour list is updated first.
  void computeAreas(const std::vector<Vector> &pos, unsigned first, unsigned stride, std::vector<double> &area);
private:
  // radius and coefficients of each atom (empty if the atom has none)
  std::vector<std::vector<double> > param_;
  // atoms with parameters, and those with an area
  std::vector<unsigned> atoms_;
  std::vector<unsigned> area_atoms_;
  // Verlet list of the candidate neighbours (CSR), with the positions at which
  // it was built
  std::vector<Vector> ref_pos_;
  std::vector<std::size_t> cand_off_;
  std::vector<unsigned> cand_idx_;
  std::vector<double> cand_d_;
  // the neighbours (CSR) with the overlap term A(i,j) of each directed pair
  std::vector<std::size_t> off_;
  std::vector<unsigned> idx_;
  std::vector<double> term_;
  void buildCandidates(const std::vector<Vector> &pos);
  void updateList(const std::vector<Vector> &pos);
};

}
}

#endif
