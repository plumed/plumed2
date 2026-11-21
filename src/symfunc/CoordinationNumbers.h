/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2017 The plumed team
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
#ifndef __PLUMED_symfunc_CoordinationNumbers_h
#define __PLUMED_symfunc_CoordinationNumbers_h

#include "core/ActionShortcut.h"

namespace PLMD {
namespace symfunc {

class CoordinationNumbers : public ActionShortcut {
public:
  static void shortcutKeywords( Keywords& keys );
  static void registerKeywords( Keywords& keys );
/// Reads the atoms (SPECIES, or SPECIESA and SPECIESB), the switching function and MASK of a
/// shortcut, creates the group lab_grp and returns the input of a contact matrix for them
/// (without the MASK, which is returned in mask)
  static std::string contactMatrixInput( const std::string& lab, const std::string& sp_str,
                                         const std::string& spa_str, const std::string& spb_str,
                                         ActionShortcut* action, std::string& mask );
/// The action that computes the sums over the neighbours (coordination numbers, spherical
/// harmonics) without the contact matrix (see NeighborSums.h): base+"_SUMS" on the CPU, or
/// base+"_SUMSACC" with USEGPU if the openaccPTM plugin is loaded. Empty if the contact matrix is
/// needed: with a MASK, or with USEGPU without the plugin
  static std::string fusedAction( bool usegpu, const std::string& mask, const std::string& base, ActionShortcut* action );
/// Creates the contact matrix lab_mat; flags (e.g. " USEGPU SERIAL") is added to its input
  static void expandMatrix( const bool& components, const std::string& lab, const std::string& sp_str,
                            const std::string& spa_str, const std::string& spb_str, ActionShortcut* action,
                            const std::string& flags="" );
  explicit CoordinationNumbers(const ActionOptions&);
};

}
}
#endif
