/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Copyright (c) 2019 Jakub Rydzewski (jr@fizyka.umk.pl). All rights reserved.

See http://www.maze-code.github.io for more information.

This file is part of maze.

maze is free software: you can redistribute it and/or modify it under the
terms of the GNU Lesser General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your option)
any later version.

maze is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE.

See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with maze. If not, see <https://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#ifndef __PLUMED_maze_Loss_h
#define __PLUMED_maze_Loss_h

/**
 * @file Loss.h
 *
 * @author J. Rydzewski (jr@fizyka.umk.pl)
 */

#include "colvar/Colvar.h"
#include "core/ActionRegister.h"
#include "Core.h"

namespace PLMD {
namespace maze {

/**
 * @class Loss Loss.h "maze/Loss.h"
 *
 * @brief Loss function desribes a score between a ligand-protein conformation.
 *
 * Loss function must be defined for an optimizer as it minimizes a loss
 * of a ligand-protein conformation in order to simulate the ligand-protein
 * dissociation process.
 */
class Loss: public colvar::Colvar {
public:
  /**
   * PLMD constructor.
   *
   * @param[in] ao PLMD::ActionOptions&.
   */
  explicit Loss(const ActionOptions& ao);

  /**
   * Destructor.
   */
  ~Loss() { /* Nothing to do. */ }

  /**
   * Register PLMD keywords.
   *
   * @param[in] keys Keywords.
   */
  static void registerKeywords(Keywords& keys);

  /**
   * Calculate a loss of a single pair of ligand-protein atoms.
   *
   * @param[in] distance Distance between atoms in the pair.
   */
  double pairing(double distance);

  // Required by the Colvar class.
  void calculate() override { /* Nothing to do. */ }

protected:
  //! Parameters of the loss function.
  std::vector<double> params_;
};

} // namespace maze
} // namespace PLMD

#endif // __PLUMED_maze_Loss_h
