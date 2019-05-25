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

/**
 * @file Random_Walk.cpp
 *
 * @author J. Rydzewski (jr@fizyka.umk.pl)
 */

#include "core/ActionRegister.h"
#include "Optimizer.h"

namespace PLMD {
namespace maze {

//+PLUMEDOC MAZE_OPTIMIZER MAZE_RANDOM_WALK
/*

Fake optimizer that can be used for debugging.

This is dummy optimizer that can be used for debugging and monitoring
purposes. It returns a random direction of biasing, changed every
OPTIMIZER_STRIDE.

Performs a random walk within the protein matrix.

\par Examples

Every optimizer implemented in the maze module needs a loss function as
an argument, and it should be passed using the \ref MAZE_LOSS keyword.

\plumedfile
MAZE_RANDOM_WALK ...
  LABEL=rw

  LOSS=l
  OPTIMIZER_STRIDE=200

  LIGAND=2635-2646
  PROTEIN=1-2634
... MAZE_RANDOM_WALK
\endplumedfile

As shown above, each optimizer should be provided with the LIGAND and
the PROTEIN keywords.

*/
//+ENDPLUMEDOC

/**
 * @class Random_Walk Random_Walk.cpp "maze/Random_Walk.cpp"
 *
 * @brief Perform a random walk within the protein matrix.
 */
class Random_Walk: public Optimizer {
public:
  /**
   * PLMD constructor.
   *
   * @param[in] ao PLMD::ActionOptions&
   */
  explicit Random_Walk(const ActionOptions& ao);

  /**
   * Destructor.
   */
  ~Random_Walk();

  /**
   * Registers PLMD keywords.
   *
   * @param[in] keys PLMD keywords
   */
  static void registerKeywords(Keywords&);

  /**
   * Each class deriving from Optimizer needs to override this function.
   */
  void optimize() override;
};

// Register MAZE_RANDOM_WALK.
PLUMED_REGISTER_ACTION(Random_Walk, "MAZE_RANDOM_WALK")

void Random_Walk::registerKeywords(Keywords& keys) {
  Optimizer::registerKeywords(keys);

  keys.remove("N_ITER");
}

Random_Walk::Random_Walk(const ActionOptions& ao)
  : PLUMED_OPT_INIT(ao)
{
  log.printf("maze> Fake optimizer that returns a next step as random,\
    can be used to monitor loss, and for debugging and regtests purposes.\n");

  set_label("RANDOM_WALK");

  start_step_0();

  checkRead();
}

Random_Walk::~Random_Walk() {
  delete neighbor_list_;
}

void Random_Walk::optimize() {
  set_opt(rnd::next_plmd_vector());
  set_opt_value(score());
}

} // namespace maze
} // namespace PLMD
