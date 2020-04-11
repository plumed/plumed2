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
 * @file Random_Acceleration_MD.cpp
 *
 * @author J. Rydzewski (jr@fizyka.umk.pl)
 */

#include "core/ActionRegister.h"
#include "Optimizer.h"

namespace PLMD {
namespace maze {

//+PLUMEDOC MAZE_OPTIMIZER MAZE_RANDOM_ACCELERATION_MD
/*

Performs random acceleration MD within the protein matrix.

\par Examples

Every optimizer implemented in the maze module needs a loss function as
an argument, and it should be passed using the \ref MAZE_LOSS keyword.

\plumedfile
MAZE_RANDOM_ACCELERATION_MD ...
  LABEL=rw

  OPTIMIZER_STRIDE=_
  LOSS=l
  RMIN=_

  LIGAND=2635-2646
  PROTEIN=1-2634
... MAZE_RANDOM_ACCELERATION_MD
\endplumedfile

As shown above, each optimizer should be provided with the LIGAND and
the PROTEIN keywords.

*/
//+ENDPLUMEDOC

/**
 * @class Random_Acceleration_MD Random_Acceleration_MD.cpp
 *  "maze/Random_Acceleration_MD.cpp"
 *
 * @brief Perform RAMD simulation.
 */
class Random_Acceleration_MD: public Optimizer {
public:
  /**
   * PLMD constructor.
   *
   * @param[in] ao PLMD::ActionOptions&
   */
  explicit Random_Acceleration_MD(const ActionOptions&);

  /**
   * Destructor.
   */
  ~Random_Acceleration_MD();

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

private:
  //! Threshold distance that the ligand needs to pass.
  double r_min_;

  //! Total distance.
  double total_dist_;

  //! Distance.
  double dist_;

  //! Ligand center of mass.
  Vector com_;

  //! PLMD value for distance.
  Value* value_dist_;

  //! PLMD value for total distance.
  Value* value_total_dist_;
};

// Register MAZE_RANDOM_ACCELERATION_MD.
PLUMED_REGISTER_ACTION(Random_Acceleration_MD, "MAZE_RANDOM_ACCELERATION_MD")

void Random_Acceleration_MD::registerKeywords(Keywords& keys) {
  Optimizer::registerKeywords(keys);

  keys.remove("N_ITER");

  keys.add(
    "compulsory",
    "R_MIN",
    "Minimal distance traveled before sampling a new direction of biasing."
  );

  keys.addOutputComponent(
    "dist",
    "default",
    "Distance traveled in one sampling interval."
  );

  keys.addOutputComponent(
    "tdist",
    "default",
    "Total distance traveled by biased atoms."
  );
}

Random_Acceleration_MD::Random_Acceleration_MD(const ActionOptions& ao)
  : PLUMED_OPT_INIT(ao),
    total_dist_(0.0),
    dist_(0.0) {
  log.printf("maze> Random accelerated molecular dynamics.\n");

  if(keywords.exists("R_MIN")) {
    parse("R_MIN", r_min_);

    plumed_massert(
      r_min_ > 0,
      "maze> R_MIN should be explicitly specified and positive.\n"
    );

    log.printf(
      "maze> R_MIN read: %f [A].\n",
      r_min_
    );
  }

  set_label("RANDOM_ACCELERATION_MD");
  set_opt(rnd::next_plmd_vector());
  set_opt_value(0.0);

  start_step_stride();

  checkRead();

  com_ = center_of_mass();

  addComponent("dist");
  componentIsNotPeriodic("dist");
  value_dist_ = getPntrToComponent("dist");

  addComponent("tdist");
  componentIsNotPeriodic("tdist");
  value_total_dist_ = getPntrToComponent("tdist");
}

Random_Acceleration_MD::~Random_Acceleration_MD() {
  delete neighbor_list_;
}

void Random_Acceleration_MD::optimize() {
  Vector c = center_of_mass();
  Vector d;

  if (pbc_) {
    d = pbcDistance(c, com_);
  }
  else {
    d = delta(c, com_);
  }

  dist_ = d.modulo();
  total_dist_ += dist_;

  if(dist_ < r_min_) {
    set_opt(rnd::next_plmd_vector());
  }

  set_opt_value(score());
  com_ = c;

  value_dist_->set(dist_);
  value_total_dist_->set(total_dist_);
}

} // namespace maze
} // namespace PLMD
