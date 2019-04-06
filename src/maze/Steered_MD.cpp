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
 * @file Steered_MD.cpp
 *
 * @author J. Rydzewski (jr@fizyka.umk.pl)
 */

#include "core/ActionRegister.h"
#include "Optimizer.h"

namespace PLMD {
namespace maze {

//+PLUMEDOC MAZE_OPTIMIZER MAZE_STEERED_MD
/*

Performs a linear unbinding along a predefined biasing direction that
needs to be provided using the PULLING keyword.

\par Examples

Every optimizer implemented in the maze module needs a loss function as
an argument, and it should be passed using the \ref MAZE_LOSS keyword.

\plumedfile
MAZE_STEERED_MD ...
  LABEL=smd

  LOSS=l
  PULLING=0.3,0.3,0.3
  OPTIMIZER_STRIDE=_

  LIGAND=2635-2646
  PROTEIN=1-2634
... MAZE_STEERED_MD
\endplumedfile

As shown above, each optimizer should be provided with the LIGAND and
the PROTEIN keywords.

*/
//+ENDPLUMEDOC

/**
 * @class Steered_MD Steered_MD.cpp "maze/Steered_MD.cpp"
 * @brief Performs steered MD simulation.
 */
class Steered_MD: public Optimizer {
public:
  /**
   * PLMD constructor.
   *
   * @param[in] ao PLMD::ActionOptions&
   */
  explicit Steered_MD(const ActionOptions& ao);

  /**
   * Destructor.
   */
  ~Steered_MD();

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
  //! Total distance traveled by the ligand.
  double total_dist_;

  //! Ligand center of mass.
  Vector com_;

  //! Constant direction of biasing.
  Vector pulling_;

  //! PLMD::Value of total distance.
  Value* value_total_dist_;
};

// Register MAZE_STEERED_MD.
PLUMED_REGISTER_ACTION(Steered_MD, "MAZE_STEERED_MD")

void Steered_MD::registerKeywords(Keywords& keys) {
  Optimizer::registerKeywords(keys);

  keys.remove("N_ITER");

  keys.add(
    "compulsory",
    "PULLING",
    "Constant biasing direction."
  );

  keys.addOutputComponent(
    "tdist",
    "default",
    "Total distance traveled by biased atoms."
  );
}

Steered_MD::Steered_MD(const ActionOptions& ao)
  : PLUMED_OPT_INIT(ao),
    total_dist_(0.0)
{
  log.printf("maze> Steered MD.\n");

  if (keywords.exists("PULLING")) {
    std::vector<double> v;
    parseVector("PULLING", v);
    pulling_ = tls::vector2Vector(v);

    log.printf("maze> PULLING read.\n");
  }

  set_label("STEERED_MD");
  set_opt(pulling_);
  set_opt_value(0.0);

  start_step_stride();

  checkRead();

  com_ = center_of_mass();

  addComponent("tdist");
  componentIsNotPeriodic("tdist");
  value_total_dist_ = getPntrToComponent("tdist");
}

Steered_MD::~Steered_MD() {
  delete neighbor_list_;
}

void Steered_MD::optimize() {
  Vector c = center_of_mass();
  Vector d;

  if (pbc_) {
    d = pbcDistance(c, com_);
  }
  else {
    d = delta(c, com_);
  }

  double dist = d.modulo();
  total_dist_ += dist;

  set_opt(pulling_);
  set_opt_value(score());
  com_ = c;

  value_total_dist_->set(total_dist_);
}

} // namespace maze
} // namespace PLMD
