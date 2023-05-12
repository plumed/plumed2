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
 * @file Optimizer_Bias.cpp
 *
 * @author J. Rydzewski (jr@fizyka.umk.pl)
 *
 * @code
 * @article{rydzewski2018finding,
 *   title={Finding Multiple Reaction Pathways of Ligand Unbinding},
 *   author={Rydzewski, J and Valsson, O},
 *   journal={arXiv preprint arXiv:1808.08089},
 *   year={2018}
 * }
 * @endcode
 */

#include "core/ActionRegister.h"
#include "core/PlumedMain.h"

#include "bias/Bias.h"

#include "Optimizer.h"
#include "Tools.h"

namespace PLMD {
namespace maze {

//+PLUMEDOC MAZE_BIAS MAZE_OPTIMIZER_BIAS
/*

Biases the ligand along the direction calculated by the chosen \ref MAZE_OPTIMIZER.

OptimizerBias is a class deriving from Bias, and it is used to adaptively
bias a ligand-protein system toward an optimal solution found by the chosen
optimizer.

Remember to define the loss function (\ref MAZE_LOSS) and the optimizer
(\ref MAZE_OPTIMIZER) prior to the adaptive bias for the optimizer.

The adaptive bias potential is the following:
\f[
  V({\bf x}_t)=\alpha
    \left(wt -
      ({\bf x} - {\bf x}^*_{t-\tau})
      \cdot
      \frac{{\bf x}^*_t - {\bf x}_t}{\|{\bf x}^*_t-{\bf x}_t\|}
    \right)^2,
\f]
where \f${\bf x}^*_t\f$ is the optimal solution at time \f$t\f$, \f$w\f$ is the
biasing rate, \f$\tau\f$ is the interval at which the loss function is minimized,
and \f$\alpha\f$ is a scaled force constant.

\par Examples

In the following example the bias potential biases a ligand atom (which have to be
given as an argument) with the biasing rate equal to 0.02 A/ps, and the biasing
constant equal to 3.6 kcal/(mol A). It also takes an optimizer (see
\ref MAZE_OPTIMIZER).

\plumedfile
UNITS LENGTH=A TIME=ps ENERGY=kcal/mol

p: POSITION ATOM=2635 NOPBC

MAZE_OPTIMIZER_BIAS ...
  LABEL=bias

  ARG=p.x,p.y,p.z

  BIASING_RATE=0.02
  ALPHA=3.6

  OPTIMIZER=opt
... MAZE_OPTIMIZER_BIAS
\endplumedfile

*/
//+ENDPLUMEDOC

/**
 * @class OptimizerBias OptimizerBias.cpp "maze/OptimizerBias.cpp"
 *
 * @brief Adaptive bias for the maze optimizers.
 *
 * OptimizerBias is a class deriving from Bias, and it is used to adaptively
 * bias a system toward an optimal solution found by an optimizer.
 */
class OptimizerBias: public bias::Bias {
public:
  /**
   * Standard PLUMED2 constructor.
   *
   * @param ao ActionOptions&.
   */
  explicit OptimizerBias(const ActionOptions& ao);

  /**
   * Destructor.
   */
  ~OptimizerBias() { /* Nothing to do. */ }

  /**
   * Register PLUMED2 keywords.
   *
   * @param keys Keywords.
   */
  static void registerKeywords(Keywords& keys);

  /**
   * Calculate the adaptive biasing potential for ligand unbinding.
   */
  void calculate();

private:
  /**
   * Biased collective variable with Cartesian components, i.e., position,
   * center of mass.
   */
  std::vector<Value*> args_;

  /**
   * Pointer to the optimizer used to minimize the collective variable for
   * ligand unbinding.
   */
  Optimizer* optimizer_;
  std::vector<Optimizer*> opt_pntrs_;

  //! Adaptive bias potential and the corresponding force.
  double bias_;
  double force_;

  /*
   * Parameters of the adaptive biasing potential:
   *  alpha_            rescaled force constant
   *  biasing_speed     biasing rate
   *  biasing_stride    biasing stride
   *  biasing_direction biasing direction
   */

  //! Rescaled force constant.
  double alpha_;
  //! Biasing rate.
  double biasing_speed_;
  //! Biasing stride.
  int biasing_stride_;

  /**
   * Biasing direction is approximated by an optimal solution found by an
   * optimizer.
   */
  Vector biasing_direction_;

  //! Total distance traveled by biased atoms.
  double total_distance_;

  //! Previous value of the collective variable.
  Vector cv0_;

  /*
   * Pointers to PLUMED2 output components.
   */

  //! Biased collective variable components.
  Value* value_dir_x_;
  Value* value_dir_y_;
  Value* value_dir_z_;

  //! Values of the bias and its force.
  Value* value_bias_;
  Value* value_force_;

  //! Total distance.
  Value* value_total_distance_;
};

// Register OPTIMIZER_BIAS as a keyword for PLUMED2 input files.
PLUMED_REGISTER_ACTION(OptimizerBias, "MAZE_OPTIMIZER_BIAS")

void OptimizerBias::registerKeywords(Keywords& keys) {
  Bias::registerKeywords(keys);

  keys.use("ARG");

  keys.add(
    "compulsory",
    "BIASING_RATE",
    "Biasing rate."
  );

  keys.add(
    "compulsory",
    "ALPHA",
    "Rescaled force constant."
  );

  keys.add(
    "compulsory",
    "OPTIMIZER",
    "Optimization technique to minimize the collective variable for ligand\
     unbinding: RANDOM_WALK,\
                STEERED_MD,\
                RANDOM_ACCELERATION_MD,\
                SIMULATED_ANNEALING,\
                MEMETIC_SAMPLING"
  );

  componentsAreNotOptional(keys);

  keys.addOutputComponent(
    "force2",
    "default",
    "Square of the biasing force."
  );

  keys.addOutputComponent(
    "x",
    "default",
    "Optimal biasing direction: x component."
  );

  keys.addOutputComponent(
    "y",
    "default",
    "Optimal biasing direction: y component."
  );

  keys.addOutputComponent(
    "z",
    "default",
    "Optimal biasing direction: z component."
  );

  keys.addOutputComponent(
    "tdist",
    "default",
    "Total distance traveled by biased atoms."
  );
}

OptimizerBias::OptimizerBias(const ActionOptions& ao)
  : PLUMED_BIAS_INIT(ao),
    bias_(0.0),
    force_(0.0),
    total_distance_(0.0)
{
  log.printf(
    "maze> You are using the maze module of PLUMED2,\
    please read and cite "
  );

  log << plumed.cite("Rydzewski J. and Valsson O., arXiv:1808.08089, 2018");
  log.printf("\n");

  args_ = getArguments();
  log.printf(
    "maze> Number of args %zu\n",
    args_.size()
  );

  if (!args_.empty()) {
    log.printf("maze> With arguments");
    for (unsigned i = 0; i < args_.size(); i++) {
      log.printf(" %s", args_[i]->getName().c_str());
    }
    log.printf("\n");
  }

  if (keywords.exists("ALPHA")) {
    parse("ALPHA", alpha_);

    plumed_massert(
      alpha_>0,
      "maze> ALPHA should be explicitly specified and positive.\n"
    );

    log.printf(
      "maze> ALPHA read: %f [kcal/mol/A].\n",
      alpha_
    );
  }

  if (keywords.exists("BIASING_RATE")) {
    parse("BIASING_RATE", biasing_speed_);

    plumed_massert(
      biasing_speed_>0,
      "maze> BIASING_RATE should be explicitly specified and positive.\n"
    );

    log.printf(
      "maze> BIASING_RATE read: %f [A/ps].\n",
      biasing_speed_
    );
  }

  if (keywords.exists("OPTIMIZER")) {
    std::vector<std::string> opt_labels(0);
    parseVector("OPTIMIZER", opt_labels);

    plumed_massert(
      opt_labels.size() > 0,
      "maze> Problem with OPTIMIZER keyword.\n"
    );

    std::string error_msg = "";
    opt_pntrs_ = tls::get_pointers_labels<Optimizer*>(
                   opt_labels,
                   plumed.getActionSet(),
                   error_msg
                 );

    if (error_msg.size() > 0) {
      plumed_merror(
        "maze> Error in keyword OPTIMIZER of " + getName() + ": " + error_msg
      );
    }

    optimizer_ = opt_pntrs_[0];
    log.printf(
      "maze> Optimizer linked: %s.\n",
      optimizer_->get_label().c_str()
    );

    biasing_stride_=optimizer_->get_optimizer_stride();
  }

  checkRead();

  addComponent("force2");
  componentIsNotPeriodic("force2");

  addComponent("x");
  componentIsNotPeriodic("x");

  addComponent("y");
  componentIsNotPeriodic("y");

  addComponent("z");
  componentIsNotPeriodic("z");

  addComponent("tdist");
  componentIsNotPeriodic("tdist");

  biasing_direction_.zero();
  cv0_.zero();

  value_bias_ = getPntrToComponent("bias");
  value_force_ = getPntrToComponent("force2");

  value_dir_x_ = getPntrToComponent("x");
  value_dir_y_ = getPntrToComponent("y");
  value_dir_z_ = getPntrToComponent("z");

  value_total_distance_=getPntrToComponent("tdist");
}

void OptimizerBias::calculate() {
  // Unpack arguments and optimizers.
  Vector cv(
    args_[0]->get(),
    args_[1]->get(),
    args_[2]->get()
  );

  Vector opt_direction(
    optimizer_->value_x_->get(),
    optimizer_->value_y_->get(),
    optimizer_->value_z_->get()
  );

  if (getStep() == 0) {
    cv0_=cv;
  }

  /*
   * For details see a paper by Rydzewski and Valsson.
   */
  double dot = dotProduct(cv - cv0_, biasing_direction_);
  double delta_cv = biasing_speed_ * getTime() - (dot + total_distance_);

  double sign = tls::sgn(delta_cv);

  bias_ = alpha_ * delta_cv * delta_cv;
  force_ = 2.0 * sign * alpha_ * fabs(delta_cv);

  if (getStep() % biasing_stride_ == 0) {
    biasing_direction_ = opt_direction;
    cv0_ = cv;
    total_distance_ += dot;
  }

  /*
   * Return the biasing force to MD engine.
   */
  setOutputForce(0, force_ * biasing_direction_[0]);
  setOutputForce(1, force_ * biasing_direction_[1]);
  setOutputForce(2, force_ * biasing_direction_[2]);

  /*
   * Set values for PLUMED2 outputs.
   */
  value_bias_->set(bias_);
  value_force_->set(force_);

  value_total_distance_->set(total_distance_);

  value_dir_x_->set(biasing_direction_[0]);
  value_dir_y_->set(biasing_direction_[1]);
  value_dir_z_->set(biasing_direction_[2]);
}

} // namespace maze
} // namespace PLMD
