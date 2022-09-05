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
 * @file Simulated_Annealing.cpp
 *
 * @author J. Rydzewski (jr@fizyka.umk.pl)
 */

#include "core/ActionRegister.h"
#include "Optimizer.h"

namespace PLMD {
namespace maze {

//+PLUMEDOC MAZE_OPTIMIZER MAZE_SIMULATED_ANNEALING
/*

Calculates the biasing direction along which the ligand unbinds by minimizing
the \ref MAZE_LOSS function. The optimal biasing direction is determined by
performing simulated annealing.

\par Examples

Every optimizer implemented in the maze module needs a loss function as an
argument, and it should be passed using the \ref MAZE_LOSS keyword.

In the following example simulated annealing is launched for 1000 iterations
as the optimizer for the loss function every 200 ps. The geometric cooling
scheme is used.

\plumedfile
UNITS LENGTH=A TIME=ps ENERGY=kcal/mol

MAZE_SIMULATED_ANNEALING ...
  LABEL=sa

  LOSS=l

  N_ITER=1000
  OPTIMIZER_STRIDE=200

  PROBABILITY_DECREASER=300
  COOLING=0.95
  COOLING_TYPE=geometric

  LIGAND=2635-2646
  PROTEIN=1-2634
... MAZE_SIMULATED_ANNEALING
\endplumedfile

As shown above, each optimizer should be provided with the LIGAND and
the PROTEIN keywords.

*/
//+ENDPLUMEDOC

/**
 * @class Simulated_Annealing Simulated_Annealing.cpp "maze/Simulated_Annealing.cpp"
 *
 * @brief Perform simulated annealing to compute an optimal bias direction.
 */
class Simulated_Annealing: public Optimizer {
public:
  /**
   * PLMD constructor.
   *
   * @param[in] ao PLMD::ActionOptions&.
   */
  explicit Simulated_Annealing(const ActionOptions& ao);

  /**
   * Destructor.
   */
  ~Simulated_Annealing();

  /**
   * Register PLMD keywords.
   *
   * @param[in] keys Keywords.
   */
  static void registerKeywords(Keywords& keys);

  /**
   * Each class deriving from Optimizer needs to override this function.
   */
  void optimize() override;

  /**
   * Reduce the temperature parameter.
   */
  void decrease_probability(unsigned int);

private:
  //! Temperature parameter.
  double probability_decreaser_;

  //! Cooling factor.
  double cooling_factor_;

  //! Cooling scheme.
  std::string cooling_scheme_;
};

// Register MAZE_SIMULATED_ANNEALING.
PLUMED_REGISTER_ACTION(Simulated_Annealing, "MAZE_SIMULATED_ANNEALING")

void Simulated_Annealing::registerKeywords(Keywords& keys) {
  Optimizer::registerKeywords(keys);

  keys.add(
    "compulsory",
    "PROBABILITY_DECREASER",
    "Temperature-like parameter that is decreased during optimization to modify "
    "the Metropolis-Hastings acceptance probability."
  );

  keys.add(
    "compulsory",
    "COOLING",
    "Reduction factor for PROBABILITY_DECREASER, should be in (0, 1]."
  );

  keys.add(
    "compulsory",
    "COOLING_SCHEME",
    "Cooling scheme: geometric."
  );
}

Simulated_Annealing::Simulated_Annealing(const ActionOptions& ao)
  : PLUMED_OPT_INIT(ao)
{
  log.printf("maze> Simulated annealing optimizer.\n");

  if(keywords.exists("COOLING")) {
    parse("COOLING", cooling_factor_);

    plumed_massert(
      cooling_factor_ > 0 && cooling_factor_ <= 1,
      "maze> COOLING should be in (0, 1]; preferably 0.95.\n"
    );
  }

  if(keywords.exists("PROBABILITY_DECREASER")) {
    parse("PROBABILITY_DECREASER", probability_decreaser_);

    plumed_massert(
      probability_decreaser_ > 0,
      "maze> PROBABILITY_DECREASER should be explicitly specified and positive.\n");
  }

  if(keywords.exists("COOLING_SCHEME")) {
    parse("COOLING_SCHEME", cooling_scheme_);

    log.printf(
      "maze> COOLING_SCHEME read: %s.\n",
      cooling_scheme_.c_str()
    );
  }

  set_label("SIMULATED_ANNEALING");

  // Calculate an optimal direction at the beginning of the MD simulation.
  start_step_0();

  checkRead();
}

Simulated_Annealing::~Simulated_Annealing() {
  delete neighbor_list_;
}

void Simulated_Annealing::decrease_probability(unsigned int time) {
  if (cooling_scheme_ == "linear") {
    probability_decreaser_ -= time * cooling_factor_;
  }
  else if (cooling_scheme_ == "exponential") {
    probability_decreaser_ *= pow(cooling_factor_, time);
  }
  else if (cooling_scheme_ == "geometric") {
    probability_decreaser_ *= cooling_factor_;
  }
  else if (cooling_scheme_ == "logarithmic") {
    probability_decreaser_ = cooling_factor_ / std::log(time + 1);
  }
  else if (cooling_scheme_ == "hoffman") {
    probability_decreaser_ = (cooling_factor_ - 1) / std::log(time);
  }
}

void Simulated_Annealing::optimize() {
  sampling_r_ = sampling_radius();
  double rad_s;
  const unsigned nl_size = neighbor_list_->size();

  Vector distance, distance_next;

  for (unsigned int iter=0; iter < get_n_iterations(); ++iter) {
    double action = 0;
    double action_next = 0;

    rad_s = rnd::next_double(sampling_r_);
    Vector dev = rnd::next_plmd_vector(rad_s);

    #pragma omp parallel num_threads(get_n_threads_openmp())
    {
      #pragma omp for reduction(+:action_next, action)
      for (unsigned int i=0; i < nl_size; i++) {
        unsigned i0 = neighbor_list_->getClosePair(i).first;
        unsigned i1 = neighbor_list_->getClosePair(i).second;

        if (getAbsoluteIndex(i0) == getAbsoluteIndex(i1)) {
          continue;
        }

        if (pbc_) {
          distance = pbcDistance(
                       getPosition(i0) + get_opt(),
                       getPosition(i1)
                     );

          distance_next = pbcDistance(
                            getPosition(i0) + dev,
                            getPosition(i1)
                          );
        }
        else {
          distance = delta(
                       getPosition(i0) + get_opt(),
                       getPosition(i1)
                     );

          distance_next = delta(
                            getPosition(i0) + dev,
                            getPosition(i1)
                          );
        }

        action += pairing(distance.modulo());
        action_next += pairing(distance_next.modulo());
      }
    }

    double p = std::min(
                 1.0,
                 std::exp(-(action_next - action) / probability_decreaser_)
               );

    double r = rnd::next_double();

    if (r < p) {
      set_opt(dev);
      set_opt_value(action_next);
    }

    decrease_probability(iter);
  }

  Vector s = get_opt() / modulo(get_opt());
  set_opt(s);
}

} // namespace maze
} // namespace PLMD
