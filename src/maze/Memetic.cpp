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
 * @file Memetic.cpp
 * @author J. Rydzewski (jr@fizyka.umk.pl)
 */

#include "Memetic.h"

namespace PLMD {
namespace maze {

//+PLUMEDOC MAZE_OPTIMIZER MAZE_MEMETIC_SAMPLING
/*

Calculates the biasing direction along which the ligand unbinds by
minimizing the \ref MAZE_LOSS function. The optimal biasing direction is
determined by performing a population-based memetic search with local
heuristics.

\par Examples

Every optimizer implemented in the maze module needs a loss function as
an argument, and it should be passed using the \ref MAZE_LOSS keyword.

\plumedfile
MAZE_MEMETIC_SAMPLING ...
  LABEL=ma

  LOSS=l

  N_ITER=10
  OPTIMIZER_STRIDE=200

  CAPACITY=20

  LOCAL_SEARCH_ON
  N_LOCAL_ITER=10
  LOCAL_SEARCH_RATE=0.1
  LOCAL_SEARCH_TYPE=stochastic_hill_climbing

  MUTATION_RATE=0.1
  MATING_RATE=0.7
  CAUCHY_ALPHA=0
  CAUCHY_BETA=0.1

  LIGAND=2635-2646
  PROTEIN=1-2634
... MAZE_MEMETIC_SAMPLING
\endplumedfile

As shown above, each optimizer should be provided with the LIGAND and
the PROTEIN keywords.

The neighbor list can be turned on by providing the NLIST keyword.

*/
//+ENDPLUMEDOC

/**
 * Register MAZE_MEMETIC_SAMPLING.
 */
PLUMED_REGISTER_ACTION(Memetic, "MAZE_MEMETIC_SAMPLING")

void Memetic::registerKeywords(Keywords& keys) {
  Optimizer::registerKeywords(keys);

  keys.add(
    "compulsory",
    "CAPACITY",
    "Sampling set size."
  );

  keys.add(
    "compulsory",
    "MUTATION_RATE",
    "Probability of mutation."
  );

  keys.add(
    "compulsory",
    "MATING_RATE",
    "Probability of mating."
  );

  keys.add(
    "compulsory",
    "CAUCHY_ALPHA",
    "Mean of cauchy distribution for sampling."
  );

  keys.add(
    "compulsory",
    "CAUCHY_BETA",
    "Spread of cauchy distribution for sampling."
  );

  keys.addFlag(
    "LOCAL_SEARCH_ON",
    false,
    "Turn local search on."
  );

  keys.add(
    "optional",
    "N_LOCAL_ITER",
    "Number of local search iterations."
  );

  keys.add(
    "optional",
    "LOCAL_SEARCH_RATE",
    "Rate of mutation in local search."
  );

  keys.add(
    "optional",
    "LOCAL_SEARCH_TYPE",
    "Type of local search."
  );
}

Memetic::Memetic(const ActionOptions& ao)
  : PLUMED_OPT_INIT(ao),
    bound_(true),
    score_worst_(0),
    score_best_(0),
    adaptation_(true),
    coding_len_(3),
    local_search_on_(false)
{
  log.printf("maze> Memetic sampling.\n");

  if (keywords.exists("CAPACITY")) {
    parse("CAPACITY", capacity_);

    plumed_massert(
      capacity_ > 0,
      "maze> CAPACITY should be explicitly specified and positive.\n"
    );

    log.printf(
      "maze> CAPACITY read: %u.\n",
      capacity_
    );
  }

  if (keywords.exists("MUTATION_RATE")) {
    parse("MUTATION_RATE", mutation_rate_);

    plumed_massert(
      mutation_rate_ > 0 && mutation_rate_ <= 1,
      "maze> MUTATION_RATE should be in [0, 1).\n"
    );

    log.printf(
      "maze> MUTATION_RATE read: %f.\n",
      mutation_rate_
    );
  }

  if (keywords.exists("MATING_RATE")) {
    parse("MATING_RATE", mating_rate_);

    plumed_massert(
      mating_rate_ > 0 && mating_rate_ <= 1,
      "maze> MATING_RATE should be in [0, 1).\n"
    );

    log.printf(
      "maze> MATING_RATE read: %f.\n",
      mating_rate_
    );
  }

  if (keywords.exists("CAUCHY_ALPHA")) {
    parse("CAUCHY_ALPHA", cauchy_mean_alpha_);

    log.printf(
      "maze> CAUCHY_ALPHA read: %f.\n",
      cauchy_mean_alpha_
    );
  }

  if (keywords.exists("CAUCHY_BETA")) {
    parse("CAUCHY_BETA", cauchy_mean_beta_);

    plumed_massert(
      cauchy_mean_beta_ > 0,
      "maze> CAUCHY_BETA should be explicitly specified and positive.\n"
    );

    log.printf(
      "maze> CAUCHY_BETA read: %f.\n",
      cauchy_mean_beta_
    );
  }

  if (keywords.exists("LOCAL_SEARCH_ON")) {
    parseFlag("LOCAL_SEARCH_ON", local_search_on_);

    log.printf("maze> LOCAL_SEARCH_ON enabled: %d.\n", local_search_on_);
  }

  if (local_search_on_) {
    parse("N_LOCAL_ITER", n_local_iterations_);

    plumed_massert(
      n_local_iterations_ > 0,
      "maze> N_LOCAL_ITER should be explicitly specified and positive.\n"
    );

    log.printf(
      "maze> N_LOCAL_ITER read: %u.\n",
      n_local_iterations_
    );

    parse("LOCAL_SEARCH_RATE", local_search_rate_);

    plumed_massert(
      local_search_rate_ > 0 && local_search_rate_ <= 1,
      "maze> LOCAL_SEARCH_RATE should be in [0, 1).\n"
    );

    log.printf(
      "maze> LOCAL_SEARCH_RATE read: %f.\n",
      local_search_rate_
    );

    parse("LOCAL_SEARCH_TYPE", local_search_type_);

    plumed_massert(
      local_search_type_ == "stochastic_hill_climbing" ||
      local_search_type_ == "adaptive_random_search",
      "maze> LOCAL_SEARCH_TYPE should be: "
      "stochastic_hill_climbing, or adaptive_random_search.\n"
    );

    log.printf(
      "maze> LOCAL_SEARCH_TYPE read: %s.\n",
      local_search_type_.c_str()
    );
  }

  set_n_global_iterations(n_iter_);

  set_label("MEMETIC_SAMPLING");

  start_step_0();

  checkRead();
}

Memetic::~Memetic() {
  delete neighbor_list_;

  members_.clear();
}

void Memetic::optimize() {
  Vector t = solve();

  set_opt(t);
  set_opt_value(score());
}

} // namespace maze
} // namespace PLMD
