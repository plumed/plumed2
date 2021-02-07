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
 * @file Optimizer.cpp
 *
 * @author J. Rydzewski (jr@fizyka.umk.pl)
 */

#include "Optimizer.h"
#include "core/PlumedMain.h"

namespace PLMD {
namespace maze {

void Optimizer::registerKeywords(Keywords& keys) {
  Colvar::registerKeywords(keys);

  keys.addFlag(
    "SERIAL",
    false,
    "Perform the simulation in serial -- used only for debugging purposes, "
    "should not be used otherwise."
  );

  keys.addFlag(
    "PAIR",
    false,
    "Pair only the 1st element of the 1st group with the 1st element in the "
    "second, etc."
  );

  keys.addFlag(
    "NLIST",
    true,
    "Use a neighbor list of ligand-protein atom pairs to speed up the "
    "calculating of the distances."
  );

  keys.add(
    "optional",
    "NL_CUTOFF",
    "Neighbor list cut-off for the distances of ligand-protein atom pairs."
  );

  keys.add(
    "optional",
    "NL_STRIDE",
    "Update stride for the ligand-protein atom pairs in the neighbor list."
  );

  keys.add(
    "compulsory",
    "N_ITER",
    "Number of optimization steps. Required only for optimizers, do not pass "
    "this keyword to the fake optimizers (results in crash) , e.g., random "
    "walk, steered MD, or random acceleration MD."
  );

  keys.add(
    "optional",
    "LOSS",
    "Loss function describing ligand-protein interactions required by every "
    "optimizer."
  );

  keys.add(
    "atoms",
    "LIGAND",
    "Indices of ligand atoms."
  );

  keys.add(
    "atoms",
    "PROTEIN",
    "Indices of protein atoms."
  );

  keys.add(
    "compulsory",
    "OPTIMIZER_STRIDE",
    "Optimizer stride. Sets up a callback function that launches the "
    "optimization process every OPTIMIZER_STRIDE."
  );

  componentsAreNotOptional(keys);

  keys.addOutputComponent(
    "x",
    "default",
    "Optimal biasing direction; x component."
  );

  keys.addOutputComponent(
    "y",
    "default",
    "Optimal biasing direction; y component."
  );

  keys.addOutputComponent(
    "z",
    "default",
    "Optimal biasing direction; z component."
  );

  keys.addOutputComponent(
    "loss",
    "default",
    "Loss function value defined by the provided pairing function."
  );

  keys.addOutputComponent(
    "sr",
    "default",
    "Sampling radius. Reduces sampling to the local proximity of the ligand "
    "position."
  );
}

Optimizer::Optimizer(const ActionOptions& ao)
  : PLUMED_COLVAR_INIT(ao),
    first_step_(true),
    opt_value_(0.0),
    pbc_(true),
    sampling_r_(0.0),
    serial_(false),
    validate_list_(true),
    first_time_(true)
{
  parseFlag("SERIAL", serial_);

  if (keywords.exists("LOSS")) {
    std::vector<std::string> loss_labels(0);
    parseVector("LOSS", loss_labels);

    plumed_massert(
      loss_labels.size() > 0,
      "maze> Something went wrong with the LOSS keyword.\n"
    );

    std::string error_msg = "";
    vec_loss_ = tls::get_pointers_labels<Loss*>(
                  loss_labels,
                  plumed.getActionSet(),
                  error_msg
                );

    if (error_msg.size() > 0) {
      plumed_merror(
        "maze> Error in the LOSS keyword " + getName() + ": " + error_msg
      );
    }

    loss_ = vec_loss_[0];
    log.printf("maze> Loss function linked to the optimizer.\n");
  }

  if (keywords.exists("N_ITER")) {
    parse("N_ITER", n_iter_);

    plumed_massert(
      n_iter_ > 0,
      "maze> N_ITER should be explicitly specified and positive.\n"
    );

    log.printf(
      "maze> Optimizer will run %u iterations once launched.\n",
      n_iter_
    );
  }

  std::vector<AtomNumber> ga_list, gb_list;
  parseAtomList("LIGAND", ga_list);
  parseAtomList("PROTEIN", gb_list);

  bool nopbc = !pbc_;
  parseFlag("NOPBC", nopbc);

  bool do_pair = false;
  parseFlag("PAIR", do_pair);

  nl_stride_ = 0;
  bool do_neigh = false;
  parseFlag("NLIST", do_neigh);

  if (do_neigh) {
    if (keywords.exists("NL_CUTOFF")) {
      parse("NL_CUTOFF", nl_cutoff_);

      plumed_massert(
        nl_cutoff_ > 0,
        "maze> NL_CUTOFF should be explicitly specified and positive.\n"
      );
    }

    if (keywords.exists("NL_STRIDE")) {
      parse("NL_STRIDE", nl_stride_);

      plumed_massert(
        nl_stride_ > 0,
        "maze> NL_STRIDE should be explicitly specified and positive.\n"
      );
    }
  }

  if (gb_list.size() > 0) {
    if (do_neigh) {
      neighbor_list_ = new NeighborList(
        ga_list,
        gb_list,
        serial_,
        do_pair,
        pbc_,
        getPbc(),
        comm,
        nl_cutoff_,
        nl_stride_
      );
    }
    else {
      neighbor_list_=new NeighborList(
        ga_list,
        gb_list,
        serial_,
        do_pair,
        pbc_,
        getPbc(),
        comm
      );
    }
  }
  else {
    if (do_neigh) {
      neighbor_list_ = new NeighborList(
        ga_list,
        serial_,
        pbc_,
        getPbc(),
        comm,
        nl_cutoff_,
        nl_stride_
      );
    }
    else {
      neighbor_list_=new NeighborList(
        ga_list,
        serial_,
        pbc_,
        getPbc(),
        comm
      );
    }
  }

  requestAtoms(neighbor_list_->getFullAtomList());

  log.printf(
    "maze> Loss will be calculated between two groups of %u and %u atoms.\n",
    static_cast<unsigned>(ga_list.size()),
    static_cast<unsigned>(gb_list.size())
  );

  log.printf(
    "maze> First group (LIGAND): from %d to %d.\n",
    ga_list[0].serial(),
    ga_list[ga_list.size()-1].serial()
  );

  if (gb_list.size() > 0) {
    log.printf(
      "maze> Second group (PROTEIN): from %d to %d.\n",
      gb_list[0].serial(),
      gb_list[gb_list.size()-1].serial()
    );
  }

  if (pbc_) {
    log.printf("maze> Using periodic boundary conditions.\n");
  }
  else {
    log.printf("maze> Without periodic boundary conditions.\n");
  }

  if (do_pair) {
    log.printf("maze> With PAIR option.\n");
  }

  if (do_neigh) {
    log.printf(
      "maze> Using neighbor lists updated every %d steps and cutoff %f.\n",
      nl_stride_,
      nl_cutoff_
    );
  }

  // OpenMP
  stride_ = comm.Get_size();
  rank_ = comm.Get_rank();

  n_threads_ = OpenMP::getNumThreads();
  unsigned int nn = neighbor_list_->size();

  if (n_threads_ * stride_ * 10 > nn) {
    n_threads_ = nn / stride_ / 10;
  }

  if (n_threads_ == 0) {
    n_threads_ = 1;
  }

  if (keywords.exists("OPTIMIZER_STRIDE")) {
    parse("OPTIMIZER_STRIDE", optimizer_stride_);

    plumed_massert(
      optimizer_stride_,
      "maze> OPTIMIZER_STRIDE should be explicitly specified and positive.\n"
    );

    log.printf(
      "maze> Launching optimization every %u steps.\n",
      optimizer_stride_
    );
  }

  rnd::randomize();

  opt_.zero();

  addComponentWithDerivatives("x");
  componentIsNotPeriodic("x");

  addComponentWithDerivatives("y");
  componentIsNotPeriodic("y");

  addComponentWithDerivatives("z");
  componentIsNotPeriodic("z");

  addComponent("loss");
  componentIsNotPeriodic("loss");

  addComponent("sr");
  componentIsNotPeriodic("sr");

  value_x_ = getPntrToComponent("x");
  value_y_ = getPntrToComponent("y");
  value_z_ = getPntrToComponent("z");
  value_action_ = getPntrToComponent("loss");
  value_sampling_radius_ = getPntrToComponent("sr");
}

double Optimizer::pairing(double distance) const {
  return loss_->pairing(distance);
}

Vector Optimizer::center_of_mass() const {
  const unsigned nl_size = neighbor_list_->size();

  Vector center_of_mass;
  center_of_mass.zero();
  double mass = 0;

  for (unsigned int i = 0; i < nl_size; ++i) {
    unsigned int i0 = neighbor_list_->getClosePair(i).first;
    center_of_mass += getPosition(i0) * getMass(i0);
    mass += getMass(i0);
  }

  return center_of_mass / mass;
}

void Optimizer::prepare() {
  if (neighbor_list_->getStride() > 0) {
    if (first_time_ || (getStep() % neighbor_list_->getStride() == 0)) {
      requestAtoms(neighbor_list_->getFullAtomList());

      validate_list_ = true;
      first_time_ = false;
    }
    else {
      requestAtoms(neighbor_list_->getReducedAtomList());

      validate_list_ = false;

      if (getExchangeStep()) {
        plumed_merror(
          "maze> Neighbor lists should be updated on exchange steps -- choose "
          "an NL_STRIDE which divides the exchange stride.\n");
      }
    }

    if (getExchangeStep()) {
      first_time_ = true;
    }
  }
}

double Optimizer::score() {
  const unsigned nl_size = neighbor_list_->size();
  Vector distance;
  double function = 0;

  #pragma omp parallel num_threads(n_threads_)
  {
    #pragma omp for reduction(+:function)
    for(unsigned int i = 0; i < nl_size; i++) {
      unsigned i0 = neighbor_list_->getClosePair(i).first;
      unsigned i1 = neighbor_list_->getClosePair(i).second;

      if (getAbsoluteIndex(i0) == getAbsoluteIndex(i1)) {
        continue;
      }

      if (pbc_) {
        distance = pbcDistance(getPosition(i0), getPosition(i1));
      }
      else {
        distance = delta(getPosition(i0), getPosition(i1));
      }

      function += pairing(distance.modulo());
    }
  }

  return function;
}

void Optimizer::update_nl() {
  if (neighbor_list_->getStride() > 0 && validate_list_) {
    neighbor_list_->update(getPositions());
  }
}

double Optimizer::sampling_radius()
{
  const unsigned nl_size=neighbor_list_->size();
  Vector d;
  double min=std::numeric_limits<int>::max();

  for (unsigned int i = 0; i < nl_size; ++i) {
    unsigned i0 = neighbor_list_->getClosePair(i).first;
    unsigned i1 = neighbor_list_->getClosePair(i).second;

    if (getAbsoluteIndex(i0) == getAbsoluteIndex(i1)) {
      continue;
    }

    if (pbc_) {
      d = pbcDistance(getPosition(i0), getPosition(i1));
    }
    else {
      d = delta(getPosition(i0), getPosition(i1));
    }

    double dist = d.modulo();

    if(dist < min) {
      min = dist;
    }
  }

  return min;
}

void Optimizer::calculate() {
  update_nl();

  if (getStep() % optimizer_stride_ == 0 && !first_step_) {
    optimize();

    value_x_->set(opt_[0]);
    value_y_->set(opt_[1]);
    value_z_->set(opt_[2]);

    value_action_->set(score());
    value_sampling_radius_->set(sampling_radius());
  }
  else {
    first_step_=false;

    value_x_->set(opt_[0]);
    value_y_->set(opt_[1]);
    value_z_->set(opt_[2]);

    value_action_->set(score());
    value_sampling_radius_->set(sampling_radius());
  }
}

} // namespace maze
} // namespace PLMD
