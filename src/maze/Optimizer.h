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
#ifndef __PLUMED_maze_Optimizer_h
#define __PLUMED_maze_Optimizer_h

/**
 * @file Optimizer.h
 *
 * @author J. Rydzewski (jr@fizyka.umk.pl)
 */

#include "colvar/Colvar.h"
#include "tools/Communicator.h"
#include "tools/OpenMP.h"
#include "tools/NeighborList.h"
#include "tools/Vector.h"

#include "Core.h"
#include "Loss.h"

#define PLUMED_OPT_INIT(ao) Action(ao), Optimizer(ao)

namespace PLMD {
namespace maze {

/**
 * @ingroup INHERIT
 *
 * @class Optimizer Optimizer.h "maze/Optimizer.h"
 *
 * @brief Base class for implementing optimizers for ligand unbinding.
 *
 * An optimizer is defined as a colvar that can be passed to Optimizer_Bias.
 */
class Optimizer: public colvar::Colvar {
public:
  /**
   * PLMD constructor.
   *
   * @param[in] ao PLMD::ActionOptions&
   */
  explicit Optimizer(const ActionOptions&);

  /**
   * Destructor.
   */
  ~Optimizer() { /* Nothing to do. */ }

  /**
   * Registers PLMD keywords.
   *
   * @param[in] keys PLMD keywords
   */
  static void registerKeywords(Keywords& keys);

  /**
   * The pairing function needs to be overridden by a specific optimizer.
   *
   * @param[in] distance distance between a pair of atoms
   */
  virtual double pairing(double distance) const;

  /**
   * Optimal values needed for biasing are computed by methods overridding the
   * optimize function.
   */
  virtual void optimize() = 0;

  /**
   * Calculate the optimal direction of pulling.
   */
  void calculate();

  /**
   * Prepare the neighbor list.
   */
  void prepare();

  /**
   * Score a ligand-protein configuration.
   *
   * @return score
   */
  double score();

  /**
   * Calculate sampling radius as the minimal distance between two groups in
   * neighbors list.
   *
   * @return minimal distance of ligand-protein atom pairs
   */
  double sampling_radius();

  /**
   * Load new positions of atoms in the neighbor list.
   */
  void update_nl();

  /**
   * Calculate the center of mass.
   *
   * @return center of mass
   */
  Vector center_of_mass() const;

public:
  /**
   * Getters and setters.
   */

  std::string get_label() const;
  void set_label(const std::string&);

  // Start optimizer at time = 0.
  void start_step_0();

  // Start optimizer at time = optimizer stride.
  void start_step_stride();

  Vector get_opt() const;
  void set_opt(Vector);

  double get_opt_value() const;
  void set_opt_value(double);

  unsigned int get_optimizer_stride() const;
  void set_optimizer_stride(unsigned int);

  bool is_pbc_on() const;
  void pbc_on();
  void pbc_off();

  unsigned int get_n_iterations() const;
  void set_n_iterations(unsigned int);

  double get_sampling_radius() const;
  void set_sampling_radius(double);

  unsigned int get_rank_openmp() const;
  void set_rank_openmp(unsigned int);

  unsigned int get_stride_openmp() const;
  void set_stride_openmp(unsigned int);

  unsigned int get_n_threads_openmp() const;
  void set_n_threads_openmp(unsigned int);

  unsigned int get_nl_stride() const;
  void set_nl_stride(unsigned int);

  double get_nl_cutofff() const;
  void set_nl_cutoff(double);

protected:
  //! Optimizer label.
  std::string label_;

  //! Start either at time =  0 or time = optimizer stride.
  bool first_step_;

  //! Biasing direction.
  Vector opt_;

  //! Current loss function value.
  double opt_value_;

  //! Optimizer stride.
  unsigned int optimizer_stride_;

  //! Periodic boundary conditions.
  bool pbc_;

  //! Number of global iterations.
  unsigned int n_iter_;

  //! Sampling radius.
  double sampling_r_;

  /**
   * OpenMP
   */
  unsigned int rank_;
  unsigned int stride_;
  unsigned int n_threads_;

  //! Neighbor list of ligand-protein atom pairs.
  NeighborList *neighbor_list_;

  //! Neighbor list cut-off.
  double nl_cutoff_;

  //! Neighbor list stride.
  int nl_stride_;

private:
  bool serial_;
  bool validate_list_;
  bool first_time_;

  //! Pointer to the loss function.
  Loss* loss_;
  std::vector<Loss*> vec_loss_;

public:
  /*
   * Pointers to PLMD components.
   */

  //! Biased cv.
  Value* value_x_;
  Value* value_y_;
  Value* value_z_;

  //! Loss value.
  Value* value_action_;
  //! Sampling radiues value.
  Value* value_sampling_radius_;
};

/*
 * Getters and setters.
 */

inline void Optimizer::set_nl_cutoff(double nl_cutoff) {
  nl_cutoff_=nl_cutoff;
}

inline double Optimizer::get_nl_cutofff() const {
  return nl_cutoff_;
}

inline void Optimizer::set_nl_stride(unsigned int nl_stride) {
  nl_stride_=nl_stride;
}

inline unsigned int Optimizer::get_nl_stride() const {
  return nl_stride_;
}

inline void Optimizer::set_n_threads_openmp(unsigned int n_threads) {
  n_threads_=n_threads;
}

inline unsigned int Optimizer::get_n_threads_openmp() const {
  return n_threads_;
}

inline void Optimizer::set_stride_openmp(unsigned int stride) {
  stride_=stride;
}

inline unsigned int Optimizer::get_stride_openmp() const {
  return stride_;
}

inline void Optimizer::set_rank_openmp(unsigned int rank) {
  rank_=rank;
}

inline unsigned int Optimizer::get_rank_openmp() const {
  return rank_;
}

inline void Optimizer::set_sampling_radius(double sampling_r) {
  sampling_r_=sampling_r;
}

inline double Optimizer::get_sampling_radius() const {
  return sampling_r_;
}

inline void Optimizer::set_n_iterations(unsigned int n_iter) {
  n_iter_=n_iter;
}

inline unsigned int Optimizer::get_n_iterations() const {
  return n_iter_;
}

inline void Optimizer::pbc_off() {
  pbc_=false;
}

inline void Optimizer::pbc_on() {
  pbc_=true;
}

inline bool Optimizer::is_pbc_on() const {
  return pbc_==true;
}

inline void Optimizer::set_optimizer_stride(unsigned int optimizer_stride) {
  optimizer_stride_=optimizer_stride;
}

inline unsigned int Optimizer::get_optimizer_stride() const {
  return optimizer_stride_;
}

inline void Optimizer::set_opt_value(double opt_value) {
  opt_value_=opt_value;
}

inline double Optimizer::get_opt_value() const {
  return opt_value_;
}

inline void Optimizer::set_opt(Vector opt) {
  opt_=opt;
}

inline Vector Optimizer::get_opt() const {
  return opt_;
}

inline void Optimizer::set_label(const std::string& label) {
  label_=label;
}

inline std::string Optimizer::get_label() const {
  return label_;
}

inline void Optimizer::start_step_0() {
  first_step_=false;
}

inline void Optimizer::start_step_stride() {
  first_step_=true;
}

} // namespace maze
} // namespace PLMD

#endif // __PLUMED_maze_Optimizer_h
