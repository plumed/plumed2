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
#ifndef __PLUMED_maze_Memetic_h
#define __PLUMED_maze_Memetic_h

/**
 * @file Memetic.h
 *
 * @author J. Rydzewski (jr@fizyka.umk.pl)
 */

#include "core/ActionRegister.h"

#include "Core.h"
#include "Member.h"
#include "Optimizer.h"

namespace PLMD {
namespace maze {

/**
 * @class Memetic Memetic.h "maze/Memetic.h"
 *
 * @brief Memetic algorithms for the optimization of the loss function.
 */
class Memetic: public Optimizer {
public:
  /**
   * PLMD constructor.
   *
   * @param[in] ao PLMD::ActionOptions&.
   */
  explicit Memetic(const ActionOptions& ao);

  /**
   * Destructor required for deriving classes.
   */
  ~Memetic();

  /**
   * Registers PLMD keywords.
   *
   * @param[in] keys Keywords.
   */
  static void registerKeywords(Keywords& keys);

  /**
   * Each class deriving from Optimizer needs to override this function.
   */
  void optimize() override;

private:
  /**
   * Create a set of translations relative to the ligand, each translation
   * encodes a probable biasing direction.
   */
  void initialize_members();

  /**
   * Score each translation by the loss function.
   */
  void score_members();

  /**
   * Calculate the mean score.
   */
  double score_mean();

  /**
   * Sort the population using heaps, required for finding a minimum of the loss
   * function.
   */
  void sort_members();

  /**
   * Encode a ligand conformation.
   */
  Vector create_coding();

  /**
   * Check if the vector length is out of bounds.
   *
   * @param[in] v Vector length.
   */
  bool out_of_bounds(double v);

  /**
   * Score a single member.
   *
   * @param[in] v Member's translation.
   */
  double score_member(const Vector& v);

  /**
   * Print a status.
   */
  void print_status() const;

  /**
   * Select a new population using the roulette selection.
   */
  void selection_roulette();

  /**
   * Perform mating in the population.
   *
   * @param[in,out] members Population.
   */
  void mating(std::vector<Member>& members);

  /**
   * Mate two members.
   *
   * @param[in,out] m1 1st member.
   * @param[in,out] m2 2nd member.
   */
  void crossover(Member& m1, Member& m2);

  /**
   * Mutate a member.
   *
   * @param[in,out] m Member.
   */
  void mutation(Member& m);

  /**
   * Mutate the population.
   *
   * @param[in,out] members Population.
   */
  void mutation(std::vector<Member>& members);

protected:
  /**
   * Local searches to improve global solutions.
   */

  /**
   * Stochastic hill climbing.
   *
   * Neighbors with better or equal cost should be accepted, allowing the
   * technique to navigate across plateaus in the response surface.
   *
   * @param[in,out] m Member.
   * @param[in] params None.
   */
  void stochastic_hill_climbing(
    Member& m,
    /* none */ const std::vector<double>& params = {}
  );

  /**
   * Random-restart hill climbing.
   *
   * The algorithm can be restarted and repeated a number of  times after it
   * converges to provide an improved result.
   *
   * @param[in,out] m Member.
   * @param[in] params Number of restarts.
   */
  void random_restart_hill_climbing(
    Member& m,
    /* n_restarts */ const std::vector<double>& params = {10}
  );

  /**
   * Solis-Wets random walk.
   *
   * Adaptive random search algorithm was designed to address the limitations of
   * the fixed step size. The strategy for adaptive random search is to
   * continually approximate the optimal step size required to reach the global
   * optimum in the search space. This is achieved by trialling and adopting
   * smaller or larger step sizes only if they result in an improvement in the
   * search performance.
   *
   * Very large step sizes are trialled with a much lower frequency. This
   * strategy of preferring large moves is intended to allow the technique to
   * escape local optima. Smaller step sizes are adopted if no improvement is
   * made for an extended period.
   *
   * @param[in,out] m Member.
   * @param[in] params
   */
  void adaptive_random_search(
    Member& m,
    /* */ const std::vector<double>& params = {1.0, 1.e-5, 2.0, 2.0, 3.0, 3.0}
  );

  /**
   * Luus–Jaakola heuristics.
   *
   * @param[in,out] m Member.
   * @param[in] params Bounds.
   */
  void luus_jaakola(
    Member& m,
    /* bounds */ const std::vector<double>& params
  );

  /**
   * Local annealing.
   *
   * @param[in,out] m Member.
   * @param[in] params T, alpha.
   */
  void annealing(
    Member& m,
    /* T, alpha */const std::vector<double>& params = {300.0, 0.95}
  );

  /**
   * Apply local search to members.
   *
   * @param[in,out] members Population.
   */
  void local_search(std::vector<Member>& members);

protected:
  /**
   * Return an optimal biasing direction.
   */
  Vector solve();

public:
  /**
   * Setters and getters.
   */

  unsigned int get_capacity() const;
  void set_capacity(unsigned int);

  unsigned int get_coding_len() const;
  void set_coding_len(unsigned int);

  unsigned int get_n_local_iterations() const;
  void set_n_local_iterations(unsigned int);

  unsigned int get_n_global_iterations() const;
  void set_n_global_iterations(unsigned int);

  double get_mutation_rate() const;
  void set_mutation_rate(double);

  double get_mating_rate() const;
  void set_mating_rate(double);

  double get_cauchy_mean() const;
  void set_cauchy_mean(double);

  double get_cauchy_spread() const;
  void set_cauchy_spread(double);

  bool is_local_search_on() const;
  void local_search_on();
  void local_search_off();

  double get_local_search_rate() const;
  void set_local_search_rate(double);

  std::string get_local_search_type() const;
  void set_local_search_type(const std::string&);

protected:
  //! Population
  std::vector<Member> members_;

  //! Bound
  double bound_;

  //! Scores
  double score_worst_;
  double score_best_;
  Member member_best_;

  //! If a local search is performed
  bool adaptation_;

protected:
  //! Size of population
  unsigned int capacity_;
  //! Length of coding
  unsigned int coding_len_;

  //! Number of local search iterations
  unsigned int n_local_iterations_;
  //! Number of global search iterations, doomsday
  unsigned int n_global_iterations_;

  //! Probability of mutation
  double mutation_rate_;
  //! Probability of mating
  double mating_rate_;

  //! Mean and spread of cauchy sampler
  double cauchy_mean_alpha_;
  double cauchy_mean_beta_;

  //! If local search is employed in sampling
  bool local_search_on_;
  //! Rate of local mutation
  double local_search_rate_;
  //! Type of local search, stochastic_hill_climbing or adaptive_random_search
  std::string local_search_type_;
};

void Memetic::initialize_members() {
  members_.clear();
  members_.resize(capacity_);

  for (size_t i = 0; i < capacity_; ++i) {
    Member m{};

    m.score=0;
    m.translation=create_coding();

    members_.at(i) = m;
  }
}

void Memetic::score_members() {
  for (size_t i = 0; i < members_.size(); ++i) {
    double s = score_member(members_[i].translation);
    members_[i].score=s;
  }
}

void Memetic::sort_members() {
  std::make_heap(
    members_.begin(),
    members_.end(),
    compare
  );

  std::sort_heap(
    members_.begin(),
    members_.end(),
    compare
  );

  member_best_ = members_[capacity_ - 1];
  score_best_ = members_[capacity_ - 1].score;
  score_worst_ = members_[0].score;
}

double Memetic::score_mean() {
  auto acc = [](double s, const Member& m) { return s + m.score; };

  return std::accumulate(
           members_.begin(),
           members_.end(),
           0.0,
           acc) / capacity_;
}

void Memetic::selection_roulette() {
  std::vector<Member> sel(members_);
  std::vector<double> rel_scores(capacity_, 0.0);

  for (std::size_t i = 0; i < capacity_; ++i) {
    double r = 1.0 / (members_[i].score + 0.01);
    rel_scores.at(i) = r;
  }

  std::vector<double> cum_sum(capacity_, 0.0);
  std::partial_sum(
    rel_scores.begin(),
    rel_scores.end(),
    cum_sum.begin(),
    std::plus<double>()
  );

  double sum = cum_sum.back();
  members_.clear();
  members_.resize(capacity_);
  int chosen = -1;

  for (size_t j = 0; j < capacity_; ++j) {
    double probability=rnd::next_double(sum);
    for (size_t i = 0; i < cum_sum.size(); ++i) {
      if (cum_sum[i] > probability) {
        chosen = i;

        members_.at(j).score = sel.at(chosen).score;
        members_.at(j).translation = sel.at(chosen).translation;

        break;
      }
    }
  }
}

void Memetic::crossover(Member& s1, Member& s2) {
  size_t i = rnd::next_int(1, coding_len_ - 1);

  Member z1(s1);
  Member z2(s2);

  for (size_t j = i; j < coding_len_; ++j) {
    z1.translation[j] = s2.translation[j];
    z2.translation[j] = s1.translation[j];
  }

  if (!out_of_bounds(z1.translation.modulo()) && !out_of_bounds(z2.translation.modulo())) {
    s1 = z1;
    s2 = z2;
  }
}

void Memetic::mutation(Member& m) {
  int which = rnd::next_int(coding_len_);
  double v = rnd::next_cauchy(cauchy_mean_alpha_, cauchy_mean_beta_);
  m.translation[which] += v;
  if (out_of_bounds(m.translation.modulo())) {
    m.translation[which] -= v;
  }
}

void Memetic::mutation(std::vector<Member>& m) {
  for (std::vector<Member>::iterator it = m.begin(), end = m.end(); it != end; ++it) {
    double r = rnd::next_double();
    if (r < mutation_rate_) {
      mutation(*it);
    }
  }
}

void Memetic::stochastic_hill_climbing(
  Member& m,
  const std::vector<double>& params)
{
  for (std::size_t i = 0; i < n_local_iterations_; ++i) {
    Member n;
    n.translation = m.translation;
    mutation(n);
    double score_n = score_member(n.translation);

    if (m.score > score_n) {
      m.translation = n.translation;
    }
  }
}

void Memetic::random_restart_hill_climbing(
  Member& m,
  const std::vector<double>& params)
{
  unsigned int n_restarts = static_cast<int>(params[0]);
  std::vector<Member> s(n_restarts);
  unsigned int ndx = 0;
  double min = m.score;

  for (std::size_t r = 0; r < n_restarts; ++r) {
    Member n = m;
    stochastic_hill_climbing(n);
    s[r] = n;

    if (min > n.score) {
      min = n.score;
      ndx = r;
    }
  }

  m = s[ndx];
}

void Memetic::annealing(
  Member& m,
  const std::vector<double>& params)
{
  double T = params[0];
  double alpha = params[1];

  for (std::size_t i = 0; i < n_local_iterations_; ++i) {
    Member n = m;
    mutation(n);
    double score_n = score_member(n.translation);

    double probability = std::min(
                           1.0,
                           std::exp(-(score_n - m.score) / T)
                         );

    double r = rnd::next_double();

    if (r < probability) {
      m = n;
    }

    T *= alpha;
  }
}

void Memetic::luus_jaakola(
  Member& m,
  const std::vector<double>& params)
{
  /* TODO */
}

void Memetic::adaptive_random_search(
  Member& m,
  const std::vector<double>& params)
{
  double rho_start = params[0];
  double rho_lower_bound = params[1];
  double ex = params[2];
  double ct = params[3];
  int s_ex = static_cast<int>(params[4]);
  int f_ct = static_cast<int>(params[5]);

  unsigned int k = 0;
  Vector xk = m.translation;
  int ns = 0;
  int nf = 0;
  Vector bk;
  bk.zero();
  double rho_k = rho_start;

  while (rho_k > rho_lower_bound && k < n_local_iterations_) {
    if (ns >= s_ex) {
      rho_k *= ex;
    }
    else if (nf >= f_ct) {
      rho_k *= ct;
    }

    std::vector<double> chiv = rnd::next_double(-rho_k, rho_k, 3);
    Vector chi = tls::vector2Vector(chiv);
    Vector tmp;

    for (unsigned int i = 0; i < 3; ++i) {
      tmp[i] = 2.0 * (xk[i] - chi[i]);
    }

    double score_chi = score_member(chi);
    double score_xk = score_member(xk);
    double score_tmp = score_member(tmp);

    if (score_chi < score_xk) {
      ns++;
      nf = 0;

      for (unsigned int i = 0; i < 3; i++) {
        bk[i] = 0.2 * bk[i] + 0.4 * (chi[i] - xk[i]);
        xk[i] = chi[i];
      }
    }
    else if (score_tmp < score_xk && score_xk <= score_chi) {
      ns++;
      nf = 0;

      for (unsigned int i = 0; i < 3; i++) {
        bk[i] = bk[i] - 0.4 * (chi[i] - xk[i]);
        xk[i] = 2.0 * xk[i] - chi[i];
      }
    }
    else {
      ns = 0;
      nf++;

      for (unsigned int i = 0; i < 3; i++) {
        bk[i] = 0.5 * bk[i];
      }
    }

    k++;
  }

  m.translation = xk;
}

void Memetic::local_search(std::vector<Member>& m) {
  adaptation_ = true;

  if (local_search_on_) {
    for (std::size_t i = 0; i < capacity_; ++i) {
      double probability = rnd::next_double();

      if (probability < local_search_rate_) {
        if (local_search_type_ == "stochastic_hill_climbing")
          stochastic_hill_climbing(m[i]);
        else if (local_search_type_ == "adaptive_random_search")
          adaptive_random_search(m[i]);
        else if (local_search_type_ == "random_restart_hill_climbing")
          random_restart_hill_climbing(m[i]);
      }
    }
  }

  adaptation_ = false;
}

void Memetic::mating(std::vector<Member>& m) {
  std::vector<Member> offspring;

  while (m.size() != 0) {
    int j = rnd::next_int(m.size());
    int i = rnd::next_int(m.size());

    while (i == j) {
      j=rnd::next_int(m.size());
    }

    Member m1 = m[i];
    Member m2 = m[j];

    if (i > j) {
      m.erase(m.begin() + i);
      m.erase(m.begin() + j);
    }
    else if (j > i) {
      m.erase(m.begin() + j);
      m.erase(m.begin() + i);
    }

    double probability = rnd::next_double();
    if (probability < mating_rate_) {
      crossover(m1, m2);
    }

    offspring.push_back(m1);
    offspring.push_back(m2);
  }

  m = offspring;
  offspring.clear();
}

Vector Memetic::create_coding() {
  double s = sampling_radius();
  double r = rnd::next_double(s);

  return rnd::next_plmd_vector(r);
}

bool Memetic::out_of_bounds(double v) {
  double s = sampling_radius();

  return v > s;
}

double Memetic::score_member(const Vector& coding) {
  double action = 0;
  Vector distance;
  const unsigned nl_size = neighbor_list_->size();
  Vector dev = coding;

  #pragma omp parallel num_threads(get_n_threads_openmp())
  {
    #pragma omp for reduction(+:action)
    for (unsigned int i = 0; i < nl_size; i++) {
      unsigned i0 = neighbor_list_->getClosePair(i).first;
      unsigned i1 = neighbor_list_->getClosePair(i).second;

      if (getAbsoluteIndex(i0) == getAbsoluteIndex(i1)) {
        continue;
      }

      if (pbc_) {
        distance = pbcDistance(getPosition(i0) + dev, getPosition(i1));
      }
      else {
        distance = delta(getPosition(i0) + dev, getPosition(i1));
      }

      action += pairing(distance.modulo());
    }
  }

  return action;
}

void Memetic::print_status() const {
  log.printf("Lowest score: %f \n", score_best_);
}

Vector Memetic::solve() {
  initialize_members();

  size_t epoch = 0;
  while (epoch < n_global_iterations_) {
    score_members();

    selection_roulette();
    mating(members_);
    mutation(members_);
    local_search(members_);

    sort_members();

    epoch++;
  }

  return member_best_.translation / member_best_.translation.modulo();
}

inline unsigned int Memetic::get_capacity() const {
  return capacity_;
}

inline void Memetic::set_capacity(unsigned int capacity) {
  capacity_ = capacity;
}

inline unsigned int Memetic::get_coding_len() const {
  return coding_len_;
}

inline void Memetic::set_coding_len(unsigned int coding_len) {
  coding_len_ = coding_len;
}

inline unsigned int Memetic::get_n_global_iterations() const {
  return n_global_iterations_;
}

inline void Memetic::set_n_global_iterations(unsigned int n_global_iterations) {
  n_global_iterations_ = n_global_iterations;
}

inline unsigned int Memetic::get_n_local_iterations() const {
  return n_local_iterations_;
}

inline void Memetic::set_n_local_iterations(unsigned int n_local_iterations) {
  n_local_iterations_ = n_local_iterations;
}

inline double Memetic::get_mating_rate() const {
  return mating_rate_;
}

inline void Memetic::set_mating_rate(double mating_rate) {
  mating_rate_ = mating_rate;
}

inline double Memetic::get_mutation_rate() const {
  return mutation_rate_;
}

inline void Memetic::set_mutation_rate(double mutation_rate) {
  mutation_rate_ = mutation_rate;
}

inline double Memetic::get_cauchy_mean() const {
  return cauchy_mean_alpha_;
}

inline void Memetic::set_cauchy_mean(double cauchy_mean_alpha) {
  cauchy_mean_alpha_ = cauchy_mean_alpha;
}

inline double Memetic::get_cauchy_spread() const {
  return cauchy_mean_beta_;
}

inline void Memetic::set_cauchy_spread(double cauchy_mean_beta) {
  cauchy_mean_beta_ = cauchy_mean_beta;
}

inline bool Memetic::is_local_search_on() const {
  return local_search_on_;
}

inline void Memetic::local_search_on() {
  local_search_on_ = true;
}

inline void Memetic::local_search_off() {
  local_search_on_ = false;
}

inline double Memetic::get_local_search_rate() const {
  return local_search_rate_;
}

inline void Memetic::set_local_search_rate(double local_search_rate) {
  local_search_rate_ = local_search_rate;
}

inline std::string Memetic::get_local_search_type() const {
  return local_search_type_;
}

inline void Memetic::set_local_search_type(const std::string& local_search_type) {
  local_search_type_ = local_search_type;
}

} // namespace maze
} // namespace PLMD

#endif // __PLUMED_maze_Memetic_h
