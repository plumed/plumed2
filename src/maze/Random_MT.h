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
#ifndef __PLUMED_maze_Random_MT_h
#define __PLUMED_maze_Random_MT_h

/**
 * @file Random_MT.h
 *
 * @author J. Rydzewski (jr@fizyka.umk.pl)
 */

#include "Core.h"

namespace PLMD {
namespace maze {

/**
 * @class rnd Random_MT.h "maze/Random_MT.h"
 *
 * @brief Mersenne Twister sampler for random variables.
 *
 * Supports generating integers, doubles, and std::vectors and PLMD::Vectors
 * within a given range.
 */
class rnd {
public:
  /**
   * Initialize MT sampler engine based on std::mt19937_64.
   */
  static std::mt19937_64& mt_eng();

  /**
   * Feed a random seed.
   */
  static void randomize();

  /**
   * Returns a random double from the Cauchy distribution.
   *
   * @param m mean
   * @param s spread
   */
  static double next_cauchy(double m, double s);

  /**
   * Returns a random int from the uniform distribution from a [f, e) range.
   *
   * @param f begin
   * @param e end
   */
  static int next_int(int f, int e);
  static int next_int(int e);

  /**
   * Returns a random double from the uniform distribution from a [f, e) range.
   */
  static double next_double(double f, double e);
  static double next_double(double e);
  static double next_double();

  /**
   * Returns a random vector<double> from the uniform distribution from a [f, e)
   * range of length n.
   */
  static std::vector<double> next_double(double f, double e, std::size_t n);

  /**
   * Returns a random PLMD::Vector of length r.
   */
  static Vector next_plmd_vector();
  static Vector next_plmd_vector(double r);

  /**
   * Returns a random std::vector of length r.
   */
  static std::vector<double> next_std_vector();
  static std::vector<double> next_std_vector(double r);
};

inline Vector rnd::next_plmd_vector(double r) {
  return next_plmd_vector() * r;
}

inline std::vector<double> rnd::next_std_vector(double r) {
  double t = next_double() * 2.0 * pi;
  double p = next_double() * pi;

  return std::vector<double> {
    r * std::sin(p) * std::cos(t),
    r * std::sin(t) * std::sin(p),
    r * std::cos(p)
  };
}

inline std::vector<double> rnd::next_double(double f, double e, size_t n) {
  std::vector<double> t(n);
  for (std::size_t i=0; i < n; ++i) {
    t[i] = next_double(f, e);
  }

  return t;
}

inline Vector rnd::next_plmd_vector() {
  double t = next_double() * 2.0 * pi;
  double p = next_double() * pi;

  return Vector(
           std::sin(p) * std::cos(t),
           std::sin(t) * std::sin(p),
           std::cos(p)
         );
}

inline std::vector<double> rnd::next_std_vector() {
  double t = next_double() * 2.0 * pi;
  double p = next_double() * pi;

  return std::vector<double> {
    std::sin(p) * std::cos(t),
    std::sin(t) * std::sin(p),
    std::cos(p)
  };
}

inline double rnd::next_double() {
  static std::uniform_real_distribution<double> dist_double(0, 1);
  std::uniform_real_distribution<double>::param_type p(0, 1);
  dist_double.param(p);

  return dist_double(mt_eng());
}

inline double rnd::next_double(double e) {
  return next_double(0, e);
}

inline double rnd::next_double(double f, double e) {
  static std::uniform_real_distribution<double> dist_double(f, e);
  std::uniform_real_distribution<double>::param_type p(f, e);
  dist_double.param(p);

  return dist_double(mt_eng());
}

inline int rnd::next_int(int e) {
  static std::uniform_int_distribution<int> dist_int(0, e-1);
  std::uniform_int_distribution<int>::param_type p(0, e-1);
  dist_int.param(p);

  return dist_int(mt_eng());
}

inline int rnd::next_int(int f, int e) {
  static std::uniform_int_distribution<int> dist_int(f, e-1);
  std::uniform_int_distribution<int>::param_type p(f, e-1);
  dist_int.param(p);

  return dist_int(mt_eng());
}

inline double rnd::next_cauchy(double m, double s) {
  static std::cauchy_distribution<double> dist_cauchy(m, s);

  return dist_cauchy(mt_eng());
}

inline std::mt19937_64& rnd::mt_eng() {
  static std::mt19937_64 mt{};

  return mt;
}

inline void rnd::randomize() {
  mt_eng().seed(time(0));
}

} // namespace maze
} // namespace PLMD

#endif // __PLUMED_maze_Random_MT_h
