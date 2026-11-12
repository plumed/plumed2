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
#ifndef __PLUMED_maze_Rnd_h
#define __PLUMED_maze_Rnd_h

/**
 * @file Rnd.h
 *
 * @author J. Rydzewski (jr@fizyka.umk.pl)
 */

#include "tools/Random.h"
#include "tools/Vector.h"

#include <vector>

namespace PLMD {
namespace maze {

/**
 * @class rnd Rnd.h "maze/Rnd.h"
 *
 * @brief Random numbers for the maze optimizers, drawn from PLMD::Random.
 *
 * Each optimizer owns its own generator, so that runs with a given SEED are reproducible.
 */
class rnd {
  PLMD::Random generator_;

public:
  /**
   * Seed the generator.
   */
  void seed(int s) {
    generator_.setSeed(s);
  }

  /**
   * Returns a random double from the Cauchy distribution.
   *
   * @param m mean
   * @param s spread
   */
  double next_cauchy(double m, double s);

  /**
   * Returns a random int from the uniform distribution from [f, e) interval.
   */
  int next_int(int f, int e) {
    return f + generator_.RandInt(e - f);
  }

  /**
   * Returns a random int from the uniform distribution from [0, e) interval.
   */
  int next_int(int e) {
    return generator_.RandInt(e);
  }

  /**
   * Returns a random double from the uniform distribution from (f, e) interval.
   */
  double next_double(double f, double e) {
    return f + (e - f) * generator_.U01();
  }

  /**
   * Returns a random double from the uniform distribution from (0, e) interval.
   */
  double next_double(double e) {
    return e * generator_.U01();
  }

  /**
   * Returns a random double from the uniform distribution from (0, 1) interval.
   */
  double next_double() {
    return generator_.U01();
  }

  /**
   * Returns a vector of n random doubles from the uniform distribution from (f, e) interval.
   */
  std::vector<double> next_double(double f, double e, std::size_t n);

  /**
   * Returns a random vector of length one, uniformly distributed on the sphere
   * (the cosine of the polar angle is uniform in (-1,1)).
   */
  Vector next_plmd_vector();

  /**
   * Returns a random vector of length r.
   */
  Vector next_plmd_vector(double r) {
    return next_plmd_vector() * r;
  }
};

} // namespace maze
} // namespace PLMD

#endif // __PLUMED_maze_Rnd_h
