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
 * @file Random_MT.cpp
 *
 * @author J. Rydzewski (jr@fizyka.umk.pl)
 *
 * @brief Dummy cpp file.
 */

#include "Random_MT.h"

namespace PLMD {
namespace maze {

std::mt19937_64& rnd::mt_eng() {
  static std::mt19937_64 mt{};

  return mt;
}

double rnd::next_double(double f, double e) {
  static std::uniform_real_distribution<double> dist_double(f, e);
  std::uniform_real_distribution<double>::param_type p(f, e);
  dist_double.param(p);

  return dist_double(mt_eng());
}

double rnd::next_double() {
  static std::uniform_real_distribution<double> dist_double(0, 1);
  std::uniform_real_distribution<double>::param_type p(0, 1);
  dist_double.param(p);

  return dist_double(mt_eng());
}

int rnd::next_int(int e) {
  static std::uniform_int_distribution<int> dist_int(0, e-1);
  std::uniform_int_distribution<int>::param_type p(0, e-1);
  dist_int.param(p);

  return dist_int(mt_eng());
}

int rnd::next_int(int f, int e) {
  static std::uniform_int_distribution<int> dist_int(f, e-1);
  std::uniform_int_distribution<int>::param_type p(f, e-1);
  dist_int.param(p);

  return dist_int(mt_eng());
}

double rnd::next_cauchy(double m, double s) {
  static std::cauchy_distribution<double> dist_cauchy(m, s);

  return dist_cauchy(mt_eng());
}

} // namespace maze
} // namespace PLMD
