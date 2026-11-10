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
 * @file Rnd.cpp
 *
 * @author J. Rydzewski (jr@fizyka.umk.pl)
 */

#include "Rnd.h"
#include "tools/Tools.h"

#include <cmath>

namespace PLMD {
namespace maze {

double rnd::next_cauchy(double m, double s) {
  return m + s * std::tan(pi * (generator_.U01() - 0.5));
}

std::vector<double> rnd::next_double(double f, double e, std::size_t n) {
  std::vector<double> t(n);
  for (std::size_t i=0; i < n; ++i) {
    t[i] = next_double(f, e);
  }
  return t;
}

Vector rnd::next_plmd_vector() {
  // cos(theta) uniform in (-1,1): uniform on the sphere
  double t = next_double() * 2.0 * pi;
  double c = 2.0 * next_double() - 1.0;
  double s = std::sqrt(1.0 - c * c);
  return Vector(
           s * std::cos(t),
           s * std::sin(t),
           c
         );
}

} // namespace maze
} // namespace PLMD
