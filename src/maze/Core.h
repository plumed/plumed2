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
#ifndef __PLUMED_maze_Core_h
#define __PLUMED_maze_Core_h

/**
 * @file Core.h
 *
 * @author J. Rydzewski (jr@fizyka.umk.pl)
 *
 * @brief Header with needed includes from std and PLMD.
 */

namespace PLMD {
namespace maze {
/* Nothing to do */
} // namespace maze
} // namespace PLMD

// std
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <utility>
#include <string>
#include <random>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <limits>
#include <numeric>
#include <algorithm>
#include <functional>
#include <memory>

// PLMD
#include "tools/Vector.h"
#include "tools/Tools.h"

// maze
#include "Random_MT.h"
#include "Tools.h"

#endif // __PLUMED_maze_Core_h
