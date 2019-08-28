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
#ifndef __PLUMED_maze_Member_h
#define __PLUMED_maze_Member_h

/**
 * @file Member.h
 *
 * @author J. Rydzewski (jr@fizyka.umk.pl)
 */

#include "Core.h"

namespace PLMD {
namespace maze {

/**
 * @class Member Member.h "maze/Member.h"
 *
 * @brief Defines the encoding for a ligand conformation.
 *
 * The Member class is required by some optimizers. Ligand conformations are
 * encoded by a Cartesian translation relative to a ligand conformation from
 * the MD simulation. Each translation has a loss (score) which tells how
 * prefered is the encoding w.r.t. a chosen loss function.
 */
struct Member {
public:
  Member()
    : score(-1),
      translation({0, 0, 0}) { /* Nothing to do. */ }

  //! Value of the loss function.
  double score;

  //! Encoding of a ligand conformation, i.e., its translation.
  Vector translation;
};

inline bool compare(Member& m, Member& n) {
  return m.score > n.score;
}

} // namespace maze
} // namespace PLMD

#endif // __PLUMED_maze_Member_h
