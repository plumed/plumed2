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
 * @file Loss.cpp
 * @author J. Rydzewski (jr@fizyka.umk.pl)
 */

#include "Loss.h"

namespace PLMD {
namespace maze {

//+PLUMEDOC MAZE_LOSS MAZE_LOSS
/*

Define a coarse-grained loss function describing interactions in a
ligand-protein complex, which is minimized during the simulation to
obtain ligand unbinding pathways.

The loss function is the following:
\f[
\mathcal{L}=
\sum_{i=1}^{N_p}
r_i^{-\alpha}\text{e}^{-\beta r_i^{-\gamma}},
\f]
where \f$N_p\f$ is the number of ligand-protein atom pairs, \f$r\f$
is a re-scaled distance between the \f$i\f$th pair, and \f$\alpha,
\beta, \gamma\f$ are the positive parameters defined in that order by
the PARAMS keyword.

\par Examples

The loss function can be defined in the following way:
\plumedfile
l: MAZE_LOSS PARAMS=1,1,1
\endplumedfile

*/
//+ENDPLUMEDOC

// Registers the LOSS action.
PLUMED_REGISTER_ACTION(Loss, "MAZE_LOSS")

void Loss::registerKeywords(Keywords& keys) {
  Colvar::registerKeywords(keys);

  keys.add(
    "compulsory",
    "PARAMS",
    "Parameters for the loss function."
  );
}

Loss::Loss(const ActionOptions& ao)
  : PLUMED_COLVAR_INIT(ao)
{
  if (keywords.exists("PARAMS")) {
    parseVector("PARAMS", params_);

    plumed_massert(
      params_.size() == 3,
      "maze> PARAMS should be of size 3: alpha, beta, gamma\n"
    );

    plumed_massert(
      params_[0] > 0 && params_[1] > 0 && params_[2] > 0,
      "maze> Each parameter should be positive\n"
    );

    log.printf("maze> \t Loss parsed with parameters: ");
    for (size_t i = 0; i < params_.size(); ++i) {
      log.printf("%f ", params_[i]);
    }
    log.printf("\n");
  }

  checkRead();
}

double Loss::pairing(double distance) {
  double alpha = params_[0];
  double beta = params_[1];
  double gamma = params_[2];

  if (atoms.getUnits().getLengthString() == "nm") {
    distance *= 10.0;
  }

  return pow(distance, -alpha) * exp(-beta * pow(distance, gamma));
}

} // namespace maze
} // namespace PLMD
