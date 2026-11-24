/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2023 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

#include "plumed/tools/Communicator.h"
#include <fstream>
#include <stdexcept>

int main() {
  // Do not initialize MPI: exercise the serial fallback in every build mode.
  PLMD::Communicator comm;
  const int input[2] = {11, 13};
  int output[4] = {-1, -1, -1, -1};
  const int counts[1] = {2};
  const int displacements[1] = {1};
  // The pointer overload cannot know the receive buffer extent. Its count is
  // supplied by counts, and the displacement is measured in datatype units.
  comm.Allgatherv(input, 2, output, counts, displacements);
  if(output[0]!=-1 || output[1]!=11 || output[2]!=13 || output[3]!=-1) {
    throw std::runtime_error("serial Allgatherv copied the wrong buffer region");
  }
  std::ofstream("result") << "OK\n";
}
