/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2026 The plumed team
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
#ifndef __PLUMED_tools_AtomDistributionFiles_h
#define __PLUMED_tools_AtomDistributionFiles_h
#include "AtomDistribution.h"

#include "Vector.h"
#include "View.h"
#include "Random.h"
#include "TrajectoryParser.h"

#include <string_view>
#include <vector>

namespace PLMD {

/// atomic distribution from a trajectory file
class fileTraj:public AtomDistribution {
  TrajectoryParser parser;
  std::vector<double> masses{};
  std::vector<double> charges{};
  std::vector<Vector> coordinates{};
  std::vector<double> cell{0.0,0.0,0.0,
        0.0,0.0,0.0,
        0.0,0.0,0.0};
  bool read=false;
  bool dont_read_pbc=false;
  void rewind();
  //read the next step
  void step(bool doRewind=true);
public:
  static constexpr auto id="file";
  void frame(View<Vector> posToUpdate,
             View<double,9> box,
             unsigned /*step*/,
             Random& /*rng*/) override;

  fileTraj(std::string_view fmt,
           std::string_view fname,
           bool useMolfile,
           int command_line_natoms);
  bool overrideNat(unsigned& natoms) override;
};
} //namespace PLMD
#endif // __PLUMED_tools_AtomDistributionFiles_h
