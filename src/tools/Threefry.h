/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2023 The plumed team
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
#ifndef __PLUMED_tools_Threefry_h
#define __PLUMED_tools_Threefry_h

#include <array>
#include <cstdint>

namespace PLMD {

/// \ingroup TOOLBOX
/// Counter-based random numbers: Threefry2x64-20 (Salmon, Moraes, Dror and Shaw,
/// "Parallel random numbers: as easy as 1, 2, 3", SC'11; the Random123 library).
///
/// Each 128-bit output block is a function of a 128-bit key and a 128-bit counter
/// only, and any keys and counters (also consecutive ones) give independent numbers.
/// This gives reproducible numbers without storing a state, e.g. one stream per step
/// with restart(step), as GROMACS does for replica exchange.
class Threefry2x64 {
public:
  using Block = std::array<std::uint64_t,2>;
/// The Threefry2x64-20 function: one block of 128 random bits
  static Block block(const Block & counter, const Block & key);
  explicit Threefry2x64(std::uint64_t key0=0, std::uint64_t key1=0);
/// Restart the stream at the counters (c0,0), (c0,1), ...
  void restart(std::uint64_t c0);
/// Next 64 random bits of the stream
  std::uint64_t next();
/// Uniform number in the open interval (0,1), 53-bit resolution (as Random::U01)
  double U01();
private:
  Block key;
  Block counter{};
  Block buffer{};
  unsigned used = 2;
};

}

#endif
