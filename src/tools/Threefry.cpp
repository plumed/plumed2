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
#include "Threefry.h"

namespace PLMD {

namespace {
inline std::uint64_t rotl(std::uint64_t x, unsigned r) {
  return (x << r) | (x >> (64 - r));
}
}

// Same algorithm as threefry2x64 with 20 rounds in Random123 (threefry.h)
Threefry2x64::Block Threefry2x64::block(const Block & in, const Block & k) {
  constexpr unsigned rotations[8] = {16,42,12,31,16,32,24,21};
  const std::uint64_t ks[3] = {k[0], k[1], UINT64_C(0x1BD11BDAA9FC1A22)^k[0]^k[1]};
  std::uint64_t x0 = in[0] + ks[0];
  std::uint64_t x1 = in[1] + ks[1];
  for(unsigned r=0; r<20; r++) {
    x0 += x1;
    x1 = rotl(x1, rotations[r%8]);
    x1 ^= x0;
    // key injection every 4 rounds
    if(r%4 == 3) {
      const unsigned j = (r+1)/4;
      x0 += ks[j%3];
      x1 += ks[(j+1)%3] + j;
    }
  }
  return {x0, x1};
}

Threefry2x64::Threefry2x64(std::uint64_t key0, std::uint64_t key1):
  key{key0,key1} {
}

void Threefry2x64::restart(std::uint64_t c0) {
  counter = {c0,0};
  used = 2;
}

std::uint64_t Threefry2x64::next() {
  if (used == 2) {
    buffer = block(counter, key);
    counter[1]++;
    used = 0;
  }
  return buffer[used++];
}

double Threefry2x64::U01() {
  // the 53 high bits, centered in their interval: the result is never 0 or 1
  return (static_cast<double>(next() >> 11) + 0.5) * 0x1.0p-53;
}

}
