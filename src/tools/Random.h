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
#ifndef __PLUMED_tools_Random_h
#define __PLUMED_tools_Random_h

#include <cstdint>
#include <string>
#include <vector>
#include <iosfwd>
#include <utility>

namespace PLMD {

/// \ingroup TOOLBOX
/// Random number generator used everywhere in PLUMED.
///
/// The engine is xoshiro256** (Blackman and Vigna), 64-bit output, period 2^256-1,
/// seeded from a single integer with splitmix64. The conversions to uniform doubles,
/// Gaussian numbers, bounded integers and shuffles are implemented here (not with the
/// std:: distributions), so that a given seed gives the same stream on every compiler
/// and standard library.
class Random {
  static const std::string noname;
  std::uint64_t s[4];
  bool switchGaussian;
  double saveGaussian;
  std::string name;
/// Next 64 random bits
  std::uint64_t next();
public:
  explicit Random(const std::string & title=noname);
/// Seed the generator. Seeds that differ only in sign give the same stream.
  void setSeed(int idum);
/// Uniform number in the open interval (0,1), 53-bit resolution
  double RandU01() {
    return U01();
  }
  double U01();
/// Same as U01() (kept for compatibility: the engine already has full double resolution)
  double U01d() {
    return U01();
  }
/// Uniform integer in [0,n), without modulo bias. n must be positive.
  int RandInt(int n);
/// Uniform 64-bit integer in [0,n), without modulo bias. n must be positive.
  std::uint64_t RandUInt64(std::uint64_t n);
/// Gaussian number with zero mean and unit variance
  double Gaussian();
/// Uniform random permutation (Fisher-Yates)
  template<typename T>
  void Shuffle(std::vector<T>& vec) {
    for(std::size_t i=vec.size(); i>1; --i) {
      std::swap(vec[i-1],vec[RandUInt64(i)]);
    }
  }
  void WriteStateFull(std::ostream &)const;
  void ReadStateFull (std::istream &);
/// State of the uniform generator only (the cached Gaussian number is not included)
  void fromString(const std::string & str);
  void toString(std::string & str)const;
/// Complete state (including the cached Gaussian number), exact to the last bit: used by the checkpoint.
/// Unlike toString/fromString, which keep only the uniform generator.
  std::string stateToString()const;
  void stateFromString(const std::string & str);
  friend std::ostream & operator<<(std::ostream & out,const Random & rng) {
    rng.WriteStateFull(out);
    return out;
  }
  friend std::istream & operator>>(std::istream & in,Random & rng) {
    rng.ReadStateFull(in);
    return in;
  }
/// No effect, kept for compatibility (U01() already has 53-bit resolution)
  void IncreasedPrecis(bool) {}
};

}

#endif
