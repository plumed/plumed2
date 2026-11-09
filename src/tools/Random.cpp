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
#include "Random.h"
#include "Exception.h"
#include <cmath>
#include <sstream>
#include <limits>
#include <iostream>

namespace PLMD {

const std::string Random::noname="noname";

namespace {
/// tag written in front of the state, so that states of other engines are rejected
const std::string engineTag="xoshiro256**";

inline std::uint64_t rotl(const std::uint64_t x, int k) {
  return (x << k) | (x >> (64 - k));
}

/// splitmix64 (Steele, Lea and Flood): used to expand the seed into the engine state
inline std::uint64_t splitmix64(std::uint64_t & z) {
  z += UINT64_C(0x9e3779b97f4a7c15);
  std::uint64_t r = z;
  r = (r ^ (r >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
  r = (r ^ (r >> 27)) * UINT64_C(0x94d049bb133111eb);
  return r ^ (r >> 31);
}

void checkTag(const std::string & tag,const std::string & name) {
  if(tag!=engineTag) {
    plumed_error() << "the state of random number generator " << name
                   << " was not written by this version of PLUMED (expected " << engineTag
                   << ", found \"" << tag << "\"): the random number engine changed in PLUMED 2.11";
  }
}
}

Random::Random(const std::string & title):
  switchGaussian(false),
  saveGaussian(0.0),
// this is required because if a Random object is created during
// initialization than Random::noname could still be initialized.
// In practice: without this it is not possible to declare
// a static Random object without enforcing the order of the
// static constructors.
  name(&title!=&noname?title:"noname") {
  setSeed(0);
}

void Random::setSeed(int idum) {
  // the sign is irrelevant (several actions pass -seed): use the magnitude
  std::uint64_t z = idum<0 ? -static_cast<std::uint64_t>(static_cast<std::int64_t>(idum))
                    : static_cast<std::uint64_t>(idum);
  for(unsigned i=0; i<4; i++) {
    s[i]=splitmix64(z);
  }
  switchGaussian=false;
  saveGaussian=0.0;
}

std::uint64_t Random::next() {
  // xoshiro256** 1.0 (Blackman and Vigna, public domain reference implementation)
  const std::uint64_t result = rotl(s[1] * 5, 7) * 9;
  const std::uint64_t t = s[1] << 17;
  s[2] ^= s[0];
  s[3] ^= s[1];
  s[1] ^= s[2];
  s[0] ^= s[3];
  s[2] ^= t;
  s[3] = rotl(s[3], 45);
  return result;
}

double Random::U01() {
  // the 53 high bits, centered in their interval: the result is never 0 or 1
  return (static_cast<double>(next() >> 11) + 0.5) * 0x1.0p-53;
}

std::uint64_t Random::RandUInt64(std::uint64_t n) {
  plumed_assert(n>0) << "Random::RandUInt64 needs a positive range";
  // reject the lowest (2^64 mod n) values, so that every residue is equally likely
  const std::uint64_t threshold = (0 - n) % n;
  while(true) {
    const std::uint64_t r = next();
    if(r >= threshold) {
      return r % n;
    }
  }
}

int Random::RandInt(int n) {
  plumed_assert(n>0) << "Random::RandInt needs a positive range, got " << n;
  return static_cast<int>(RandUInt64(static_cast<std::uint64_t>(n)));
}

void Random::WriteStateFull(std::ostream & out)const {
  out<<name<<std::endl;
  out<<engineTag;
  for(unsigned i=0; i<4; i++) {
    out<<" "<<s[i];
  }
  out<<" "<<switchGaussian;
  out<<" "<<saveGaussian;
  out<<std::endl;
}

void Random::ReadStateFull (std::istream & in) {
  getline(in,name);
  std::string tag;
  in>>tag;
  checkTag(tag,name);
  for(unsigned i=0; i<4; i++) {
    in>>s[i];
  }
  in>>switchGaussian;
  in>>saveGaussian;
  if(in.fail()) {
    plumed_error() << "cannot read the state of random number generator " << name;
  }
}

void Random::toString(std::string & str)const {
  std::ostringstream ostr;
  ostr<<engineTag;
  for(unsigned i=0; i<4; i++) {
    ostr<<"|"<<s[i];
  }
  str=ostr.str();
}

void Random::fromString(const std::string & str) {
  std::string st=str;
  for(unsigned i=0; i<st.length(); i++)
    if(st[i]=='|') {
      st[i]=' ';
    }
  std::istringstream istr(st);
  std::string tag;
  istr>>tag;
  checkTag(tag,name);
  for(unsigned i=0; i<4; i++) {
    istr>>s[i];
  }
  if(istr.fail()) {
    plumed_error() << "cannot read the state of random number generator " << name;
  }
}

std::string Random::stateToString()const {
  std::ostringstream ostr;
  ostr.precision(std::numeric_limits<double>::max_digits10);
  WriteStateFull(ostr);
  return ostr.str();
}

void Random::stateFromString(const std::string & str) {
  std::istringstream istr(str);
  ReadStateFull(istr);
}

// This allows to have the same stream of random numbers
// with different compilers:
#ifdef __INTEL_COMPILER
#pragma intel optimization_level 0
#endif

double Random::Gaussian() {
  double v1,v2,rsq;
  if(switchGaussian) {
    switchGaussian=false;
    return saveGaussian;
  }
  while(true) {
    v1=2.0*U01()-1.0;
    v2=2.0*U01()-1.0;
    rsq=v1*v1+v2*v2;
    if(rsq<1.0 && rsq>0.0) {
      break;
    }
  }
  double fac=std::sqrt(-2.*std::log(rsq)/rsq);
  saveGaussian=v1*fac;
  switchGaussian=true;
  return v2*fac;
}

}
