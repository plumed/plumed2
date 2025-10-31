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

#include <string>
#include <vector>
#include <iosfwd>

namespace PLMD {

/// \ingroup TOOLBOX
class Random {
  static constexpr int IA=16807;
  static constexpr int IM=2147483647;
  static constexpr int IQ=127773;
  static constexpr int IR=2836;
  static constexpr int NTAB=32;
  static constexpr int NDIV=(1+(IM-1)/NTAB);
  static constexpr double fact=5.9604644775390625e-8;     /* 1 / 2^24  */
  static constexpr double EPS=3.0e-16;
  static constexpr double AM=1.0/IM;
  static constexpr double RNMX=1.0-EPS;
  static const std::string noname;
  bool incPrec;
  bool switchGaussian;
  double saveGaussian;
  int iy;
  int iv[NTAB];
  int idum;
  std::string name;
public:
  explicit Random(const std::string & title=noname);
  void setSeed(int idum);
  double RandU01();
  double U01();
  double U01d();
  int RandInt(int i);
  void Shuffle(std::vector<unsigned>& vec);
  void WriteStateFull(std::ostream &)const;
  void ReadStateFull (std::istream &);
  void fromString(const std::string & str);
  void toString(std::string & str)const;
  friend std::ostream & operator<<(std::ostream & out,const Random & rng) {
    rng.WriteStateFull(out);
    return out;
  }
  friend std::istream & operator>>(std::istream & in,Random & rng) {
    rng.ReadStateFull(in);
    return in;
  }
  double Gaussian();
  void IncreasedPrecis(bool i) {
    incPrec=i;
  }
};

}

#endif

