/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2020 The plumed team
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
#ifndef __PLUMED_core_ExchangePatterns_h
#define __PLUMED_core_ExchangePatterns_h

#include "tools/ForwardDecl.h"

namespace PLMD {
class Random;

class ExchangePatterns {
  int    PatternFlag;
  int    NumberOfReplicas;
  ForwardDecl<Random> random_fwd;
  Random& random=*random_fwd;
public:
  ExchangePatterns();
  ~ExchangePatterns();
  enum PatternFlags { NONE, RANDOM, NEIGHBOR, TOTAL };
  void setNofR(const int);
  void setSeed(const int);
  void setFlag(const int);
  void getList(int *ind);
  void getFlag(int&);
};
}
#endif
