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
#include "ExchangePatterns.h"
#include "tools/Random.h"

namespace PLMD {

ExchangePatterns::ExchangePatterns():
  PatternFlag(NONE),
  NumberOfReplicas(1)
{}

ExchangePatterns::~ExchangePatterns() {}

void ExchangePatterns::setNofR(const int nrepl) {
  NumberOfReplicas=nrepl;
}

void ExchangePatterns::setFlag(const int flag) {
  PatternFlag=flag;
}

void ExchangePatterns::getFlag(int &flag) {
  flag=PatternFlag;
}

void ExchangePatterns::setSeed(const int seed)
{
  random.setSeed(seed);
}

void ExchangePatterns::getList(const TypesafePtr & ind)
{
  auto iind=ind.get<int*>(NumberOfReplicas);
  switch(PatternFlag)
  {
  case RANDOM:
    for(int i=0; i<NumberOfReplicas; i++) {
      int stat=1;
      while(stat) {
        stat=0;
        iind[i] = (int) (random.U01()*NumberOfReplicas);
        for(int j=0; j<i; j++) if(iind[i]==iind[j]) stat=1;
      }
    }
    break;
  case NEIGHBOR:
    for(int i=0; i<NumberOfReplicas; i++) iind[i]=i;
    break;
  }
}

}
