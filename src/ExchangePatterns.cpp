/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

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

using namespace std;

namespace PLMD{

void ExchangePatterns::setFlag(const int flag){
  PatternFlag=flag;
}

void ExchangePatterns::getFlag(int &flag){
  flag=PatternFlag;
}

void ExchangePatterns::setSeed(int seed)
{
  random.setSeed(seed);
}

void ExchangePatterns::getList(int *ind, int nrepl)
{
  /* in principle here we can add a switch(patter) case in order to get a list of exchanges dependent on a specific pattern */
  switch(PatternFlag)
  {
    case RANDOM:
      for(int i=0;i<nrepl;i++) {
        int stat=1;
        while(stat) {
          stat=0;
          ind[i] = random.RandU01()*nrepl;
          for(int j=0;j<i;j++) if(ind[i]==ind[j]) stat=1;
        }
      }
      break;
    case NEIGHBOR:
      for(int i=0;i<nrepl;i++) ind[i]=i; 
      break; 
  }
}

}
