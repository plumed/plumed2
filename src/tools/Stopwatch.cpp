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

#include "Stopwatch.h"
#include "Exception.h"
#include "Log.h"

#include <cstdio>
#include <iostream>
#include <vector>
#include <algorithm>

namespace PLMD {

// this is needed for friend operators
std::ostream& operator<<(std::ostream&os,const Stopwatch&sw) {
  return sw.log(os);
}

Stopwatch::~Stopwatch() {
  if(mylog && mylog->isOpen()) {
// Make sure paused watches are stopped.
// this is necessary e.g. to make sure the main watch present in PlumedMain
// is stopped correctly.
    for(auto & w : watches) {
      if(w.second.state==Watch::State::paused) w.second.start().stop();
    }
    *mylog << *this;
  }
}

std::ostream& Stopwatch::log(std::ostream&os)const {
  char buffer[1000];
  buffer[0]=0;
  for(unsigned i=0; i<40; i++) os<<" ";
  os<<"      Cycles        Total      Average      Minimum      Maximum\n";

  std::vector<std::string> names;
  for(const auto & it : watches) names.push_back(it.first);
  std::sort(names.begin(),names.end());

  const double frac=1.0/1000000000.0;

  for(const auto & name : names) {
    const Watch&t(watches.find(name)->second);
    os<<name;
    for(unsigned i=name.length(); i<40; i++) os<<" ";
    std::sprintf(buffer,"%12u %12.6f %12.6f %12.6f %12.6f\n", t.cycles, frac*t.total, frac*t.total/t.cycles, frac*t.min,frac*t.max);
    os<<buffer;
  }
  return os;
}

}




