/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2020 The plumed team
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
#include "MetricRegister.h"
#include <iostream>

namespace PLMD {

MetricRegister::~MetricRegister() {
  if(m.size()>0) {
    std::string names="";
    for(const auto & p : m) names+=p.first+" ";
    std::cerr<<"WARNING: ReferenceConfiguration "+ names +" has not been properly unregistered. This might lead to memory leak!!\n";
  }
}

MetricRegister& metricRegister() {
  static MetricRegister ans;
  return ans;
}

void MetricRegister::remove(creator_pointer f) {
  for(auto p=m.begin(); p!=m.end(); ++p) {
    if((*p).second==f) {
      m.erase(p); break;
    }
  }
}

void MetricRegister::add( std::string type, creator_pointer f ) {
  plumed_massert(m.count(type)==0,"type has already been registered");
  m.insert(std::pair<std::string,creator_pointer>(type,f));
}

bool MetricRegister::check(std::string type) {
  if( m.count(type)>0 ) return true;
  return false;
}

}
