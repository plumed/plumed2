/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2019 The plumed team
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
#include "VesselRegister.h"
#include "Vessel.h"
#include <iostream>

namespace PLMD {
namespace vesselbase {

VesselRegister::~VesselRegister() {
  if(m.size()>0) {
    std::string names="";
    for(const auto & p : m) names+=p.first+" ";
    std::cerr<<"WARNING: Vessel "+ names +" has not been properly unregistered. This might lead to memory leak!!\n";
  }
}

VesselRegister& vesselRegister() {
  static VesselRegister ans;
  return ans;
}

void VesselRegister::remove(creator_pointer f) {
  for(auto p=m.begin(); p!=m.end(); ++p) {
    if((*p).second==f) {
      m.erase(p); break;
    }
  }
}

void VesselRegister::add(std::string keyword,creator_pointer f,keyword_pointer k,keyword_pointer ik) {
  plumed_massert(m.count(keyword)==0,"keyword has already been registered");
  m.insert(std::pair<std::string,creator_pointer>(keyword,f));
  k( keywords );   // Store the keywords for all the things
  // Store a pointer to the function that creates keywords
  // A pointer is stored and not the keywords because all
  // Vessels must be dynamically loaded before the actions.
  mk.insert(std::pair<std::string,keyword_pointer>(keyword,ik));
}

bool VesselRegister::check(std::string key) {
  if( m.count(key)>0 ) return true;
  return false;
}

std::unique_ptr<Vessel> VesselRegister::create(std::string keyword, const VesselOptions&da) {
  std::unique_ptr<Vessel> df;
  if(check(keyword)) {
    Keywords keys; mk[keyword](keys);
    VesselOptions nda( da,keys );
    df=m[keyword](nda);
  }
  return df;
}

Keywords VesselRegister::getKeywords() {
  return keywords;
}

}
}
