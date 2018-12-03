/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2018 The plumed team
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
#include "TypesafePtr.h"
#include <iostream>

namespace PLMD {


void TypesafePtrPool::add(const void*ptr) {
// lock:
//    std::lock_guard<std::mutex> lock(mtx);
    refcount[ptr]++;
  }
  void TypesafePtrPool::remove(const void*ptr) {
// lock:
//    std::lock_guard<std::mutex> lock(mtx);
    auto f=refcount.find(ptr);
    if(f!=refcount.end()) {
      f->second--;
      if(f->second<=0) refcount.erase(f);
    }
  }

  void TypesafePtrPool::print(std::ostream & os) {
    for(const auto & f : refcount) std::cout<<f.first<<" "<<f.second<<"\n";
  }

void TypesafePtrPool::forget(const void*ptr) {
    auto f=refcount.find(ptr);
    if(f!=refcount.end()) {
      refcount.erase(f);
    }
}

}


