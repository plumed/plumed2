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
#include "Citations.h"
#include "Exception.h"
#include "Tools.h"
#include <iostream>

using namespace std;
namespace PLMD {

std::string Citations::cite(const std::string & item) {
  unsigned i;
  for(i=0; i<items.size(); ++i) if(items[i]==item) break;
  if(i==items.size()) items.push_back(item);
  plumed_assert(i<items.size());
  string ret;
  Tools::convert(i+1,ret);
  ret="["+ret+"]";
  return ret;
}

std::ostream & operator<<(std::ostream &log,const Citations&cit) {
  for(unsigned i=0; i<cit.items.size(); ++i)
    log<<"  ["<<i+1<<"] "<<cit.items[i]<<"\n";
  return log;
}

void Citations::clear() {
  items.clear();
}

bool Citations::empty()const {
  return items.empty();
}

}


