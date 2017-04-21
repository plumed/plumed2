/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2017 The plumed team
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
#include "ActionSet.h"

using namespace std;
namespace PLMD {

ActionSet::ActionSet(PlumedMain&p):
  plumed(p) {
  (void) plumed; // to suppress warning about "unused plumed"
}

ActionSet::~ActionSet()
{
  for(int i=size()-1; i>=0; i--) delete (*this)[i];
}

void ActionSet::clearDelete() {
  for(int i=size()-1; i>=0; i--) delete (*this)[i];
  clear();
}


std::string ActionSet::getLabelList() const {
  std::string outlist;
  for(const auto & p : (*this)) {
    outlist+=dynamic_cast<Action*>(p)->getLabel()+" ";
  };
  return  outlist;
}

std::vector<std::string> ActionSet::getLabelVector() const {
  std::vector<std::string> outlist;
  for(const auto & p : (*this)) {
    outlist.push_back(dynamic_cast<Action*>(p)->getLabel());
  };
  return  outlist;
}




}
