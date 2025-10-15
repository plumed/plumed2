
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2025 The plumed team
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
#include "TokenizedLine.h"

#include <iterator>
#include <string>
#include <algorithm>

namespace PLMD {
TokenizedLine::line wholeLineReplicas(const std::string& myLine) {
//NOTE: as now (exactly like before this PR) a vector will not accept only the first element as replicas:
// - @replicas:{1,2,3},4,5 do not work (KEY={@replicas:{1,2,3},4,5})
// - @replicas:{1,4,5},{2,4,5},{3,4,5} works (KEY=@replicas:{{1,4,5},{2,4,5},{3,4,5}}), it should be equivalent to the line above
// - 1,@replicas:{2,3,4},5 works (KEY={1,@replicas{2,3,4},5})
// the problem will be solved by passing to TokenizedLine the line with the first parentheses not removed
//
  if(Tools::startWith(myLine,Tools::replicaToken)) {
    auto tr = Tools::getWords(
                myLine.substr(Tools::replicaToken.length()),
                "\t\n ,");
    tr.insert(tr.begin(),myLine);
    return tr;
  } else {
    return {myLine};
  }
}

inline auto render(const TokenizedLine::line& l) {
  return l[0];
}

template <typename IT>
auto mapCreator(IT k, IT const end) {
  TokenizedLine::mapType toret;
  for (; k!=end; ++k) {
    auto eqpos=k->find('=');
    // We do not want to modify the original line
    std::string key = k->substr(0,eqpos);
    std::transform(key.begin(),key.end(),key.begin(),::toupper);
    if(eqpos != std::string::npos) {
      toret[key]=wholeLineReplicas(k->substr(eqpos+1));
    } else {
      //is a flag
      //maybe giving it a special value to confirm that it is indeed a flag?
      toret[key].emplace_back("");
    }
  }
  return toret;
}

TokenizedLine::TokenizedLine(vectorIt begin, vectorIt end):
  tokens(mapCreator(begin, end)) {}

TokenizedLine::TokenizedLine(const_vectorIt begin,
                             const_vectorIt end):
  tokens(mapCreator(begin, end)) {}

TokenizedLine::TokenizedLine(const std::vector<std::string>& dataline):
  PLMD::TokenizedLine(dataline.begin(),dataline.end()) {}

std::string TokenizedLine::convertToString(bool alsoClear) {
  std::string output;
  for(auto p=tokens.begin(); p!=tokens.end(); ++p) {
    auto tmp = render(p->second);
    if( tmp.find(" " )!=std::string::npos ) {
      output += " " + p->first+ "={" + tmp + "}";
    } else {
      output += " "+ p->first+ "=" + tmp;
    }
  }
  if(alsoClear) {
    tokens.clear();
  }
  return output;
}

std::string TokenizedLine::keyList(std::string_view sep ) {
  std::string mylist="";
  int i=0;
  std::string separator = "";
  for(const auto & l:tokens) {
    mylist = mylist + separator + l.first;
    if(i==0) {
      separator = std::string(sep);
      ++i;
    }
  }
  return mylist;
}

std::size_t TokenizedLine::size() const {
  return tokens.size();
}

bool TokenizedLine::empty() const {
  return tokens.empty();
}

std::string TokenizedLine::getKeyword(std::string_view key) const {
  auto keyArg = tokens.find(key);
  if( keyArg != tokens.end()) {
    return std::string(key) +"="+ render(keyArg->second);
  }
  return "";
}

bool TokenizedLine::readAndRemoveFlag(std::string_view key) {
  auto keytext = tokens.find(key);
  if(keytext != tokens.end()) {
    tokens.erase(keytext);
    return true;
  }
  return false;
}

std::string TokenizedLine::getValue(std::string_view key, int rep) const {
  auto keyArg = tokens.find(key);
  rep = (rep<0) ? -1:rep;
  if( keyArg != tokens.end()) {
    if (keyArg->second.size() == 1) {
      return keyArg->second[0];
    } else {
      if (rep < static_cast<long long int>(keyArg->second.size())-1) {
        return keyArg->second[rep+1];
      }
    }
  }
  return "";
}
std::string_view TokenizedLine::replica(const line& l,int rep) {
  rep = (rep<0) ? -1:rep;
  if (l.size() == 1) {
    return l[0];
  } else {
    //TODO: check that this does not go rep> l.size()
    return l[rep+1];

  }
}
} // namespace PLMD

