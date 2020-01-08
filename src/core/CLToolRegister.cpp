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
#include "CLToolRegister.h"
#include "tools/Tools.h"
#include "CLTool.h"
#include <algorithm>
#include <iostream>


using namespace std;
namespace PLMD {

CLToolRegister::~CLToolRegister() {
  if(m.size()>0) {
    string names="";
    for(const auto & p : m) names+=p.first+" ";
    std::cerr<<"WARNING: CLTools "+ names +" has not been properly unregistered. This might lead to memory leak!!\n";
  }
}

CLToolRegister& cltoolRegister() {
  static CLToolRegister ans;
  return ans;
}

void CLToolRegister::remove(creator_pointer f) {
  for(auto p=m.begin(); p!=m.end(); ++p) {
    if((*p).second==f) {
      m.erase(p); break;
    }
  }
}

void CLToolRegister::add(string key,creator_pointer f,keywords_pointer kf) {
  if(m.count(key)) {
    m.erase(key);
    disabled.insert(key);
  } else {
    m.insert(pair<string,creator_pointer>(key,f));
    Keywords keys; kf(keys);
    mk.insert(pair<string,Keywords>(key,keys));
  };
}

bool CLToolRegister::check(string key)const {
  if(m.count(key)>0) return true;
  return false;
}

std::unique_ptr<CLTool> CLToolRegister::create(const CLToolOptions&ao) {
  if(ao.line.size()<1)return NULL;
  std::unique_ptr<CLTool> cltool;
  if(check(ao.line[0])) {
    CLToolOptions nao( ao,mk[ao.line[0]] );
    cltool=m[ao.line[0]](nao);
  }
  return cltool;
}


std::ostream & operator<<(std::ostream &log,const CLToolRegister&ar) {
  vector<string> s(ar.list());
  for(unsigned i=0; i<s.size(); i++) log<<"  "<<s[i]<<"\n";
  if(!ar.disabled.empty()) {
    s.assign(ar.disabled.size(),"");
    copy(ar.disabled.begin(),ar.disabled.end(),s.begin());
    sort(s.begin(),s.end());
    log<<"+++++++ WARNING +++++++\n";
    log<<"The following keywords have been registered more than once and will be disabled:\n";
    for(unsigned i=0; i<s.size(); i++) log<<"  - "<<s[i]<<"\n";
    log<<"+++++++ END WARNING +++++++\n";
  };
  return log;
}

bool CLToolRegister::printManual( const std::string& cltool, const bool& spelling ) {
  if( spelling && check(cltool) ) {
    mk[cltool].print_spelling();
    return true;
  } else if ( check(cltool) ) {
    mk[cltool].print_html();
    return true;
  } else {
    return false;
  }
}

std::vector<std::string> CLToolRegister::getKeys(const std::string& cltool)const {
  if ( check(cltool) ) {
    return mk.find(cltool)->second.getKeys();
  } else {
    std::vector<std::string> empty;
    return empty;
  }
}


vector<string> CLToolRegister::list()const {
  vector<string> s;
  for(const auto & it : m) s.push_back(it.first);
  sort(s.begin(),s.end());
  return s;
}



}
