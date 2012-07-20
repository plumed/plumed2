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
#include "ActionRegister.h"
#include "Tools.h"
#include "Action.h"
#include <algorithm>
#include <iostream>


using namespace std;
using namespace PLMD;

ActionRegister::~ActionRegister(){
  if(m.size()>0){
    string names="";
    for(mIterator p=m.begin();p!=m.end();++p)names+=p->first+" ";
    plumed_merror("Directive "+ names +" has not been properly unregistered");
  }
}

ActionRegister& PLMD::actionRegister(){
  static ActionRegister ans;
  return ans;
}

void ActionRegister::remove(creator_pointer f){
  for(mIterator p=m.begin();p!=m.end();++p){
    if((*p).second==f){
      m.erase(p); break;
    }
  }
}

void ActionRegister::add(string key,creator_pointer f,keywords_pointer k){
  if(m.count(key)){
    m.erase(key);
    disabled.insert(key);
  }else{
    m.insert(pair<string,creator_pointer>(key,f));
    // Create keywords using a function pointer
    Keywords kk; k(kk);
    // Store array of keywords inside an associative array 
    // we can then find them using the keyword
    mk.insert(pair<string,Keywords>(key,kk));
  };
}

bool ActionRegister::check(string key){
  if(m.count(key)>0 && mk.count(key)>0) return true;
  return false;
}

Action* ActionRegister::create(const ActionOptions&ao){
  if(ao.line.size()<1)return NULL;
  Action* action; ActionOptions nao( ao,mk[ao.line[0]] );
  if(check(ao.line[0])) action=m[ao.line[0]](nao);
  else action=NULL;
  return action;
}

bool ActionRegister::printManual( const std::string& action ){
  if ( check(action) ){
     mk[action].print_html();
     return true;
  } else {
     return false;
  } 
}

std::ostream & PLMD::operator<<(std::ostream &log,const ActionRegister&ar){
  vector<string> s;
  for(ActionRegister::const_mIterator it=ar.m.begin();it!=ar.m.end();++it)
    s.push_back((*it).first);
  sort(s.begin(),s.end());
  for(unsigned i=0;i<s.size();i++) log<<"  "<<s[i]<<"\n";
  if(ar.disabled.size()>0){
    s.assign(ar.disabled.size(),"");
    copy(ar.disabled.begin(),ar.disabled.end(),s.begin());
    sort(s.begin(),s.end());
    log<<"+++++++ WARNING +++++++\n";
    log<<"The following keywords have been registered more than once and will be disabled:\n";
    for(unsigned i=0;i<s.size();i++) log<<"  - "<<s[i]<<"\n";
    log<<"+++++++ END WARNING +++++++\n";
  };
  return log;
}


