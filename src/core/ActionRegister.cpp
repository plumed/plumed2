/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2020 The plumed team
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
#include "ActionRegister.h"
#include "tools/Tools.h"
#include "Action.h"
#include <algorithm>
#include <iostream>

namespace PLMD {

ActionRegister::~ActionRegister() {
  if(m.size()>0) {
    std::string names="";
    for(const auto & p : m) names+=p.first+" ";
    std::cerr<<"WARNING: Directive "+ names +" has not been properly unregistered. This might lead to memory leak!!\n";
  }
}

ActionRegister& actionRegister() {
  static ActionRegister ans;
  return ans;
}

void ActionRegister::remove(creator_pointer f) {
  for(auto p=m.begin(); p!=m.end(); ++p) {
    if((*p).second==f) {
      m.erase(p); break;
    }
  }
}

void ActionRegister::add(std::string key,creator_pointer f,keywords_pointer k) {
  // this force each action to be registered as an uppercase string
  if ( std::any_of( std::begin( key ), std::end( key ), []( char c ) { return ( std::islower( c ) ); } ) ) plumed_error() << "Action: " + key + " cannot be registered, use only UPPERCASE characters";
  if(m.count(key)) {
    m.erase(key);
    disabled.insert(key);
  } else {
    m.insert(std::pair<std::string,creator_pointer>(key,f));
    // Store a pointer to the function that creates keywords
    // A pointer is stored and not the keywords because all
    // Vessels must be dynamically loaded before the actions.
    mk.insert(std::pair<std::string,keywords_pointer>(key,k));
  };
}

bool ActionRegister::check(std::string key) {
  if(m.count(key)>0 && mk.count(key)>0) return true;
  return false;
}

std::unique_ptr<Action> ActionRegister::create(const ActionOptions&ao) {
  if(ao.line.size()<1)return NULL;
  // Create a copy of the manual locally. The manual is
  // then added to the ActionOptions. This allows us to
  // ensure during construction that all the keywords for
  // the action have been documented. In addition, we can
  // generate the documentation when the user makes an error
  // in the input.
  std::unique_ptr<Action> action;
  if( check(ao.line[0]) ) {
    Keywords keys; mk[ao.line[0]](keys);
    ActionOptions nao( ao,keys );
    action=m[ao.line[0]](nao);
  }
  return action;
}

bool ActionRegister::getKeywords(const std::string& action, Keywords& keys) {
  if ( check(action) ) {  mk[action](keys); return true; }
  return false;
}

bool ActionRegister::printManual(const std::string& action, const bool& vimout, const bool& spellout) {
  if ( check(action) ) {
    Keywords keys; getKeywords( action, keys );
    if( vimout ) {
      printf("%s",action.c_str()); keys.print_vim(); printf("\n");
    } else if( spellout ) {
      keys.print_spelling();
    } else {
      keys.print_html();
    }
    return true;
  } else {
    return false;
  }
}

bool ActionRegister::printTemplate(const std::string& action, bool include_optional) {
  if( check(action) ) {
    Keywords keys; mk[action](keys);
    keys.print_template(action, include_optional);
    return true;
  } else {
    return false;
  }
}

std::ostream & operator<<(std::ostream &log,const ActionRegister&ar) {
  std::vector<std::string> s;
  for(const auto & it : ar.m) s.push_back(it.first);
  std::sort(s.begin(),s.end());
  for(unsigned i=0; i<s.size(); i++) log<<"  "<<s[i]<<"\n";
  if(!ar.disabled.empty()) {
    s.assign(ar.disabled.size(),"");
    std::copy(ar.disabled.begin(),ar.disabled.end(),s.begin());
    std::sort(s.begin(),s.end());
    log<<"+++++++ WARNING +++++++\n";
    log<<"The following keywords have been registered more than once and will be disabled:\n";
    for(unsigned i=0; i<s.size(); i++) log<<"  - "<<s[i]<<"\n";
    log<<"+++++++ END WARNING +++++++\n";
  };
  return log;
}


}
