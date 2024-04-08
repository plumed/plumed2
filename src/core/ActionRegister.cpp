/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2023 The plumed team
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
#include <sstream> //for std::stringstream 
#include <algorithm>

#ifdef __PLUMED_HAS_DLADDR
#include <dlfcn.h>
#endif

namespace PLMD {

namespace {
std::string imageToString(void* image) {
  std::stringstream ss;
  ss << image;
  return ss.str();
}
}


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
    if((*p).second.create==f) {
      m.erase(p); break;
    }
  }
}

void ActionRegister::add(std::string key,creator_pointer f,keywords_pointer k) {
  // this force each action to be registered as an uppercase string
  if ( std::any_of( std::begin( key ), std::end( key ), []( char c ) { return ( std::islower( c ) ); } ) ) plumed_error() << "Action: " + key + " cannot be registered, use only UPPERCASE characters";

  if(registeringCounter) {
    plumed_assert(!staged_m.count(key)) << "cannot registed action twice with the same name "<< key<<"\n";
    // Store a pointer to the function that creates keywords
    // A pointer is stored and not the keywords because all
    // Vessels must be dynamically loaded before the actions.
    staged_m.insert(std::pair<std::string,Item>(key, {f,k}));
  } else {
    plumed_assert(!m.count(key)) << "cannot registed action twice with the same name "<< key<<"\n";
    m.insert(std::pair<std::string,Item>(key, {f,k}));
  }
}

bool ActionRegister::check(const std::string & key) {
  std::vector<void*> images; // empty vector
  return check(images,key);
}

bool ActionRegister::check(const std::vector<void*> & images,const std::string & key) {
  if(m.count(key)>0) return true;
  for(auto image : images) {
    std::string k=imageToString(image)+":"+key;
    if(m.count(k)>0) return true;
  }
  return false;
}

std::unique_ptr<Action> ActionRegister::create(const ActionOptions&ao) {
  std::vector<void*> images; // empty vector
  return create(images,ao);
}

std::unique_ptr<Action> ActionRegister::create(const std::vector<void*> & images,const ActionOptions&ao) {
  if(ao.line.size()<1)return NULL;

  // Create a copy of the manual locally. The manual is
  // then added to the ActionOptions. This allows us to
  // ensure during construction that all the keywords for
  // the action have been documented. In addition, we can
  // generate the documentation when the user makes an error
  // in the input.
  std::unique_ptr<Action> action;
  if( check(images,ao.line[0]) ) {
    std::string found_key=ao.line[0];
    for(auto image = images.rbegin(); image != images.rend(); ++image) {
      auto key=imageToString(*image) + ":" + ao.line[0];
      if(m.count(key)>0) {
        found_key=key;
        break;
      }
    }
    Keywords keys; m[found_key].keys(keys);
    ActionOptions nao( ao,keys );
    action=m[found_key].create(nao);
  }
  return action;
}

bool ActionRegister::getKeywords(const std::string& action, Keywords& keys) {
  if ( check(action) ) {  m[action].keys(keys); return true; }
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
    Keywords keys; m[action].keys(keys);
    keys.print_template(action, include_optional);
    return true;
  } else {
    return false;
  }
}

std::vector<std::string> ActionRegister::getActionNames() const {
  std::vector<std::string> s;
  for(const auto & it : m) s.push_back(it.first);
  std::sort(s.begin(),s.end());
  return s;
}

std::ostream & operator<<(std::ostream &log,const ActionRegister&ar) {
  std::vector<std::string> s(ar.getActionNames());
  for(unsigned i=0; i<s.size(); i++) log<<"  "<<s[i]<<"\n";
  return log;
}

void ActionRegister::pushDLRegistration() {

  registeringMutex.lock();
  if(registeringCounter>0) {
    registeringMutex.unlock();
    plumed_error()<<"recursive registrations are technically possible but disabled at this stage";
  }
  registeringCounter++;
}

void ActionRegister::popDLRegistration() noexcept {
  try {
    staged_m.clear();
  } catch(...) {
    // should never happen
    std::fprintf(stderr,"Unexpected error in popDLRegistration\n");
    std::terminate();
  }

  registeringCounter--;
  registeringMutex.unlock();
}

void ActionRegister::completeRegistration(void* handle) {
  for(auto & it : staged_m) {
    auto key=imageToString(handle) + ":" + it.first;
    plumed_assert(!m.count(key)) << "cannot registed action twice with the same name "<< key<<"\n";
    m.insert(std::pair<std::string,Item>(key,it.second));
  }
  staged_m.clear();
}

ActionRegister::RegistrationLock::RegistrationLock(ActionRegister* ar):
  ar(ar)
{
  ar->pushDLRegistration();
}
ActionRegister::RegistrationLock::RegistrationLock(RegistrationLock&& other) noexcept:
  ar(other.ar)
{
  other.ar=nullptr;
}
ActionRegister::RegistrationLock::~RegistrationLock() noexcept {
  if(ar) ar->popDLRegistration();
}

ActionRegister::RegistrationLock ActionRegister::registrationLock() {
  return RegistrationLock(this);
}


}
