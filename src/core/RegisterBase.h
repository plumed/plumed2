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
#ifndef __PLUMED_core_RegisterBase_h
#define __PLUMED_core_RegisterBase_h

#include "tools/Exception.h"
#include <string>
#include <map>
#include <memory>
#include <mutex>
#include <iostream>
#include <vector>
#include <algorithm>

namespace PLMD {

/// General register.
template<class Register,class Content>
class RegisterBase {
/// Main register map
  std::map<std::string,std::unique_ptr<Content>> m;
/// Map of staged keys
  std::map<std::string,std::unique_ptr<Content>> staged_m;
/// Mutex to avoid simultaneous registrations from multiple threads
/// It is a recursive mutex so that recursive calls will be detected and throw.
/// (a non recursive mutex would lead to a lock instead)
  std::recursive_mutex registeringMutex;
  unsigned registeringCounter=0;

  /// initiate registration
  /// all keys registered after this call will be staged
  /// Better use the RAII interface as registrationLock()
  void pushDLRegistration() {
    registeringMutex.lock();
    if(registeringCounter>0) {
      registeringMutex.unlock();
      plumed_error()<<"recursive registrations are technically possible but disabled at this stage";
    }
    registeringCounter++;
  }

  /// finish registration
  /// all keys that were staged will be removed.
  /// Better use the RAII interface as registrationLock()
  void popDLRegistration() noexcept {
    staged_m.clear();
    registeringCounter--;
    registeringMutex.unlock();
  }

  /// Internal tool
  static std::string imageToString(void* image) {
    std::stringstream ss;
    ss << image;
    return ss.str();
  }
public:
  struct ID {
    Content* ptr{nullptr};
  };
/// Register a new class.
/// \param key The name of the directive to be used in the input file
/// \param content The registered content
  ID add(std::string key,Content content) {

    auto ptr=std::make_unique<Content>(content);
    ID id{ptr.get()};
    if(registeringCounter) {
      plumed_assert(!staged_m.count(key)) << "cannot stage key twice with the same name "<< key<<"\n";
      staged_m.insert({key, std::move(ptr)});
    } else {
      plumed_assert(!m.count(key)) << "cannot register key twice with the same name "<< key<<"\n";
      m.insert({key, std::move(ptr)});
    }
    return id;
  }

/// Verify if a key is present in the register

  bool check(const std::vector<void*> & images,const std::string & key) const {
    if(m.count(key)>0) return true;
    for(auto image : images) {
      std::string k=imageToString(image)+":"+key;
      if(m.count(k)>0) return true;
    }
    return false;
  }

  bool check(const std::string & key) const {
    return m.count(key)>0;
  }

  const Content & get(const std::vector<void*> & images,const std::string & key) const {
    for(auto image = images.rbegin(); image != images.rend(); ++image) {
      auto qualified_key=imageToString(*image) + ":" + key;
      if(m.count(qualified_key)>0) return *(m.find(qualified_key)->second);
    }
    plumed_assert(m.count(key)>0);
    return *(m.find(key)->second);
  }

  const Content & get(const std::string & key) const {
    plumed_assert(m.count(key)>0);
    return *(m.find(key)->second);
  }


  void remove(ID id) {
    if(id.ptr) {
      for(auto p=m.begin(); p!=m.end(); ++p) {
        if(p->second.get()==id.ptr) {
          m.erase(p); break;
        }
      }
    }
  }
/// Get a list of keys
  std::vector<std::string> getKeys() const {
    std::vector<std::string> s;
    for(const auto & it : m) s.push_back(it.first);
    std::sort(s.begin(),s.end());
    return s;
  }

  ~RegisterBase() noexcept {
    if(m.size()>0) {
      std::string names="";
      for(const auto & p : m) names+=p.first+" ";
      std::cerr<<"WARNING: Directive "+ names +" has not been properly unregistered. This might lead to memory leak!!\n";
    }
  }
  /// complete registration
  /// all staged keys will be enabled
  /// Should be called after dlopen has been completed correctly.
  void completeRegistration(void*handle) {
    for (auto iter = staged_m.begin(); iter != staged_m.end(); ) {
      auto key = imageToString(handle) + ":" + iter->first;
      plumed_assert(!m.count(key)) << "cannot registed key twice with the same name "<< key<<"\n";
      m[key] = std::move(iter->second);
      // Since we've moved out the value, we can safely erase the element from the original map
      // This also avoids invalidating our iterator since erase returns the next iterator
      iter = staged_m.erase(iter);
    }
    plumed_assert(staged_m.empty());
  }

  /// small class to manage registration lock
  class RegistrationLock {
    RegisterBase* reg{nullptr};
  public:
    RegistrationLock(RegisterBase* reg) :
      reg(reg)
    {
      plumed_assert(reg);
      reg->pushDLRegistration();
    }
    RegistrationLock(const RegistrationLock&) = delete;
    RegistrationLock(RegistrationLock&& other) noexcept:
      reg(other.reg)
    {
      other.reg=nullptr;
    }
    ~RegistrationLock() noexcept {
      if(reg) reg->popDLRegistration();
    }
  };

  /// return a registration lock
  RegistrationLock registrationLock() {
    return RegistrationLock(this);
  }
};


template<class Register,class Content>
std::ostream & operator<<(std::ostream &log,const RegisterBase<Register,Content> &reg) {
  std::vector<std::string> s(reg.getKeys());
  for(unsigned i=0; i<s.size(); i++) log<<"  "<<s[i]<<"\n";
  return log;
}

}

#endif


