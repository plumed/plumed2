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
#include "RegisterBase.h"

#include <mutex>
#include "tools/Tools.h"

namespace PLMD {

namespace {
class Singleton {
public:
  /// Mutex to avoid simultaneous registrations from multiple threads
  std::mutex registeringMutex;

  /// Count simultaneous registrations
  /// This is here if in the future we want to expand to recursive registrations
  unsigned registeringCounter=0;

  /// Take care of all exisiting registers
  std::vector<Register*> registers;

  /// Full path of the registering library
  std::string fullPath;
};

Singleton & getSingleton() {
  static Singleton singleton;
  return singleton;
}

}

Register::RegistrationLock::RegistrationLock(const std::string & fullPath):
  active(true) {
  pushDLRegistration(fullPath);
}

Register::RegistrationLock::~RegistrationLock() noexcept {
  if(active) {
    popDLRegistration();
  }
}

Register::RegistrationLock::RegistrationLock(RegistrationLock&& other) noexcept:
  active(other.active) {
  other.active=false;
}

Register::RegistrationLock Register::registrationLock(const std::string & fullPath) {
  return RegistrationLock(fullPath);
}

void Register::pushDLRegistration(const std::string & fullPath) {
  auto & singleton=getSingleton();
  singleton.registeringMutex.lock();
  singleton.fullPath=fullPath;
  if(singleton.registeringCounter>0) {
    singleton.registeringMutex.unlock();
    plumed_error()<<"recursive registrations are technically possible but disabled at this stage "<<singleton.registeringCounter;
  }
  singleton.registeringCounter++;
}

void Register::popDLRegistration() noexcept {
  auto & singleton=getSingleton();
  for(auto & reg : singleton.registers) {
    reg->clearStaged();
  }
  singleton.registeringCounter--;
  singleton.registeringMutex.unlock();
}

void Register::completeAllRegistrations(void* image) {
  auto & singleton=getSingleton();
  for(auto & reg : singleton.registers) {
    reg->completeRegistration(image);
  }
}

std::string Register::imageToString(void* image) {
  std::stringstream ss;
  ss << image;
  return ss.str();
}

bool Register::isDLRegistering() noexcept {
  auto & singleton=getSingleton();
  return singleton.registeringCounter>0;
}

const std::string Register::getRegisteringFullPath() noexcept {
  auto & singleton=getSingleton();
  return singleton.fullPath;
}

Register::Register() {
  auto & singleton=getSingleton();
  // this is to protect insertion
  std::unique_lock lock(singleton.registeringMutex);
  singleton.registers.push_back(this);
}

Register::~Register() noexcept {
  auto & singleton=getSingleton();
  // this is to protect removal
  std::unique_lock lock(singleton.registeringMutex);
  auto it=std::find(singleton.registers.begin(),singleton.registers.end(),this);
  if(it!=singleton.registers.end()) {
    singleton.registers.erase(it);
  }
}

std::vector<std::string> Register::getKeysWithDLHandle(void* image) const {
  std::vector<std::string> res;
  const auto prefix=imageToString(image)+":";
  for(auto & k : getKeys()) {
    if(Tools::startWith(k,prefix)) {
      if(!std::getenv("PLUMED_LOAD_ACTION_DEBUG")) {
        k=k.substr(prefix.length());
      }
      res.push_back(k);
    }
  }
  return res;
}

std::ostream & operator<<(std::ostream &log,const Register &reg) {
  std::vector<std::string> s(reg.getKeys());
  for(unsigned i=0; i<s.size(); i++) {
    log<<"  "<<s[i]<<"\n";
  }
  return log;
}


}
