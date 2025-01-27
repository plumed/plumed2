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
#include <string_view>
#include <map>
#include <memory>
#include <iostream>
#include <vector>
#include <algorithm>
#include <mutex>
#include <shared_mutex>

namespace PLMD {

/// Base class, with type independent information.
/// Actual registers should inherit through the RegisterBase class below
class Register {
  /// Initialize registration - only used by registrationLock()
  static void pushDLRegistration(const std::string & fullpath);
  /// Finalize registration - only used by registrationLock()
  static void popDLRegistration() noexcept;

protected:
  /// Mutex protecting access to map
  mutable std::shared_mutex mutex;
  /// Internal tool to format image addresses
  static std::string imageToString(void* image);
  /// Check if we are in a dlopen section
  static bool isDLRegistering() noexcept;
  /// Return the path of the currently-loading library
  static const std::string getRegisteringFullPath() noexcept;
  /// Save all staged objects from a register
  virtual void completeRegistration(void* image)=0;
  /// Clear staged objects.
  /// Should be used when leaving the dlopen section to remove
  /// any dangling object.
  virtual void clearStaged() noexcept =0;
  /// Get all registered keys.
  /// These are the keys in the map, not the plumed keywords!
  virtual std::vector<std::string> getKeys() const =0;

public:
  /// Constructor.
  /// This keeps track of all created instances.
  Register();
  /// Destructor.
  virtual ~Register() noexcept;
  /// Disable move
  Register(Register &&) = delete;
  /// Disable copy
  Register(const Register &) = delete;

  /// Small class to manage registration lock
  /// This is used during dlopen, to avoid data races in registrations
  class RegistrationLock {
    bool active;
  public:
    RegistrationLock(const std::string & fullpath);
    RegistrationLock(const RegistrationLock&) = delete;
    RegistrationLock(RegistrationLock&& other) noexcept;
    ~RegistrationLock() noexcept;
  };

  /// return a registration lock
  static RegistrationLock registrationLock(const std::string & fullpath);

  /// Save all staged objects in all registers
  static void completeAllRegistrations(void* image);

  /// Get only keys registered specifically by a given image
  std::vector<std::string> getKeysWithDLHandle(void* handle) const;

  friend std::ostream & operator<<(std::ostream &log,const Register &reg);
};

/// Class representing an error in a register
class ExceptionRegisterError :
  public Exception {
  /// The missing key
  std::string missingKey;
public:
  using Exception::Exception;
  /// Sets the missing key
  /// \param key The missing key
  /// \return This exception
  ///
  /// ExceptionRegisterError can be used as a builder pattern:
  /// `throw ExceptionRegisterError().setMissingKey(key);`
  ///
  /// the key can be retrieved with ExceptionRegisterError::getMissingKey()
  ExceptionRegisterError& setMissingKey (std::string_view key) {
    missingKey=key;
    return *this;
  }
  /// Returns the missing key
  const std::string& getMissingKey() const {
    return missingKey;
  }
  template<typename T>
  ExceptionRegisterError& operator<<(const T & x) {
    *static_cast<Exception*>(this) <<x;
    return *this;
  }
};

/// General register.
/// This class provide a generic implementation based on the content of the Register
template<class Content>
class RegisterBase :
  public Register {

public:
/// auxiliary class
  struct ContentAndFullPath {
    Content content;
    std::string fullPath;
  };

private:
/// Main register map
  std::map<std::string,std::unique_ptr<ContentAndFullPath>> m;
/// Map of staged keys
  std::map<std::string,std::unique_ptr<ContentAndFullPath>> staged_m;

public:

  struct ID {
    ContentAndFullPath* ptr{nullptr};
  };
/// Register a new class.
/// \param key The name of the directive to be used in the input file
/// \param content The registered content
/// \param ID A returned ID that can be used to remove the directive later
  ID add(std::string key,const Content & content);

/// Verify if a key is present in the register, accessing to registered images
  bool check(const std::vector<void*> & images,const std::string & key) const;

/// Verify if a key is present in the register, only considering the default image
  bool check(const std::string & key) const;

/// Return the content associated to a key in the register, accessing to registerd images
  const Content & get(const std::vector<void*> & images,const std::string & key) const;

/// Return the full path associated to a key in the register, accessing to registerd images
  const std::string & getFullPath(const std::vector<void*> & images,const std::string & key) const;

/// Return the content associated to a key in the register, only considering the default image
  const Content & get(const std::string & key) const;

/// Remove a registered keyword.
/// Use the ID returned by add().
  void remove(ID id);

/// Get a list of keys
/// Notice that these are the keys in the map, not the plumed keywords!
/// Also notice that this list includes keys from all images, including the
/// textual version of the image void*
  std::vector<std::string> getKeys() const override;

  ~RegisterBase() noexcept override;

  /// complete registration
  /// all staged keys will be enabled
  /// Should be called after dlopen has been completed correctly.
  void completeRegistration(void*handle) override;

  void clearStaged() noexcept override;

};

template<class Content>
typename RegisterBase<Content>::ID RegisterBase<Content>::add(std::string key,const Content & content) {

  auto ptr=std::make_unique<ContentAndFullPath>(ContentAndFullPath{content,getRegisteringFullPath()});
  ID id{ptr.get()};

  // lock map for writing
  std::unique_lock<std::shared_mutex> lock(mutex);

  if(isDLRegistering()) {
    plumed_assert(!staged_m.count(key)) << "cannot stage key twice with the same name "<< key<<"\n";
    staged_m.insert({key, std::move(ptr)});
  } else {
    plumed_assert(!m.count(key)) << "cannot register key twice with the same name "<< key<<"\n";
    m.insert({key, std::move(ptr)});
  }
  return id;
}
std::ostream & operator<<(std::ostream &log,const Register &reg);

template<class Content>
bool RegisterBase<Content>::check(const std::vector<void*> & images,const std::string & key) const {
  // lock map for reading
  std::shared_lock<std::shared_mutex> lock(mutex);
  if(m.count(key)>0) {
    return true;
  }
  for(auto image : images) {
    std::string k=imageToString(image)+":"+key;
    if(m.count(k)>0) {
      return true;
    }
  }
  return false;
}

template<class Content>
bool RegisterBase<Content>::check(const std::string & key) const {
  // lock map for reading
  std::shared_lock<std::shared_mutex> lock(mutex);
  return m.count(key)>0;
}

template<class Content>
const Content & RegisterBase<Content>::get(const std::vector<void*> & images,const std::string & key) const {
  // lock map for reading
  std::shared_lock<std::shared_mutex> lock(mutex);
  for(auto image = images.rbegin(); image != images.rend(); ++image) {
    auto qualified_key=imageToString(*image) + ":" + key;
    if(m.count(qualified_key)>0) {
      return m.find(qualified_key)->second->content;
    }
  }
  if (m.count(key) == 0 ) {
    throw ExceptionRegisterError().setMissingKey(key);
  }
  return m.find(key)->second->content;
}

template<class Content>
const std::string & RegisterBase<Content>::getFullPath(const std::vector<void*> & images,const std::string & key) const {
  // lock map for reading
  std::shared_lock<std::shared_mutex> lock(mutex);
  for(auto image = images.rbegin(); image != images.rend(); ++image) {
    auto qualified_key=imageToString(*image) + ":" + key;
    if(m.count(qualified_key)>0) {
      return m.find(qualified_key)->second->fullPath;
    }
  }
  if (m.count(key) == 0 ) {
    throw ExceptionRegisterError().setMissingKey(key);
  }
  return m.find(key)->second->fullPath;
}

template<class Content>
const Content & RegisterBase<Content>::get(const std::string & key) const {
  // lock map for reading
  std::shared_lock<std::shared_mutex> lock(mutex);
  if (m.count(key) == 0 ) {
    throw ExceptionRegisterError().setMissingKey(key);
  }
  return m.find(key)->second->content;
}

template<class Content>
void RegisterBase<Content>::remove(ID id) {
  // lock map for writing
  std::unique_lock<std::shared_mutex> lock(mutex);
  if(id.ptr) {
    for(auto p=m.begin(); p!=m.end(); ++p) {
      if(p->second.get()==id.ptr) {
        m.erase(p);
        break;
      }
    }
  }
}

template<class Content>
std::vector<std::string> RegisterBase<Content>::getKeys() const {
  // lock map for reading
  std::shared_lock<std::shared_mutex> lock(mutex);
  std::vector<std::string> s;
  for(const auto & it : m) {
    s.push_back(it.first);
  }
  std::sort(s.begin(),s.end());
  return s;
}

template<class Content>
RegisterBase<Content>::~RegisterBase() noexcept {
  if(m.size()>0) {
    std::string names="";
    for(const auto & p : m) {
      names+=p.first+" ";
    }
    std::cerr<<"WARNING: Directive "+ names +" has not been properly unregistered. This might lead to memory leak!!\n";
  }
}

template<class Content>
void RegisterBase<Content>::completeRegistration(void*handle) {
  // lock map for writing
  std::unique_lock<std::shared_mutex> lock(mutex);
  for (auto iter = staged_m.begin(); iter != staged_m.end(); ) {
    auto key = imageToString(handle) + ":" + iter->first;
    plumed_assert(!m.count(key)) << "cannot register key twice with the same name "<< key<<"\n";
    m[key] = std::move(iter->second);
    // Since we've moved out the value, we can safely erase the element from the original map
    // This also avoids invalidating our iterator since erase returns the next iterator
    iter = staged_m.erase(iter);
  }
  plumed_assert(staged_m.empty());
}

template<class Content>
void RegisterBase<Content>::clearStaged() noexcept {
  // lock map for writing
  std::unique_lock<std::shared_mutex> lock(mutex);
  staged_m.clear();
}
}

#endif
