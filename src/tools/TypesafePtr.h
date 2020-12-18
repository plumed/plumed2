/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2018 The plumed team
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
#ifndef __PLUMED_tools_TypesafePtr_h
#define __PLUMED_tools_TypesafePtr_h

#include "Exception.h"

#include <memory>
#include <iosfwd>
#include <map>
#include <utility>
#include <mutex>
#include <cstdio>

namespace PLMD {

template<class T>
std::size_t typesafePtrSizeof() {
  return sizeof(T);
}

template<>
inline
std::size_t typesafePtrSizeof<void>() {
  return 0;
}

/**
\ingroup TOOLBOX
Class to deal with propoagation of typesafe pointers.

*/
class TypesafePtr {
  TypesafePtr(void* ptr, std::size_t nelem, unsigned long int flags):
    ptr(ptr),
    nelem(nelem),
    flags(flags)
  {}

public:

  static TypesafePtr fromSafePtr(void* safe);
  static constexpr unsigned short is_integral=3;
  static constexpr unsigned short is_floating_point=4;
  static constexpr unsigned short is_file=5;

  TypesafePtr() {}

  TypesafePtr(const void*ptr) :
    ptr(const_cast<void*>(ptr))
  {}

  ~TypesafePtr() {
  }


  TypesafePtr(const TypesafePtr&other) = delete;

  TypesafePtr & operator=(const TypesafePtr & other) = delete;

  TypesafePtr(TypesafePtr&&other):
    ptr(other.ptr),
    nelem(other.nelem),
    flags(other.flags)
  {
    other.ptr=nullptr;
  }

  TypesafePtr copy() const;

  TypesafePtr & operator=(TypesafePtr && other) {
    ptr=other.ptr;
    flags=other.flags;
    nelem=other.nelem;
    other.ptr=nullptr;
    return *this;
  }

  std::string type_str() const {
    auto type=(flags>>16)&0xff;
    if(type==0) return "wildcard";
    if(type==1) return "void";
    if(type==2) return "integral";
    if(type==3) return "integral";
    if(type==4) return "floating point";
    if(type==5) return "FILE";
    return "unknown";
  }

  template<typename T>
  T* get_priv(std::size_t nelem, bool byvalue) const {
    typedef typename std::remove_const<T>::type T_noconst;
    typedef typename std::remove_pointer<T>::type T_noptr;
    if(flags==0) return (T*) ptr; // no check
    auto size=flags&0xffff;
    auto type=(flags>>16)&0xff;
    // auto unsi=(flags>>24)&0x1; // ignored
    auto cons=(flags>>25)&0x7;

    // type=0: ignore check
    // type=1: void -> ignore check
    // type>5: undefined yet
    if(type!=0 && type!=1 && type<=5) {
      if(std::is_integral<T_noptr>::value && (type!=is_integral)) {
        throw ExceptionTypeError() <<"This command expects an integer type. Received a " << type_str() << " instead";
      }
      if(std::is_floating_point<T_noptr>::value && (type!=is_floating_point)) {
        throw ExceptionTypeError() <<"This command expects a floating point type. Received a " << type_str() << " instead";
      }
      if(std::is_same<std::remove_const<FILE>,T_noconst>::value && (type!=is_file)) {
        throw ExceptionTypeError() <<"This command expects a FILE. Received a " << type_str() << " instead";
      }
    }

    if(size>0 && typesafePtrSizeof<T_noptr>() >0 && typesafePtrSizeof<T_noptr>()!=size) {
      throw ExceptionTypeError() << "This command expects a type with size " << typesafePtrSizeof<T_noptr>() << ". Received type has size " << size << " instead";
    }

    if(!byvalue) if(cons==1) {
        throw ExceptionTypeError() << "This command is trying to take the address of an argument that was passed by value";
      }

    // cons==1 (by value) is here treated as cons==3 (const type*)
    if(!std::is_pointer<T>::value) {
      if(cons!=1 && cons!=2 && cons!=3) {
        throw ExceptionTypeError() << "This command expects a pointer or an value. It received a pointer-to-pointer instead";
      }
      if(!std::is_const<T>::value) {
        if(cons==3) {
          throw ExceptionTypeError() << "This command expects a modifiable pointer (T*). It received a non modifiable pointer instead (const T*)";
        } else if(cons==1) {
          throw ExceptionTypeError() << "This command expects a modifiable pointer (T*). It received a value instead (T)";
        }
      }
    } else {
      if(!std::is_const<T>::value) {
        if(cons==1) throw ExceptionTypeError() << "This command expects a pointer-to-pointer. It received a value intead";
        if(cons==2 || cons==3) throw ExceptionTypeError() << "This command expects a pointer-to-pointer. It received a pointer intead";
        if(!std::is_const<T_noptr>::value) {
          if(cons!=4) throw ExceptionTypeError() << "This command expects a modifiable pointer-to-pointer (T**)";
        } else {
          if(cons!=6) throw ExceptionTypeError() << "This command expects a modifiable pointer to unmodifiable pointer (const T**)";
        }
      } else {
        if(!std::is_const<T_noptr>::value) {
          if(cons!=4 && cons!=5) throw ExceptionTypeError() << "This command expects T*const* pointer, and can only receive T**  or T*const* pointers";
        }
      }
    }
    if(nelem>0 && this->nelem>0) if(!(nelem<=this->nelem)) {
        throw ExceptionTypeError() << "This command wants to access " << nelem << " from this pointer, but only " << this->nelem << " have been passed";
      }
    return (T*) ptr;
  }

  template<typename T>
  T* get(std::size_t nelem=0) const {
    return get_priv<T>(nelem,false);
  }

  template<typename T>
  T getVal() const {
    return *get_priv<const T>(1,true);
  }

  operator bool() const noexcept {
    return ptr;
  }

private:
  void* ptr=nullptr;
  std::size_t nelem=0;
  unsigned long int flags=0;

};

}


#endif
