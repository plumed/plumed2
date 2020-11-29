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

class TypesafePtrPool {
  std::mutex mtx;
  std::map<const void*,int> refcount;
public:
  void add(const void*ptr);
  void remove(const void*ptr);
  bool check(const void*ptr) {
    return refcount.find(ptr)!=refcount.end();
  }
  int useCount(const void*ptr) {
    auto f=refcount.find(ptr);
    if(f!=refcount.end()) return f->second;
    else return 0;
  }
  void forget(const void*ptr);
  void print(std::ostream & os);
};

/**
\ingroup TOOLBOX
Class to deal with propoagation of typesafe pointers.

*/
class TypesafePtr {
public:
  static constexpr unsigned short is_integral=3;
  static constexpr unsigned short is_floating_point=4;
  static constexpr unsigned short is_file=5;

  TypesafePtr() {}

  TypesafePtr(const void*ptr) :
    ptr(const_cast<void*>(ptr))
  {}

  TypesafePtr(TypesafePtrPool*pool,const void*ptr,std::size_t nelem=0,unsigned long flags=0) :
    pool(pool),
    ptr(const_cast<void*>(ptr)),
    nelem(nelem),
    flags(flags)
  {
    if(pool && ptr) pool->add(ptr);
  }

  ~TypesafePtr() {
    if(pool && ptr) pool->remove(ptr);
  }

  TypesafePtr(const TypesafePtr&other):
    pool(other.pool),
    ptr(other.ptr),
    nelem(other.nelem),
    flags(other.flags)
  {
    if(pool && ptr) pool->add(ptr);
  }

  TypesafePtr(TypesafePtr&&other):
    pool(other.pool),
    ptr(other.ptr),
    nelem(other.nelem),
    flags(other.flags)
  {
    other.pool=nullptr;
    other.ptr=nullptr;
  }

  TypesafePtr & operator=(const TypesafePtr & other) {
    plumed_assert(((other.flags>>25)&0x7)!=1);
    if(this==&other) return *this;
    if(pool && ptr) pool->remove(ptr);
    pool=other.pool;
    ptr=other.ptr;
    flags=other.flags;
    nelem=other.nelem;
    if(pool && ptr) pool->add(ptr);
    return *this;
  }

  TypesafePtr & operator=(TypesafePtr && other) {
    if(pool && ptr) pool->remove(ptr);
    pool=other.pool;
    ptr=other.ptr;
    flags=other.flags;
    nelem=other.nelem;
    other.pool=nullptr;
    other.ptr=nullptr;
    return *this;
  }

  template<typename T>
  T* get_priv(std::size_t nelem, bool byvalue) const {
    typedef typename std::remove_const<T>::type T_noconst;
    typedef typename std::remove_pointer<T>::type T_noptr;
    if(flags==0) return (T*) ptr; // no check
    if(pool && ptr) if(!pool->check(ptr)) throw ExceptionTypeError();
    auto size=flags&0xffff;
    auto type=(flags>>16)&0xff;
    // auto unsi=(flags>>24)&0x1; // ignored
    auto cons=(flags>>25)&0x7;

    // type=0: ignore check
    // type>5: undefined yet
    if(type!=0 && type<=5) {
      if(std::is_integral<T_noptr>::value && (type!=is_integral)) throw ExceptionTypeError() <<"Expecting integer type";
      if(std::is_floating_point<T_noptr>::value && (type!=is_floating_point)) throw ExceptionTypeError() <<"Expecting floating point type";
      if(std::is_same<std::remove_const<FILE>,T_noconst>::value && (type!=is_file)) throw ExceptionTypeError() <<"Expecting FILE type";
      // if T is void*, then we do no check on type of TypesafePtr
    }

    if(size>0 && typesafePtrSizeof<T_noptr>()!=size) throw ExceptionTypeError() << "Incorrect sizeof";

    if(!byvalue) if(cons==1) throw ExceptionTypeError() << "Cannot take the address of an argument passed by value";

    // cons==1 (by value) is here treated as cons==3 (const type*)
    if(!std::is_pointer<T>::value) {
      if(cons!=1 && cons!=2 && cons!=3) throw ExceptionTypeError() << "Cannot convert non-pointer to pointer";
      if(!std::is_const<T>::value) {
        if(cons!=2) throw ExceptionTypeError() << "Cannot convert const T* to T*";
      }
    } else {
      if(!std::is_const<T>::value) {
        if(!std::is_const<T_noptr>::value) {
          if(cons!=4) throw ExceptionTypeError() << "Only T** can be passed to T**";
        } else {
          if(cons!=6) throw ExceptionTypeError() << "Only const T** can be passed to const T**";
        }
      } else {
        if(!std::is_const<T_noptr>::value) {
          if(cons!=4 && cons!=5) throw ExceptionTypeError() << "Only T** and T*const* can be passed to T*const*";
        } else {
          if(cons!=4 && cons!=5 && cons!=6 && cons!=7) throw ExceptionTypeError() << "Only pointer-to-pointer can be passed to const T*const*";
        }
      }
    }
    if(nelem>0 && this->nelem>0) if(!(nelem<=this->nelem)) throw ExceptionTypeError() << "Incorrect number of elements";
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
  TypesafePtrPool*pool=nullptr;
  void* ptr=nullptr;
  std::size_t nelem=0;
  unsigned long int flags=0;
};

}


#endif
